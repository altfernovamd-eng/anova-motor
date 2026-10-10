# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Fernando Ramírez Nova
# Referencia de R para los modelos mixtos logísticos (glmm.ts): proporciones (binomial-normal) y odds ratio (hipergeométrica-normal).
#
# Qué se guarda para cada conjunto (semilla fija; 40 de proporciones y 40 de tablas 2×2 simuladas, con eventos raros, estudios sin
# eventos o con todos los eventos, tamaños y heterogeneidad variados):
#  1. El óptimo de la verosimilitud marginal con una implementación INDEPENDIENTE en R (Gauss-Hermite adaptativa propia, `optim` BFGS con
#     tolerancia estricta, Hessiano de `numDeriv`): μ, ln τ, log-verosimilitud y error estándar de μ (inversa del Hessiano completo).
#  2. metafor::rma.glmm como CAJA NEGRA: `PLO` (binomial-normal; en el fondo lme4) y `OR` con `model = "CM.EL"` (hipergeométrica-normal),
#     solo si converge (CM.EL falla con frecuencia: es una limitación conocida de metafor, no de la validación).
#  3. La verosimilitud marginal por FUERZA BRUTA (`integrate` por tramos) en un conjunto pequeño de cada familia y 5 puntos del espacio de
#     parámetros; las probabilidades de la hipergeométrica no central salen de BiasedUrn::dFNCHypergeo (otro código).
#
#   Rscript validacion/generar-referencia-glmm.R
# Escribe validacion/referencias/referencia-glmm.json (lo lee validacion/comparar-glmm.mts).

lib <- file.path(Sys.getenv("HOME"), "Library", "R-anova")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
pdf(NULL)
suppressPackageStartupMessages({ library(metafor); library(jsonlite); library(numDeriv); library(BiasedUrn) })
set.seed(20261009)
options(digits = 17)

lsig <- function(x) -log1p(exp(-abs(x))) + pmin(x, 0)
gh <- function(q) { J <- matrix(0, q, q); for (i in 1:(q - 1)) J[i, i + 1] <- J[i + 1, i] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values); list(z = e$values[o], w = e$vectors[1, o]^2) }
nod <- gh(21)

# Familias: logP(eta), d1, d2
fam_bin <- function(x, n) list(logP = function(eta) lchoose(n, x) + x * lsig(eta) + (n - x) * lsig(-eta),
                               d1 = function(eta) x - n * plogis(eta), d2 = function(eta) { p <- plogis(eta); -n * p * (1 - p) })
fam_hyp <- function(a, n1, c, n2) { s <- a + c; u <- max(0, s - n2):min(n1, s); lw <- lchoose(n1, u) + lchoose(n2, s - u); base <- lchoose(n1, a) + lchoose(n2, s - a)
  mom <- function(eta) { t <- lw + eta * u; m <- max(t); w <- exp(t - m); z <- sum(w); mu <- sum(w * u) / z; list(lnZ = m + log(z), mu = mu, v = sum(w * u^2) / z - mu^2) }
  list(logP = function(eta) base + eta * a - mom(eta)$lnZ, d1 = function(eta) a - mom(eta)$mu, d2 = function(eta) -mom(eta)$v) }

ll_estudio <- function(f, mu, tau) {
  g <- function(z) f$logP(mu + tau * z) - 0.5 * z^2
  z <- 0
  for (it in 1:100) { d1 <- tau * f$d1(mu + tau * z) - z; d2 <- tau^2 * f$d2(mu + tau * z) - 1; st <- -d1 / d2; s <- 1; g0 <- g(z)
    while (s > 1e-10 && !(g(z + s * st) >= g0 - 1e-14)) s <- s / 2
    zn <- z + s * st; if (abs(zn - z) < 1e-13) { z <- zn; break }; z <- zn }
  sc <- 1 / sqrt(1 - tau^2 * f$d2(mu + tau * z))
  t <- sapply(seq_along(nod$z), function(k) g(z + sc * nod$z[k]) + 0.5 * nod$z[k]^2 + log(nod$w[k]))
  m <- max(t); log(sc) + m + log(sum(exp(t - m)))
}
LIMT <- c(log(1e-3), log(30))
loglik <- function(par, fams) { mu <- min(max(par[1], -30), 30); tau <- exp(min(max(par[2], LIMT[1]), LIMT[2])); sum(vapply(fams, ll_estudio, 0, mu = mu, tau = tau)) }
ajustar <- function(fams, mu0) {
  ob <- NULL
  for (lt in c(log(0.3), log(1), log(0.05))) {
    o <- tryCatch(optim(c(mu0, lt), function(p) -loglik(p, fams), method = "BFGS", control = list(reltol = 1e-14, maxit = 500)), error = function(e) NULL)
    if (!is.null(o) && (is.null(ob) || o$value < ob$value)) ob <- o
  }
  par <- c(min(max(ob$par[1], -30), 30), min(max(ob$par[2], LIMT[1]), LIMT[2]))
  limite <- par[2] <= LIMT[1] + 1e-3
  hs <- tryCatch(numDeriv::hessian(function(p) -loglik(p, fams), par), error = function(e) NULL)
  se <- NA
  if (!limite && !is.null(hs)) { iv <- tryCatch(solve(hs), error = function(e) NULL); if (!is.null(iv) && iv[1, 1] > 0) se <- sqrt(iv[1, 1]) }
  if (is.na(se)) { h1 <- tryCatch(numDeriv::hessian(function(m) -loglik(c(m, par[2]), fams), par[1]), error = function(e) NULL); if (!is.null(h1) && h1[1, 1] > 0) se <- sqrt(1 / h1[1, 1]) }
  list(mu = par[1], tau = exp(par[2]), loglik = -ob$value, se = se, limite = limite)
}

# ---------------- conjuntos de proporciones ----------------
conj <- list()
reusar <- nzchar(Sys.getenv("GLMM_REUSAR")) && file.exists("/tmp/glmm-conj.rds")
if (reusar) conj <- readRDS("/tmp/glmm-conj.rds")
for (i in if (reusar) integer(0) else 1:40) {
  k <- sample(c(3, 5, 8, 12, 20, 40), 1); mu <- sample(c(-4, -2.5, -1, 0, 1.2, 2.8), 1); tau <- sample(c(0.05, 0.3, 0.7, 1.3), 1)
  n <- if (i %% 3 == 0) round(runif(k, 8, 40)) else round(runif(k, 30, 400)); x <- rbinom(k, n, plogis(mu + rnorm(k, 0, tau)))
  fams <- lapply(seq_len(k), function(j) fam_bin(x[j], n[j]))
  cat(sprintf("prop%02d (k=%d)\n", i, k)); flush.console()
  r <- list(id = sprintf("prop%02d", i), datos = data.frame(x = x, n = n), indep = ajustar(fams, mean(log((x + 0.5) / (n - x + 0.5)))))
  mf <- tryCatch(suppressWarnings(rma.glmm(measure = "PLO", xi = x, ni = n)), error = function(e) NULL)
  if (!is.null(mf)) r$metafor <- list(mu = as.numeric(mf$beta), se = mf$se, tau2 = mf$tau2)
  conj[[length(conj) + 1]] <- r
}
# ---------------- conjuntos de tablas 2x2 ----------------
for (i in if (reusar) integer(0) else 1:40) {
  k <- sample(c(3, 5, 8, 12, 20, 35), 1); th <- sample(c(-1, -0.4, 0, 0.5, 1.2), 1); tau <- sample(c(0, 0.2, 0.6, 1.1), 1); b0 <- sample(c(-5, -3.5, -2, -0.8), 1)
  n1 <- if (i %% 3 == 0) round(runif(k, 10, 60)) else round(runif(k, 40, 400)); n2 <- pmax(8, round(n1 * runif(k, 0.6, 1.6)))
  a0 <- b0 + rnorm(k, 0, 0.6); thi <- th + rnorm(k, 0, tau)
  a <- rbinom(k, n1, plogis(a0 + thi)); c <- rbinom(k, n2, plogis(a0))
  fams_idx <- which((pmin(n1, a + c) - pmax(0, a + c - n2)) > 0)
  if (length(fams_idx) < 2) next
  fams <- lapply(fams_idx, function(j) fam_hyp(a[j], n1[j], c[j], n2[j]))
  cat(sprintf("or%02d (k=%d, informativos=%d)\n", i, k, length(fams_idx))); flush.console()
  mu0 <- mean(log(((a[fams_idx] + 0.5) * (n2[fams_idx] - c[fams_idx] + 0.5)) / ((n1[fams_idx] - a[fams_idx] + 0.5) * (c[fams_idx] + 0.5))))
  r <- list(id = sprintf("or%02d", i), datos = data.frame(a = a, n1 = n1, c = c, n2 = n2), k_informativos = length(fams_idx), indep = ajustar(fams, mu0))
  mf <- tryCatch(suppressWarnings(rma.glmm(measure = "OR", ai = a, n1i = n1, ci = c, n2i = n2, model = "CM.EL")), error = function(e) NULL)
  if (!is.null(mf)) r$metafor <- list(mu = as.numeric(mf$beta), se = mf$se, tau2 = mf$tau2)
  conj[[length(conj) + 1]] <- r
}

saveRDS(conj, file.path(tempdir(), "glmm-conj.rds")); saveRDS(conj, "/tmp/glmm-conj.rds")
# ---------------- fuerza bruta en conjuntos pequeños ----------------
integ_tramos <- function(f, lo, hi, n = 40) {
  fmax <- max(f(seq(lo, hi, length.out = 4001))); if (!(fmax > 0)) return(0)   # se normaliza el pico a 1 para fijar la tolerancia absoluta
  g <- function(z) f(z) / fmax
  cortes <- seq(lo, hi, length.out = n + 1)
  fmax * sum(vapply(1:n, function(i) tryCatch(integrate(g, cortes[i], cortes[i + 1], rel.tol = 1e-11, abs.tol = 1e-17, subdivisions = 300)$value,
                                              error = function(e) integrate(g, cortes[i], cortes[i + 1], rel.tol = 1e-9, abs.tol = 1e-15, subdivisions = 600)$value), 0))
}
brute <- function(par, pmf_list) { mu <- par[1]; tau <- exp(par[2])
  sum(vapply(pmf_list, function(pm) log(integ_tramos(function(z) vapply(z, function(zz) pm(mu + tau * zz) * dnorm(zz), 0), -12, 12)), 0)) }
xb <- c(3, 9, 0, 14, 5, 11); nb <- c(20, 31, 12, 40, 18, 25)
pm_bin <- lapply(seq_along(xb), function(j) function(eta) dbinom(xb[j], nb[j], plogis(eta)))
ab <- c(2, 0, 5, 1, 7, 3); n1b <- c(25, 18, 40, 22, 35, 30); cb <- c(6, 3, 4, 0, 12, 8); n2b <- c(28, 20, 38, 24, 33, 31)
pm_hyp <- lapply(seq_along(ab), function(j) function(eta) dFNCHypergeo(ab[j], n1b[j], n2b[j], ab[j] + cb[j], exp(eta)))
pts <- list(c(-1, log(0.4)), c(0.3, log(1.1)), c(-2.2, log(0.08)), c(0.9, log(2)), c(0, log(0.7)))
cat("fuerza bruta\n"); flush.console()
fams_bin_b <- lapply(seq_along(xb), function(j) fam_bin(xb[j], nb[j]))
fams_hyp_b <- lapply(seq_along(ab), function(j) fam_hyp(ab[j], n1b[j], cb[j], n2b[j]))
bruta <- list(
  proporcion = list(datos = data.frame(x = xb, n = nb), puntos = lapply(pts, function(p) list(par = p, bruta = brute(p, pm_bin), gh = loglik(p, fams_bin_b)))),
  odds = list(datos = data.frame(a = ab, n1 = n1b, c = cb, n2 = n2b), puntos = lapply(pts, function(p) list(par = p, bruta = brute(p, pm_hyp), gh = loglik(p, fams_hyp_b)))))

escribir <- toJSON(list(version = R.version.string, paquetes = list(metafor = as.character(packageVersion("metafor")), BiasedUrn = as.character(packageVersion("BiasedUrn"))), conjuntos = conj, bruta = bruta),
                   digits = NA, auto_unbox = TRUE, na = "null", null = "null")
write(escribir, file.path("validacion", "referencias", "referencia-glmm.json"))
cat("Escrito: validacion/referencias/referencia-glmm.json (", length(conj), "conjuntos )\n")
