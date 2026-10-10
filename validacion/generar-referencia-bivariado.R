# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Fernando Ramírez Nova
# Referencia de R para el modelo bivariado de precisión diagnóstica (sroc-bivariado.ts).
#
# Qué se guarda para cada conjunto (semilla fija; 30 simulados con heterogeneidad, correlación, tamaños y ceros variados + AuditC de `mada`):
#  1. mada::reitsma (aproximación normal, REML) como CAJA NEGRA: coeficientes, covarianza, Psi, IC de sensibilidad y FPR, AUC, curvas y elipse.
#  2. El óptimo EXACTO del REML normal, con una implementación propia en R de la misma función objetivo (optim, tolerancia estricta),
#     para medir la precisión del optimizador de aNova (el de mvmeta se detiene antes).
#  3. El óptimo de la verosimilitud BINOMIAL-NORMAL exacta con una implementación INDEPENDIENTE en R (cuadratura de Gauss-Hermite adaptativa
#     con nodos de `eigen()`, optimizador `optim` con tolerancia estricta, Hessiano de `numDeriv`): la verdad contra la que se compara aNova.
#  4. lme4::glmer (Laplace) como caja negra, solo para comprobar que el modelo está definido igual que el de la literatura.
#  5. La log-verosimilitud marginal calculada por FUERZA BRUTA (integrate anidado) en algunos puntos, para verificar la cuadratura.
#
#   Rscript validacion/generar-referencia-bivariado.R
# Escribe validacion/referencias/referencia-bivariado.json (lo lee validacion/comparar-bivariado.mts).

lib <- file.path(Sys.getenv("HOME"), "Library", "R-anova")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
pdf(NULL)
suppressPackageStartupMessages({ library(mada); library(lme4); library(jsonlite); library(numDeriv) })
set.seed(20261007)
options(digits = 17)

# ---------------- verosimilitud binomial-normal: implementación propia ----------------
gh <- function(q) { J <- matrix(0, q, q); for (i in 1:(q - 1)) J[i, i + 1] <- J[i + 1, i] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values); list(z = e$values[o], w = e$vectors[1, o]^2) }
lsig <- function(x) -log1p(exp(-abs(x))) + pmin(x, 0)
ll_estudio <- function(vp, fn, vn, fp, mu, S, nod) {
  P <- solve(S); ld <- log(det(S)); n1 <- vp + fn; n0 <- vn + fp
  g <- function(a, b) { d1 <- a - mu[1]; d2 <- b - mu[2]
    vp * lsig(a) + fn * lsig(-a) + vn * lsig(b) + fp * lsig(-b) - log(2 * pi) - 0.5 * ld - 0.5 * (P[1, 1] * d1^2 + 2 * P[1, 2] * d1 * d2 + P[2, 2] * d2^2) }
  th <- mu
  for (it in 1:100) {
    pa <- plogis(th[1]); pb <- plogis(th[2]); d <- th - mu
    gr <- c(vp - n1 * pa, vn - n0 * pb) - as.numeric(P %*% d)
    H <- -diag(c(n1 * pa * (1 - pa), n0 * pb * (1 - pb))) - P
    stp <- -solve(H, gr); s <- 1; g0 <- g(th[1], th[2])
    while (s > 1e-10 && !(g(th[1] + s * stp[1], th[2] + s * stp[2]) >= g0 - 1e-14)) s <- s / 2
    th2 <- th + s * stp; if (sum(abs(th2 - th)) < 1e-12) { th <- th2; break }; th <- th2
  }
  pa <- plogis(th[1]); pb <- plogis(th[2])
  H <- -diag(c(n1 * pa * (1 - pa), n0 * pb * (1 - pb))) - P
  L <- t(chol(solve(-H))); q <- length(nod$z); zz <- expand.grid(i = 1:q, j = 1:q)
  z1 <- nod$z[zz$i]; z2 <- nod$z[zz$j]
  a <- th[1] + L[1, 1] * z1; b <- th[2] + L[2, 1] * z1 + L[2, 2] * z2
  t <- g(a, b) + (z1^2 + z2^2) / 2 + log(nod$w[zz$i]) + log(nod$w[zz$j])
  m <- max(t); log(2 * pi) + log(L[1, 1] * L[2, 2]) + m + log(sum(exp(t - m)))
}
limitar <- function(par) { par[3:4] <- pmin(pmax(par[3:4], log(1e-3)), log(30)); par[5] <- pmin(pmax(par[5], -6.9), 6.9); par[1:2] <- pmin(pmax(par[1:2], -14), 14); par }
loglik <- function(par, dat, nod) {
  par <- limitar(par); mu <- par[1:2]; sA <- exp(par[3]); sB <- exp(par[4]); r <- tanh(par[5])
  S <- matrix(c(sA^2, r * sA * sB, r * sA * sB, sB^2), 2)
  r <- tryCatch(sum(mapply(function(vp, fn, vn, fp) ll_estudio(vp, fn, vn, fp, mu, S, nod) + lchoose(vp + fn, vp) + lchoose(vn + fp, vn),
                           dat$vp, dat$fn, dat$vn, dat$fp)), error = function(e) -1e12)
  if (!is.finite(r)) -1e12 else r
}
# Fuerza bruta: integral doble con integrate anidado, sobre tramos (un solo intervalo largo pierde el pico estrecho y subestima).
integ_tramos <- function(f, lo, hi, n = 30) { cortes <- seq(lo, hi, length.out = n + 1)
  sum(vapply(1:n, function(i) integrate(f, cortes[i], cortes[i + 1], rel.tol = 1e-13, abs.tol = 0, subdivisions = 200)$value, 0)) }
ll_bruta <- function(par, dat) {
  mu <- par[1:2]; sA <- exp(par[3]); sB <- exp(par[4]); r <- tanh(par[5])
  una <- function(vp, fn, vn, fp) {
    f_ext <- function(a) vapply(a, function(ai) {
      mc <- mu[2] + r * sB / sA * (ai - mu[1]); sc <- sB * sqrt(1 - r^2)
      int <- integ_tramos(function(b) exp(vn * lsig(b) + fp * lsig(-b)) * dnorm(b, mc, sc), mc - 12 * sc, mc + 12 * sc)
      exp(vp * lsig(ai) + fn * lsig(-ai)) * dnorm(ai, mu[1], sA) * int }, 0)
    log(integ_tramos(f_ext, mu[1] - 12 * sA, mu[1] + 12 * sA)) + lchoose(vp + fn, vp) + lchoose(vn + fp, vn)
  }
  sum(mapply(una, dat$vp, dat$fn, dat$vn, dat$fp))
}

# ---------------- REML normal: óptimo exacto con implementación propia ----------------
reml_normal <- function(par, y, v) {
  sA <- exp(min(max(par[1], log(1e-3)), log(30))); sB <- exp(min(max(par[2], log(1e-3)), log(30))); r <- tanh(min(max(par[3], -6.9), 6.9))
  S <- matrix(c(sA^2, r * sA * sB, r * sA * sB, sB^2), 2); k <- nrow(y); ld <- 0; sW <- matrix(0, 2, 2); sWy <- c(0, 0); Ws <- vector("list", k)
  for (i in 1:k) { V <- diag(v[i, ]) + S; dt <- det(V); if (!(dt > 0)) return(1e12); ld <- ld + log(dt); W <- solve(V); Ws[[i]] <- W; sW <- sW + W; sWy <- sWy + W %*% y[i, ] }
  mu <- solve(sW, sWy); q <- sum(sapply(1:k, function(i) { r2 <- y[i, ] - mu; as.numeric(t(r2) %*% Ws[[i]] %*% r2) }))
  0.5 * (ld + log(det(sW)) + q)
}

generar_datos <- function(k, muA, muB, sA, sB, rho, nEnf, nSan) {
  S <- matrix(c(sA^2, rho * sA * sB, rho * sA * sB, sB^2), 2); Z <- matrix(rnorm(2 * k), k) %*% chol(S)
  a <- muA + Z[, 1]; b <- muB + Z[, 2]; n1 <- pmax(5, round(runif(k, nEnf[1], nEnf[2]))); n0 <- pmax(5, round(runif(k, nSan[1], nSan[2])))
  vp <- rbinom(k, n1, plogis(a)); vn <- rbinom(k, n0, plogis(b))
  data.frame(vp = vp, fn = n1 - vp, vn = vn, fp = n0 - vn)
}

escenarios <- list()
cfg <- expand.grid(k = c(5, 8, 12, 20, 35), muA = c(0.8, 2, 3), muB = c(1, 2.5), sA = c(0.4, 1.3), sB = c(0.4, 1.1), rho = c(-0.8, -0.2, 0.5), stringsAsFactors = FALSE)
cfg <- cfg[sample(nrow(cfg), 30), ]
for (i in 1:30) {
  cc <- cfg[i, ]
  chico <- runif(1) < 0.4   # algunos conjuntos con estudios pequeños (producen celdas en cero)
  dat <- generar_datos(cc$k, cc$muA, cc$muB, cc$sA, cc$sB, cc$rho, if (chico) c(8, 40) else c(25, 300), if (chico) c(10, 60) else c(40, 500))
  escenarios[[i]] <- list(id = sprintf("biv%02d", i), datos = dat, tipo = if (chico) "pequeños" else "grandes")
}
data(AuditC)
escenarios[[31]] <- list(id = "auditc", datos = data.frame(vp = AuditC$TP, fn = AuditC$FN, vn = AuditC$TN, fp = AuditC$FP), tipo = "real (mada::AuditC)")

nod <- gh(21)
fprs <- c(0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.95)
salida <- list()
for (e in escenarios) {
  dat <- e$datos; k <- nrow(dat); cat(sprintf("%s (k=%d)\n", e$id, k)); flush.console()
  r <- list(id = e$id, tipo = e$tipo, k = k, datos = dat)
  # 1) mada
  fit <- tryCatch(reitsma(data.frame(TP = dat$vp, FN = dat$fn, FP = dat$fp, TN = dat$vn), correction.control = "single"), error = function(er) NULL)
  if (!is.null(fit)) {
    s <- summary(fit); cf <- s$coefficients
    r$mada <- list(coef = as.numeric(coef(fit)), vcov = as.numeric(vcov(fit)), Psi = as.numeric(fit$Psi),
                   sens = as.numeric(cf["sensitivity", c("Estimate", "95%ci.lb", "95%ci.ub")]), fpr = as.numeric(cf["false pos. rate", c("Estimate", "95%ci.lb", "95%ci.ub")]),
                   auc_fino = as.numeric(AUC(fit, fpr = (1:19999) / 20000)$AUC), fprs = fprs,
                   sroc_rg = as.numeric(sroc(fit, fpr = fprs)[, 2]), sroc_ingenua = as.numeric(sroc(fit, fpr = fprs, type = "naive")[, 2]),
                   elipse = unname(ROCellipse(fit, level = 0.95)[[1]]))
  }
  # 2) REML normal exacto (parámetros en (A = logit sens, B = logit esp))
  tp <- dat$vp; fnn <- dat$fn; fpp <- dat$fp; tn <- dat$vn
  c0 <- ifelse(tp == 0 | fnn == 0 | fpp == 0 | tn == 0, 0.5, 0); tp <- tp + c0; fnn <- fnn + c0; fpp <- fpp + c0; tn <- tn + c0
  y <- cbind(log(tp / fnn), log(tn / fpp)); v <- cbind(1 / tp + 1 / fnn, 1 / tn + 1 / fpp)
  mejor <- NULL
  for (st in list(c(log(0.8), log(0.8), 0), c(log(0.4), log(0.4), -1), c(log(1.2), log(0.6), 1))) {
    o <- optim(st, reml_normal, y = y, v = v, method = "BFGS", control = list(reltol = 1e-15, maxit = 1000))
    if (is.null(mejor) || o$value < mejor$value) mejor <- o
  }
  r$reml_exacto <- list(par = as.numeric(mejor$par), objetivo = mejor$value)
  # 3) binomial exacta (implementación propia independiente)
  ini <- c(if (!is.null(fit)) c(coef(fit)[1], -coef(fit)[2]) else c(1.5, 1.5), mejor$par)
  ob <- NULL
  for (st in list(ini, c(ini[1:2], log(0.8), log(0.8), 0), c(ini[1:2], log(0.5), log(0.5), -1))) {
    o <- tryCatch(optim(st, function(p) -loglik(p, dat, nod), method = "BFGS", control = list(reltol = 1e-14, maxit = 400)), error = function(er) NULL)
    if (!is.null(o) && (is.null(ob) || o$value < ob$value)) ob <- o
  }
  par <- limitar(ob$par)
  hs <- tryCatch(numDeriv::hessian(function(p) -loglik(p, dat, nod), par), error = function(er) NULL)
  vcovMu <- NULL
  if (!is.null(hs)) { inv <- tryCatch(solve(hs), error = function(er) NULL); if (!is.null(inv) && all(diag(inv) > 0)) vcovMu <- as.numeric(inv[1:2, 1:2]) }
  r$binomial <- list(par = as.numeric(par), loglik = -ob$value, vcovMu = vcovMu, limite = any(par[3:4] <= log(1e-3) + 1e-3) || abs(par[5]) >= 6.9 - 1e-3)
  # 4) glmer (Laplace)
  L <- rbind(data.frame(id = 1:k, y = dat$vp, n = dat$vp + dat$fn, sens = 1, spec = 0), data.frame(id = 1:k, y = dat$vn, n = dat$vn + dat$fp, sens = 0, spec = 1))
  m <- tryCatch(suppressWarnings(suppressMessages(glmer(cbind(y, n - y) ~ 0 + sens + spec + (0 + sens + spec | id), data = L, family = binomial, nAGQ = 1,
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))))), error = function(er) NULL)
  if (!is.null(m)) { vc <- VarCorr(m)$id; r$glmer <- list(fixef = as.numeric(fixef(m)), sd = as.numeric(attr(vc, "stddev")), cor = as.numeric(attr(vc, "correlation")[1, 2]), loglik = as.numeric(logLik(m))) }
  salida[[length(salida) + 1]] <- r
}
# 5) Fuerza bruta (integrate anidado) en un conjunto PEQUEÑO y varios puntos del espacio de parámetros: con recuentos grandes el pico de la
#    verosimilitud es tan estrecho que `integrate` lo pierde y subestima, así que no sirve como verdad ahí (sí lo es con pocos pacientes).
dat_b <- data.frame(vp = c(7, 12, 4, 15, 9, 20), fn = c(3, 5, 2, 9, 4, 6), vn = c(18, 25, 9, 30, 12, 22), fp = c(2, 6, 3, 4, 1, 9))
puntos_b <- list(c(1.2, 1.5, log(0.8), log(0.6), -0.5), c(2, 0.5, log(1.3), log(0.4), 0.6), c(0.5, 2.2, log(0.3), log(1.1), -0.9),
                 c(1.7, 1.7, log(0.05), log(0.05), 0), c(-0.4, 0.3, log(2), log(1.5), 0.2))
cat("fuerza bruta\n"); flush.console()
bruta <- lapply(puntos_b, function(p) list(par = p, loglik_bruta = ll_bruta(p, dat_b), loglik_gh21 = loglik(p, dat_b, nod)))

escribir <- toJSON(list(version = R.version.string, paquetes = list(mada = as.character(packageVersion("mada")), lme4 = as.character(packageVersion("lme4"))), conjuntos = salida,
                        bruta = list(datos = dat_b, puntos = bruta)),
                   digits = NA, auto_unbox = TRUE, na = "null", null = "null")
write(escribir, file.path("validacion", "referencias", "referencia-bivariado.json"))
cat("Escrito: validacion/referencias/referencia-bivariado.json (", length(salida), "conjuntos )\n")
