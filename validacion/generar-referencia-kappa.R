# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Fernando Ramírez Nova
# Referencia de R (paquete `irr`) para κ de Cohen y κ de Fleiss: conjuntos de calificaciones sintéticos con semilla fija.
# `irr` se usa como CAJA NEGRA: solo se guardan sus resultados (kappa, z, p), no se lee ni se copia su código.
# `pExacto` = 2·pnorm(−|z|) de R con el z de irr: el p de irr se calcula con una resta que pierde cifras cuando p < 1e-8.
#
#   Rscript validacion/generar-referencia-kappa.R
#
# Escribe validacion/referencias/referencia-kappa.json, que lee pruebas/test-kappa-r.mts.

suppressMessages(library(irr, lib.loc = file.path(Sys.getenv("HOME"), "Library/R-anova")))
options(digits = 17, scipen = 999)
pdf(NULL)
set.seed(20261006)

fmt <- function(x) sprintf("%.17g", x)
cat_a_texto <- function(m) apply(m, c(1, 2), function(v) paste0("c", v))

casos <- list()

# --- Cohen: 2 evaluadores, 2 a 6 categorías, 8 a 250 sujetos; acuerdo variable (de bajo a casi perfecto) ---
for (i in 1:90) {
  k <- sample(2:6, 1)
  n <- sample(c(8:30, 40, 60, 100, 250), 1)
  acuerdo <- runif(1, 0.1, 0.98)
  probs <- rexp(k); probs <- probs / sum(probs)
  a <- sample(1:k, n, replace = TRUE, prob = probs)
  b <- ifelse(runif(n) < acuerdo, a, sample(1:k, n, replace = TRUE, prob = probs))
  m <- cbind(a, b)
  r <- tryCatch(kappa2(m), error = function(e) NULL)
  if (!is.null(r) && is.finite(r$value) && is.finite(r$statistic)) {
    casos[[length(casos) + 1]] <- list(tipo = "cohen", datos = m, kappa = r$value, z = r$statistic, p = r$p.value)
  }
}

# --- Fleiss: 3 a 8 evaluadores, 2 a 6 categorías, 6 a 150 sujetos ---
for (i in 1:90) {
  k <- sample(2:6, 1)
  n <- sample(c(6:30, 50, 80, 150), 1)
  m_ev <- sample(3:8, 1)
  acuerdo <- runif(1, 0.1, 0.95)
  probs <- rexp(k); probs <- probs / sum(probs)
  verdad <- sample(1:k, n, replace = TRUE, prob = probs)
  m <- sapply(1:m_ev, function(j) ifelse(runif(n) < acuerdo, verdad, sample(1:k, n, replace = TRUE, prob = probs)))
  r <- tryCatch(kappam.fleiss(m), error = function(e) NULL)
  if (!is.null(r) && is.finite(r$value) && is.finite(r$statistic)) {
    casos[[length(casos) + 1]] <- list(tipo = "fleiss", datos = m, kappa = r$value, z = r$statistic, p = r$p.value)
  }
}

# --- Casos límite fijos ---
perfecto <- cbind(c(1, 2, 3, 1, 2, 3), c(1, 2, 3, 1, 2, 3))
casos[[length(casos) + 1]] <- list(tipo = "cohen", datos = perfecto, kappa = kappa2(perfecto)$value, z = kappa2(perfecto)$statistic, p = kappa2(perfecto)$p.value)
pf <- cbind(c(1, 2, 3, 1, 2), c(1, 2, 3, 1, 2), c(1, 2, 3, 1, 2))
casos[[length(casos) + 1]] <- list(tipo = "fleiss", datos = pf, kappa = kappam.fleiss(pf)$value, z = kappam.fleiss(pf)$statistic, p = kappam.fleiss(pf)$p.value)
contrario <- cbind(c(1, 1, 1, 2, 2, 2, 1, 2), c(2, 2, 2, 1, 1, 1, 2, 1))   # peor que el azar
casos[[length(casos) + 1]] <- list(tipo = "cohen", datos = contrario, kappa = kappa2(contrario)$value, z = kappa2(contrario)$statistic, p = kappa2(contrario)$p.value)

mat <- function(m) paste0("[", paste(apply(m, 1, function(f) paste0("[", paste(f, collapse = ","), "]")), collapse = ","), "]")
txt <- vapply(casos, function(c) {
  paste0("{\"tipo\":\"", c$tipo, "\",\"datos\":", mat(c$datos), ",\"kappa\":", fmt(c$kappa), ",\"z\":", fmt(c$z), ",\"p\":", fmt(c$p),
         ",\"pExacto\":", fmt(2 * pnorm(-abs(c$z))), "}")
}, "")
writeLines(paste0("{\"irr\":\"", as.character(packageVersion("irr")), "\",\"casos\":[", paste(txt, collapse = ","), "]}"),
           "validacion/referencias/referencia-kappa.json")
cat(sprintf("Escrito: %d casos (%d Cohen, %d Fleiss)\n", length(casos),
            sum(vapply(casos, function(c) c$tipo == "cohen", TRUE)), sum(vapply(casos, function(c) c$tipo == "fleiss", TRUE))))
