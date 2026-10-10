# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Fernando Ramírez Nova
# Referencia de R para el modelo multinivel con varianza robusta (multinivel.ts).
#
# Qué se guarda para cada conjunto (semilla fija; 60 simulados con estudios de 1 a 12 efectos, correlación supuesta, varianzas entre estudios
# y entre efectos, escalas de varianza y moderadores variados):
#  1. metafor::rma.mv (REML) como CAJA NEGRA: componentes de varianza, coeficientes, covarianza del modelo y log-verosimilitud restringida.
#     El modelo es el CHE: random = ~ 1 | estudio/efecto, con la matriz de muestreo de metafor::vcalc (correlación ρ supuesta entre efectos).
#  2. clubSandwich::coef_test y conf_int (CR2, Satterthwaite) como CAJA NEGRA: error estándar robusto, grados de libertad, p e IC.
#     `comparar-multinivel.mts` evalúa además la parte robusta de aNova con los componentes de varianza de metafor, para aislar el
#     estimador sándwich del optimizador.
#
#   Rscript validacion/generar-referencia-multinivel.R
# Escribe validacion/referencias/referencia-multinivel.json (lo lee validacion/comparar-multinivel.mts).

lib <- file.path(Sys.getenv("HOME"), "Library", "R-anova")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({ library(metafor); library(clubSandwich); library(jsonlite) })
set.seed(20261008)
options(digits = 17)

simular <- function(m, tam, rho, su, sw, escalaV, mod) {
  d <- data.frame(estudio = rep(seq_len(m), tam))
  d$efecto <- seq_len(nrow(d))
  d$v <- runif(nrow(d), escalaV[1], escalaV[2])
  u <- rnorm(m, 0, su)[d$estudio]
  # error de muestreo con la correlación rho dentro del estudio
  e <- numeric(nrow(d))
  for (s in seq_len(m)) { i <- which(d$estudio == s); n <- length(i)
    R <- matrix(rho, n, n); diag(R) <- 1; Sg <- diag(sqrt(d$v[i]), n) %*% R %*% diag(sqrt(d$v[i]), n)
    e[i] <- as.numeric(t(chol(Sg)) %*% rnorm(n)) }
  d$x1 <- if (mod == "efecto") rnorm(nrow(d)) else if (mod == "estudio") (d$estudio %% 2) else NA
  b1 <- if (mod == "ninguno") 0 else 0.25
  d$y <- 0.3 + ifelse(is.na(d$x1), 0, b1 * d$x1) + u + rnorm(nrow(d), 0, sw) + e
  d
}

cfg <- expand.grid(m = c(4, 6, 10, 20, 40), rho = c(0, 0.3, 0.6, 0.8), su = c(0, 0.2, 0.5), sw = c(0, 0.15, 0.4), mod = c("ninguno", "efecto", "estudio"), stringsAsFactors = FALSE)
cfg <- cfg[sample(nrow(cfg), 130), ]
salida <- list(); n_ok <- 0
for (i in seq_len(nrow(cfg))) {
  if (n_ok >= 60) break
  cc <- cfg[i, ]
  tam <- sample(c(1, 1, 2, 2, 3, 4, 6), cc$m, replace = TRUE); if (all(tam == 1)) tam[1] <- 3
  if (i %% 15 == 0) tam[1] <- 12
  escala <- if (i %% 4 == 0) c(0.002, 0.02) else if (i %% 4 == 1) c(0.02, 0.15) else if (i %% 4 == 2) c(0.05, 0.4) else c(0.01, 0.08)
  d <- simular(cc$m, tam, cc$rho, cc$su, cc$sw, escala, cc$mod)
  if (cc$mod == "estudio" && length(unique(d$x1)) < 2) next
  V <- vcalc(v, cluster = estudio, obs = efecto, data = d, rho = cc$rho)
  m <- tryCatch(suppressWarnings(if (cc$mod == "ninguno") rma.mv(y, V, random = ~ 1 | estudio/efecto, data = d, method = "REML")
                                  else rma.mv(y, V, mods = ~ x1, random = ~ 1 | estudio/efecto, data = d, method = "REML")), error = function(e) NULL)
  if (is.null(m)) next
  ct <- tryCatch(coef_test(m, vcov = "CR2", cluster = d$estudio), error = function(e) NULL)
  ci <- tryCatch(conf_int(m, vcov = "CR2", cluster = d$estudio), error = function(e) NULL)
  if (is.null(ct) || is.null(ci)) next
  n_ok <- n_ok + 1
  cat(sprintf("mn%02d (m=%d, N=%d, rho=%.1f, mod=%s)\n", n_ok, cc$m, nrow(d), cc$rho, cc$mod)); flush.console()
  salida[[n_ok]] <- list(id = sprintf("mn%02d", n_ok), rho = cc$rho, moderador = cc$mod,
    datos = data.frame(estudio = d$estudio, y = d$y, v = d$v, x1 = d$x1),
    sigma2 = as.numeric(m$sigma2), beta = as.numeric(coef(m)), vcov = as.numeric(vcov(m)), ll = as.numeric(logLik(m)),
    se_cr2 = ct$SE, gl = ct$df_Satt, p = ct$p_Satt, ic_inf = ci$CI_L, ic_sup = ci$CI_U)
}
escribir <- toJSON(list(version = R.version.string, paquetes = list(metafor = as.character(packageVersion("metafor")), clubSandwich = as.character(packageVersion("clubSandwich"))), conjuntos = salida),
                   digits = NA, auto_unbox = TRUE, na = "null", null = "null")
write(escribir, file.path("validacion", "referencias", "referencia-multinivel.json"))
cat("Escrito: validacion/referencias/referencia-multinivel.json (", length(salida), "conjuntos )\n")
