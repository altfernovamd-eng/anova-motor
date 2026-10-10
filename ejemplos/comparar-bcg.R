# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Fernando Ramírez Nova
# Comprueba en R (metafor) el mismo análisis que `ejemplos/ejemplo-bcg.mts`. Requiere R con los paquetes `metafor` y `metadat`.
#   Rscript ejemplos/comparar-bcg.R
library(metafor)
d <- read.csv("ejemplos/bcg.csv")
e <- escalc(measure = "OR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = d)
r <- rma(yi, vi, data = e, method = "REML")
# Intervalo de predicción como lo calcula aNova: t con k-2 gl (IntHout et al. 2014; en `meta`: method.predict = "HTS").
# Ojo: predict(r) de metafor usa la normal, y `meta` >= 8 usa por defecto t con k-1 gl; ver README, «Diferencias de convención».
pi <- exp(as.numeric(r$b) + c(-1, 1) * qt(0.975, r$k - 2) * sqrt(r$tau2 + r$se^2))
cat(sprintf("Estudios: %d\n", r$k))
cat(sprintf("OR combinado (efectos aleatorios, REML): %.4f  IC 95 %% [%.4f, %.4f]\n", exp(r$b), exp(r$ci.lb), exp(r$ci.ub)))
cat(sprintf("Intervalo de predicción (t, k-2 gl): [%.4f, %.4f]\n", pi[1], pi[2]))
cat(sprintf("Heterogeneidad: Q = %.2f (gl %d), tau2 = %.4f, I2 de metafor = %.1f %%  (aNova muestra el I2 clásico (Q-gl)/Q = %.1f %%)\n",
            r$QE, r$k - 1, r$tau2, r$I2, 100 * (r$QE - (r$k - 1)) / r$QE))
