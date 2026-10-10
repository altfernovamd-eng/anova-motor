# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Fernando Ramírez Nova
# Ayudantes compartidos por los generadores de referencia (`generar-referencia.R` y `generar-referencia-avanzada.R`).
# Se cargan con source() DESPUÉS de library(metafor), library(meta) y library(jsonlite).

# Estimadores iterativos con tolerancia muy estricta: se compara el ÓPTIMO matemático, no la regla de parada de cada
# programa (metafor se detiene por defecto con una tolerancia de 1e-5).
ctrl <- list(threshold = 1e-12, maxiter = 5000, stepadj = 0.5)

num <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) NULL else as.numeric(x)

# Q generalizado de Paule-Mandel/Q-profile: Σ (yᵢ − μ̂(τ²))² / (vᵢ + τ²), decreciente en τ². Su raíz se busca con una
# tolerancia de 1e-15: la referencia es el valor matemático exacto, no el que devuelve cada paquete tras su regla de
# parada (el `uniroot` interno de metafor para PM y para `confint` se detiene antes y puede diferir en la 4ª–5ª cifra).
qgen <- function(yi, vi, t2) { w <- 1 / (vi + t2); mu <- sum(w * yi) / sum(w); sum(w * (yi - mu)^2) }
raiz_qgen <- function(yi, vi, objetivo) {
  if (qgen(yi, vi, 0) <= objetivo) return(0)
  hi <- 1
  while (qgen(yi, vi, hi) > objetivo && hi < 1e12) hi <- hi * 2
  uniroot(function(t2) qgen(yi, vi, t2) - objetivo, c(0, hi), tol = 1e-15, maxiter = 1000)$root
}

# ------------------------------------------------------------------------------------------------
# Combinación (pooling) de efectos con su error estándar: todo lo que hace `calcularMetaAnalisis`.
# ------------------------------------------------------------------------------------------------
agrupar <- function(yi, sei) {
  vi <- sei^2
  k <- length(yi)
  fe <- rma(yi, vi, method = "FE")
  q <- fe$QE
  out <- list(
    k = k,
    fijo = list(est = num(fe$beta[1]), se = num(fe$se), lb = num(fe$ci.lb), ub = num(fe$ci.ub)),
    q = num(q), gl = k - 1, pq = num(fe$QEp),
    i2 = if (q > 0) max(0, (q - (k - 1)) / q * 100) else 0
  )
  for (m in c("DL", "PM", "REML")) {
    # PM: τ² exacto (ver `raiz_qgen`) y el resto se calcula con ese τ² fijo; DL y REML: el propio estimador de R.
    if (m == "PM") {
      t2 <- raiz_qgen(yi, vi, k - 1)
      r <- tryCatch(rma(yi, vi, tau2 = t2), error = function(e) NULL)
      kn <- tryCatch(rma(yi, vi, tau2 = t2, test = "knha"), error = function(e) NULL)
      pred <- tryCatch(metagen(TE = yi, seTE = sei, tau.preset = sqrt(t2), prediction = TRUE, method.predict = "HTS", hakn = FALSE),
                       error = function(e) NULL)
    } else {
      r <- tryCatch(rma(yi, vi, method = m, control = ctrl), error = function(e) NULL)
      kn <- tryCatch(rma(yi, vi, method = m, test = "knha", control = ctrl), error = function(e) NULL)
      pred <- tryCatch(metagen(TE = yi, seTE = sei, method.tau = m, prediction = TRUE, method.predict = "HTS", hakn = FALSE, control = ctrl),
                       error = function(e) NULL)
    }
    if (is.null(r)) { out[[m]] <- NULL; next }
    o <- list(
      tau2 = num(r$tau2),
      est = num(r$beta[1]), se = num(r$se), lb = num(r$ci.lb), ub = num(r$ci.ub),
      hk_lb = if (is.null(kn)) NULL else num(kn$ci.lb),
      hk_ub = if (is.null(kn)) NULL else num(kn$ci.ub),
      pred_lb = if (is.null(pred) || k < 3) NULL else num(pred$lower.predict),
      pred_ub = if (is.null(pred) || k < 3) NULL else num(pred$upper.predict)
    )
    out[[m]] <- o
  }
  # IC de τ² e I² por Q-profile (Viechtbauer 2007), con raíces exactas
  if (k >= 3) {
    w <- 1 / vi
    s2 <- (k - 1) * sum(w) / (sum(w)^2 - sum(w^2))
    lo <- raiz_qgen(yi, vi, qchisq(0.975, k - 1))
    hi <- raiz_qgen(yi, vi, qchisq(0.025, k - 1))
    out$tau2_ic <- c(lo, hi)
    out$i2_ic <- c(100 * lo / (lo + s2), 100 * hi / (hi + s2))
  }
  out
}


# mada::madauni (caja negra) para el DOR combinado de una tabla de precisión diagnóstica: ln(DOR) y su EE, τ² de DerSimonian-Laird y Q,
# con la corrección de continuidad «single» (solo el estudio con una celda en cero), que es la misma de aNova. `madauni` solo admite
# DerSimonian-Laird ("DSL") y Mantel-Haenszel para el DOR; el DOR de efectos fijos de aNova es de inverso-varianza, así que no se compara.
mada_dor <- function(vp, fn, fp, vn) {
  m <- tryCatch(suppressWarnings(mada::madauni(data.frame(TP = vp, FN = fn, FP = fp, TN = vn), type = "DOR", method = "DSL", correction.control = "single")),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  list(lnDOR = as.numeric(coef(m)), se = sqrt(as.numeric(vcov(m))), tau2 = as.numeric(m$tau_sq), Q = as.numeric(m$CQ[1]), gl = as.numeric(m$CQ[3]))
}
