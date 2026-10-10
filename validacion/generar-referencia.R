# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Fernando Ramírez Nova
# Genera los datos de REFERENCIA para validar el motor estadístico de aNova contra R (`metafor` y `meta`).
#
# Qué hace: crea cientos de conjuntos de datos sintéticos con semilla fija (efectos genéricos, tablas 2×2,
# medias, proporciones), los analiza con R y guarda entradas + resultados en `referencias/referencia.json`.
# `comparar.mts` alimenta esos mismos datos al motor de aNova y compara cada cantidad.
#
# Para qué sirve guardar el resultado: la comparación (`npm run validacion`) NO necesita R. R solo hace falta para
# regenerar la referencia (`npm run validacion:regenerar`), p. ej. con otra semilla o con una versión nueva de R.
#
# Requisitos: R ≥ 4.3 con los paquetes metafor, meta y jsonlite. Se buscan primero en ~/Library/R-anova.
#   Rscript validacion/generar-referencia.R [semilla]

lib <- file.path(Sys.getenv("HOME"), "Library", "R-anova")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
pdf(NULL)  # baujat() dibuja aunque se pida plot = FALSE; así no deja Rplots.pdf
suppressPackageStartupMessages({
  library(metafor)
  library(meta)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
semilla <- if (length(args) >= 1) as.integer(args[1]) else 20261005L
set.seed(semilla)

# Ayudantes compartidos (ctrl, num, qgen, raiz_qgen, agrupar).
source(file.path("validacion", "ayudas.R"))

# ------------------------------------------------------------------------------------------------
# Sesgo de publicación, influencia y subgrupos, sobre un conjunto (yi, sei) con año y grupo.
# ------------------------------------------------------------------------------------------------
adicionales <- function(yi, sei, anio, grupo) {
  vi <- sei^2
  k <- length(yi)
  out <- list()
  fe <- rma(yi, vi, method = "FE")
  if (k >= 3) {
    rt <- tryCatch(regtest(fe, model = "lm", predictor = "sei"), error = function(e) NULL)
    if (!is.null(rt)) {
      cf <- coef(summary(rt$fit))
      out$egger <- list(intercepto = num(cf["Xsei", "Estimate"]), se = num(cf["Xsei", "Std. Error"]),
                        t = num(cf["Xsei", "t value"]), p = num(cf["Xsei", "Pr(>|t|)"]))
    }
    bg <- tryCatch(ranktest(fe), error = function(e) NULL)
    if (!is.null(bg)) out$begg <- list(tau = num(bg$tau), p = num(bg$pval))
    bj <- tryCatch(suppressWarnings(baujat(fe, plot = FALSE)), error = function(e) NULL)
    if (!is.null(bj)) out$baujat <- list(x = as.numeric(bj$x), y = as.numeric(bj$y))
    # Diagnóstico de influencia (Viechtbauer y Cheung 2010) y gráfico radial/Galbraith, con los tres modelos que ofrece aNova.
    # `influence()` y `radial()` se usan como caja negra: solo se guardan sus resultados.
    infl <- list()
    for (modelo in c("FE", "DL", "REML")) {
      r <- tryCatch(rma(yi, vi, method = modelo, control = ctrl), error = function(e) NULL)
      if (is.null(r)) next
      inf <- tryCatch(suppressWarnings(influence(r)), error = function(e) NULL)
      if (is.null(inf)) next
      d <- as.data.frame(inf$inf)
      rd <- tryCatch(radial(r), error = function(e) NULL)
      infl[[modelo]] <- list(rstudent = as.numeric(d$rstudent), dffits = as.numeric(d$dffits), cook = as.numeric(d$cook.d),
                             covr = as.numeric(d$cov.r), tau2del = as.numeric(d$tau2.del), qedel = as.numeric(d$QE.del),
                             hat = as.numeric(d$hat), weight = as.numeric(d$weight),
                             dfbs = as.numeric(as.data.frame(inf$dfbs)[[1]]), marcado = as.character(d$inf) == "*",
                             radial_x = if (is.null(rd)) NULL else as.numeric(rd$x), radial_z = if (is.null(rd)) NULL else as.numeric(rd$y))
    }
    out$influencia <- infl
    # Trim-and-fill: aNova estima k0 con la media de efectos fijos (Duval-Tweedie) y combina con el modelo elegido.
    tfs <- list()
    for (modelo in c("FE", "DL", "REML")) {
      r <- tryCatch(rma(yi, vi, method = modelo, control = ctrl), error = function(e) NULL)
      if (is.null(r)) next
      for (lado in c("left", "right")) {
        tf <- tryCatch(suppressWarnings(trimfill(r, estimator = "L0", side = lado)), error = function(e) NULL)
        if (!is.null(tf)) {
          # trimfill() reajusta el modelo con el control por defecto de metafor (REML: umbral 1e-5), por eso su efecto ajustado
          # con REML solo es fiable a ~1e-6; se reajusta con el control estricto sobre los mismos datos aumentados.
          ref_estricto <- tryCatch(rma(as.numeric(tf$yi), as.numeric(tf$vi), method = modelo, control = ctrl), error = function(e) NULL)
          tfs[[paste0(modelo, "_", lado)]] <- list(k0 = tf$k0, est = num(tf$beta[1]), lb = num(tf$ci.lb), ub = num(tf$ci.ub),
                                                   est_estricto = if (is.null(ref_estricto)) NULL else num(ref_estricto$beta[1]),
                                                   lb_estricto = if (is.null(ref_estricto)) NULL else num(ref_estricto$ci.lb),
                                                   ub_estricto = if (is.null(ref_estricto)) NULL else num(ref_estricto$ci.ub))
        }
      }
    }
    out$trimfill <- tfs
    # Leave-one-out
    lo <- list()
    for (modelo in c("FE", "DL", "REML")) {
      r <- tryCatch(rma(yi, vi, method = modelo, control = ctrl), error = function(e) NULL)
      if (is.null(r)) next
      l <- tryCatch(leave1out(r), error = function(e) NULL)
      if (!is.null(l)) lo[[modelo]] <- list(est = as.numeric(l$estimate), lb = as.numeric(l$ci.lb), ub = as.numeric(l$ci.ub),
                                             tau2 = as.numeric(l$tau2), i2 = as.numeric(l$I2), q = as.numeric(l$Q))
    }
    out$leave1out <- lo
  }
  # Acumulado (ordenado por año; empates por el orden de captura)
  ord <- order(anio, seq_along(anio))
  cu <- list()
  for (modelo in c("FE", "DL", "REML")) {
    r <- tryCatch(rma(yi, vi, method = modelo, control = ctrl), error = function(e) NULL)
    if (is.null(r)) next
    cm <- tryCatch(cumul(r, order = anio), error = function(e) NULL)
    if (!is.null(cm)) cu[[modelo]] <- list(orden = as.integer(ord), est = as.numeric(cm$estimate), lb = as.numeric(cm$ci.lb), ub = as.numeric(cm$ci.ub))
  }
  out$acumulado <- cu
  # Subgrupos (tau² propio por subgrupo; Q entre subgrupos)
  if (length(unique(grupo)) >= 2) {
    sg <- list()
    for (m in c("DL", "REML")) {
      s <- tryCatch(
        metagen(TE = yi, seTE = sei, subgroup = grupo, tau.common = FALSE, method.tau = m, hakn = FALSE,
                common = TRUE, random = TRUE, control = ctrl),
        error = function(e) NULL
      )
      if (is.null(s)) next
      sg[[m]] <- list(
        niveles = as.character(s$subgroup.levels), k = as.integer(s$k.w),
        est_aleatorio = as.numeric(s$TE.random.w), lb_aleatorio = as.numeric(s$lower.random.w), ub_aleatorio = as.numeric(s$upper.random.w),
        est_fijo = as.numeric(s$TE.common.w), lb_fijo = as.numeric(s$lower.common.w), ub_fijo = as.numeric(s$upper.common.w),
        qb_aleatorio = num(s$Q.b.random), qb_fijo = num(s$Q.b.common), gl = num(s$df.Q.b),
        p_aleatorio = num(s$pval.Q.b.random), p_fijo = num(s$pval.Q.b.common)
      )
    }
    out$subgrupos <- sg
  }
  out
}

# ------------------------------------------------------------------------------------------------
# Generadores de conjuntos de datos
# ------------------------------------------------------------------------------------------------
disenos <- c("ECA", "Cohortes", "Casos y Controles")
nuevo_id <- local({ n <- 0; function(pref) { n <<- n + 1; sprintf("%s%03d", pref, n) } })

# 1) Efectos genéricos con su error estándar (sirve para validar el pooling aislado de cómo se calculó cada efecto)
escenarios_iv <- list()
for (k in c(2, 3, 4, 5, 8, 12, 20, 40)) {
  for (tau2 in c(0, 0.02, 0.2, 1)) {
    for (rep in 1:2) {
      sei <- runif(k, 0.05, 0.7)
      mu <- runif(1, -1, 1)
      yi <- rnorm(k, mu, sqrt(tau2 + sei^2))
      anio <- sample(2000:2020, k, replace = TRUE)
      grupo <- sample(disenos, k, replace = TRUE)
      escenarios_iv[[length(escenarios_iv) + 1]] <- list(
        id = nuevo_id("iv"), k = k, tau2_verdadero = tau2, yi = yi, sei = sei, anio = anio, grupo = grupo,
        ref = c(list(pool = agrupar(yi, sei)), adicionales(yi, sei, anio, grupo))
      )
    }
  }
}

# 2) Tablas 2×2 (OR, RR, RD; Mantel-Haenszel; Peto; Peters)
escenarios_bin <- list()
configs <- expand.grid(k = c(3, 6, 12, 25), pc = c(0.01, 0.05, 0.2, 0.5), lor = c(-0.7, 0, 0.5), tau = c(0, 0.4))
configs <- configs[sample(nrow(configs), 70), ]
for (i in seq_len(nrow(configs))) {
  cfg <- configs[i, ]
  k <- cfg$k
  n1 <- round(runif(k, 25, 400)); n2 <- pmax(10, n1 + round(rnorm(k, 0, 20)))
  logor <- rnorm(k, cfg$lor, cfg$tau)
  pc <- pmin(0.95, pmax(0.003, plogis(qlogis(cfg$pc) + rnorm(k, 0, 0.3))))
  pt <- plogis(qlogis(pc) + logor)
  ci <- rbinom(k, n2, pc); ai <- rbinom(k, n1, pt)
  anio <- sample(2000:2020, k, replace = TRUE); grupo <- sample(disenos, k, replace = TRUE)
  bi <- n1 - ai; di <- n2 - ci
  ref <- list()
  for (medida in c("OR", "RR", "RD")) {
    es <- escalc(measure = medida, ai = ai, bi = bi, ci = ci, di = di, add = 1/2, to = "only0", drop00 = FALSE)
    yi <- as.numeric(es$yi); sei <- sqrt(as.numeric(es$vi))
    # Cochrane Handbook 10.4.4.2: las tablas sin eventos en ningún brazo (o con el evento en todos) no informan el OR ni el RR y se
    # excluyen de la combinación; en la diferencia de riesgos sí entran (RD = 0). Es el `drop00` de metafor. `yi`/`sei` por estudio
    # se guardan completos (para comparar el cálculo de cada tabla); lo que se combina es lo que queda tras excluirlas.
    sin_info <- (ai == 0 & ci == 0) | (bi == 0 & di == 0)
    conserva <- if (medida == "RD") rep(TRUE, k) else !sin_info
    ref[[medida]] <- c(list(yi = yi, sei = sei, conserva = conserva, pool = agrupar(yi[conserva], sei[conserva])),
                       adicionales(yi[conserva], sei[conserva], anio[conserva], grupo[conserva]))
    # Mantel-Haenszel (OR/RR/RD). aNova NO corrige celdas cero en el estimador MH y descarta las tablas doble-cero,
    # que es la forma estándar del método; se configura R igual (add = 0, drop00 = TRUE) y se documenta.
    mh <- tryCatch(suppressWarnings(rma.mh(ai = ai, bi = bi, ci = ci, di = di, measure = medida)), error = function(e) NULL)
    if (!is.null(mh)) ref[[medida]]$mh <- list(est = num(mh$beta[1]), se = num(mh$se), lb = num(mh$ci.lb), ub = num(mh$ci.ub),
                                              pz = num(mh$pval), k = mh$k, q = num(mh$QE))
    # Diferencia de riesgos: metafor usa el estimador de varianza de Sato-Greenland-Robins (1989) y excluye del Q las tablas
    # doble-cero (drop00 = TRUE). RevMan 5 / Cochrane Handbook usan Greenland-Robins (1985) y mantienen las doble-cero con +0,5.
    # aNova sigue a RevMan, así que se transcriben aquí las fórmulas literales de «Statistical methods in RevMan 5» (ecs. 5-6):
    #   SE = sqrt(J)/K ; J = Σ (a b n2³ + c d n1³)/(n1 n2 N²) ; K = Σ n1 n2/N ; Q = Σ (RDi − RD_MH)² / SEi²  (+0,5 si hay celda cero)
    if (medida == "RD") {
      n1 <- ai + bi; n2 <- ci + di; N <- n1 + n2
      K <- sum(n1 * n2 / N); J <- sum((ai * bi * n2^3 + ci * di * n1^3) / (n1 * n2 * N^2))
      est <- sum(n1 * n2 / N * (ai / n1 - ci / n2)) / K; se <- sqrt(J) / K
      cero <- ai == 0 | bi == 0 | ci == 0 | di == 0
      a2 <- ai + 0.5 * cero; b2 <- bi + 0.5 * cero; c2 <- ci + 0.5 * cero; d2 <- di + 0.5 * cero
      m1 <- a2 + b2; m2 <- c2 + d2
      rdi <- a2 / m1 - c2 / m2; vi <- a2 * b2 / m1^3 + c2 * d2 / m2^3
      ref[[medida]]$mh_revman <- list(est = est, se = se, lb = est - qnorm(0.975) * se, ub = est + qnorm(0.975) * se,
                                      q = sum((rdi - est)^2 / vi))
    }
  }
  pe <- tryCatch(suppressWarnings(rma.peto(ai = ai, bi = bi, ci = ci, di = di)), error = function(e) NULL)
  if (!is.null(pe)) ref$peto <- list(est = num(pe$beta[1]), se = num(pe$se), lb = num(pe$ci.lb), ub = num(pe$ci.ub), k = pe$k, q = num(pe$QE))
  # Peters: regresión de ln(OR) sobre 1/n con pesos 1/(1/(a+c)+1/(b+d)) (meta::metabias, method.bias = "Peters")
  mb <- tryCatch(suppressWarnings(metabin(event.e = ai, n.e = n1, event.c = ci, n.c = n2, sm = "OR", method = "Inverse", incr = 0.5)),
                 error = function(e) NULL)
  if (!is.null(mb) && k >= 3) {
    pt <- tryCatch(metabias(mb, method.bias = "Peters", k.min = 3), error = function(e) NULL)
    if (!is.null(pt)) ref$peters <- list(pendiente = num(pt$estimate[["bias"]]), se = num(pt$estimate[["se.bias"]]), p = num(pt$p.value),
                                        t = num(pt$statistic))
  }
  escenarios_bin[[length(escenarios_bin) + 1]] <- list(id = nuevo_id("bin"), k = k, ai = ai, n1 = n1, ci = ci, n2 = n2,
                                                       anio = anio, grupo = grupo, ref = ref)
}

# 3) Medias (DM y g de Hedges)
escenarios_cont <- list()
for (i in 1:40) {
  k <- sample(c(3, 5, 10, 20), 1)
  n1 <- round(runif(k, 10, 200)); n2 <- round(runif(k, 10, 200))
  sd1 <- runif(k, 0.5, 5); sd2 <- sd1 * runif(k, 0.8, 1.25)
  mu <- runif(1, -1, 1); tau <- sample(c(0, 0.3), 1)
  m2 <- rnorm(k, 10, 2); m1 <- m2 + rnorm(k, mu, tau + 0.2) * sd1
  ref <- list()
  for (medida in c("MD", "SMD")) {
    es <- escalc(measure = medida, m1i = m1, sd1i = sd1, n1i = n1, m2i = m2, sd2i = sd2, n2i = n2)
    yi <- as.numeric(es$yi); sei <- sqrt(as.numeric(es$vi))
    ref[[medida]] <- list(yi = yi, sei = sei, pool = agrupar(yi, sei))
  }
  escenarios_cont[[length(escenarios_cont) + 1]] <- list(id = nuevo_id("cont"), k = k, m1 = m1, sd1 = sd1, n1 = n1, m2 = m2, sd2 = sd2, n2 = n2, ref = ref)
}

# 4) Proporciones de un solo brazo (logit)
escenarios_prop <- list()
for (i in 1:40) {
  k <- sample(c(3, 6, 12, 25), 1)
  ni <- round(runif(k, 15, 600))
  p <- plogis(rnorm(k, qlogis(sample(c(0.02, 0.1, 0.4, 0.8), 1)), sample(c(0, 0.5), 1)))
  xi <- rbinom(k, ni, p)
  es <- escalc(measure = "PLO", xi = xi, ni = ni, add = 1/2, to = "only0")
  yi <- as.numeric(es$yi); sei <- sqrt(as.numeric(es$vi))
  escenarios_prop[[length(escenarios_prop) + 1]] <- list(id = nuevo_id("prop"), k = k, xi = xi, ni = ni,
                                                         ref = list(PLO = list(yi = yi, sei = sei, pool = agrupar(yi, sei))))
}

# 5) Tasas de incidencia (razón y diferencia de tasas): metafor::escalc(IRR / IRD, add = 1/2, to = "only0")
escenarios_tasas <- list()
for (i in 1:40) {
  k <- sample(c(3, 5, 10, 20), 1)
  t1 <- round(runif(k, 200, 6000)); t2 <- round(runif(k, 200, 6000))
  lam0 <- sample(c(0.0002, 0.001, 0.005, 0.02), 1); rr <- rnorm(k, sample(c(-0.5, 0, 0.4), 1), sample(c(0, 0.3), 1))
  x2 <- rpois(k, lam0 * t2 * exp(rnorm(k, 0, 0.3))); x1 <- rpois(k, lam0 * t1 * exp(rr))
  anio <- sample(2000:2020, k, replace = TRUE); grupo <- sample(disenos, k, replace = TRUE)
  ref <- list()
  for (medida in c("IRR", "IRD")) {
    es <- escalc(measure = medida, x1i = x1, t1i = t1, x2i = x2, t2i = t2, add = 1/2, to = "only0")
    yi <- as.numeric(es$yi); sei <- sqrt(as.numeric(es$vi))
    # Sin eventos en ningún brazo no se informa la razón: se excluye del OR-like (IRR); la diferencia de tasas (IRD) lo conserva.
    conserva <- if (medida == "IRR") !(x1 == 0 & x2 == 0) else rep(TRUE, k)
    # La diferencia de tasas se combina «por 1000 personas-tiempo» (factor 1000): con efectos del orden de 0,001 el REML de R usa una
    # tolerancia ABSOLUTA y se detiene antes de converger (su τ² difiere del óptimo exacto en ~6e-4 relativo). A esa escala R es fiable,
    # y el REML a escala pequeña se comprueba aparte contra el óptimo exacto (test-regresion-r).
    f <- if (medida == "IRD") 1000 else 1
    if (sum(conserva) >= 2) {
      ref[[medida]] <- c(list(yi = yi, sei = sei, factor = f, conserva = conserva, pool = agrupar(f * yi[conserva], f * sei[conserva])),
                         adicionales(f * yi[conserva], f * sei[conserva], anio[conserva], grupo[conserva]))
    } else {
      ref[[medida]] <- list(yi = yi, sei = sei, factor = f, conserva = conserva)  # con menos de 2 estudios no hay meta-análisis
    }
  }
  escenarios_tasas[[length(escenarios_tasas) + 1]] <- list(id = nuevo_id("tasa"), k = k, x1 = x1, t1 = t1, x2 = x2, t2 = t2, anio = anio, grupo = grupo, ref = ref)
}

# 6) Razón de medias. Fórmula del Handbook (Friedrich 2008), sin la corrección de sesgo de metafor: se calcula con aritmética base y se
#    guarda aparte el yi de metafor (con corrección) para medir la diferencia de convención.
escenarios_rom <- list()
for (i in 1:40) {
  k <- sample(c(3, 5, 10, 20), 1)
  n1 <- round(runif(k, 10, 200)); n2 <- round(runif(k, 10, 200))
  mu2 <- runif(k, 5, 50); cv <- runif(k, 0.15, 0.6)
  theta <- rnorm(k, sample(c(-0.3, 0, 0.25), 1), sample(c(0, 0.2), 1))   # efecto verdadero de cada estudio (τ = 0 → todos iguales)
  sd1 <- mu2 * exp(theta) * cv; sd2 <- mu2 * cv * runif(k, 0.8, 1.25)
  # Medias observadas = media verdadera + error de muestreo (así los efectos varían aunque τ = 0, y siguen siendo positivas)
  m2 <- pmax(0.2 * mu2, rnorm(k, mu2, sd2 / sqrt(n2))); m1 <- pmax(0.2 * mu2, rnorm(k, mu2 * exp(theta), sd1 / sqrt(n1)))
  anio <- sample(2000:2020, k, replace = TRUE); grupo <- sample(disenos, k, replace = TRUE)
  yi <- log(m1 / m2); sei <- sqrt(sd1^2 / (n1 * m1^2) + sd2^2 / (n2 * m2^2))
  yim <- as.numeric(escalc(measure = "ROM", m1i = m1, sd1i = sd1, n1i = n1, m2i = m2, sd2i = sd2, n2i = n2)$yi)
  escenarios_rom[[length(escenarios_rom) + 1]] <- list(id = nuevo_id("rom"), k = k, m1 = m1, sd1 = sd1, n1 = n1, m2 = m2, sd2 = sd2, n2 = n2,
                                                       anio = anio, grupo = grupo, yi_metafor = yim,
                                                       ref = list(ROM = c(list(yi = yi, sei = sei, pool = agrupar(yi, sei)), adicionales(yi, sei, anio, grupo))))
}

# 7) Correlaciones (z de Fisher): escalc(ZCOR)
escenarios_cor <- list()
for (i in 1:40) {
  k <- sample(c(3, 5, 10, 20), 1)
  ni <- round(runif(k, 12, 300)); mu <- sample(c(-0.4, 0.1, 0.5, 0.8), 1)
  ri <- tanh(rnorm(k, atanh(mu), sample(c(0.05, 0.2), 1) + 1 / sqrt(ni - 3)))
  anio <- sample(2000:2020, k, replace = TRUE); grupo <- sample(disenos, k, replace = TRUE)
  es <- escalc(measure = "ZCOR", ri = ri, ni = ni)
  yi <- as.numeric(es$yi); sei <- sqrt(as.numeric(es$vi))
  escenarios_cor[[length(escenarios_cor) + 1]] <- list(id = nuevo_id("cor"), k = k, ri = ri, ni = ni, anio = anio, grupo = grupo,
                                                       ref = list(ZCOR = c(list(yi = yi, sei = sei, pool = agrupar(yi, sei)), adicionales(yi, sei, anio, grupo))))
}

# 8) Tiempo hasta el evento con O−E y V: ln HR = (O−E)/V, EE = 1/√V; el modelo fijo con varianza inversa debe dar Σ(O−E)/ΣV (Peto)
escenarios_oev <- list()
for (i in 1:40) {
  k <- sample(c(3, 5, 10, 20), 1)
  v <- runif(k, 3, 60); th <- rnorm(k, sample(c(-0.5, 0, 0.3), 1), sample(c(0, 0.25), 1))
  oe <- rnorm(k, v * th, sqrt(v))
  anio <- sample(2000:2020, k, replace = TRUE); grupo <- sample(disenos, k, replace = TRUE)
  yi <- oe / v; sei <- 1 / sqrt(v)
  escenarios_oev[[length(escenarios_oev) + 1]] <- list(id = nuevo_id("oev"), k = k, oe = oe, v = v, anio = anio, grupo = grupo,
                                                       peto = list(est = sum(oe) / sum(v), se = 1 / sqrt(sum(v))),
                                                       ref = list(HR = c(list(yi = yi, sei = sei, pool = agrupar(yi, sei)), adicionales(yi, sei, anio, grupo))))
}

resultado <- list(
  meta = list(
    fecha = format(Sys.time(), "%Y-%m-%d"), semilla = semilla,
    R = R.version.string,
    paquetes = list(metafor = as.character(packageVersion("metafor")), meta = as.character(packageVersion("meta")),
                    jsonlite = as.character(packageVersion("jsonlite"))),
    escenarios = list(iv = length(escenarios_iv), bin = length(escenarios_bin), cont = length(escenarios_cont), prop = length(escenarios_prop),
                      tasas = length(escenarios_tasas), rom = length(escenarios_rom), cor = length(escenarios_cor), oev = length(escenarios_oev))
  ),
  iv = escenarios_iv, bin = escenarios_bin, cont = escenarios_cont, prop = escenarios_prop,
  tasas = escenarios_tasas, rom = escenarios_rom, cor = escenarios_cor, oev = escenarios_oev
)

destino <- file.path("validacion", "referencias")
dir.create(destino, recursive = TRUE, showWarnings = FALSE)
write(toJSON(resultado, digits = NA, auto_unbox = TRUE, na = "null", null = "null"), file.path(destino, "referencia.json"))
cat("Referencia escrita:", file.path(destino, "referencia.json"), "\n")
print(resultado$meta$escenarios)
