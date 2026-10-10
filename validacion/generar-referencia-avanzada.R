# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Fernando Ramírez Nova
# Referencia de R para los métodos «avanzados» de aNova: META-REGRESIÓN, SROC (precisión diagnóstica) y META-ANÁLISIS EN RED.
#
# Igual que `generar-referencia.R`: conjuntos sintéticos con semilla fija, analizados con R, guardados en
# `referencias/referencia-avanzada.json`; `comparar.mts` los pasa por el motor real de aNova y compara número a número.
# `npm run validacion` no necesita R; R solo se usa para regenerar (`npm run validacion:regenerar`).
#
# Requisitos: R ≥ 4.3 con metafor, meta, netmeta, mada y jsonlite (se buscan primero en ~/Library/R-anova).
#   Rscript validacion/generar-referencia-avanzada.R [semilla]
#
# MODO MASIVO (para comprobar una reescritura o un ajuste del motor con muchas pruebas): con la variable de entorno AVANZADA_MASIVA=1 se
# omite la meta-regresión, se generan 400 conjuntos de SROC y 300 redes (con estudios de 2, 3 y 4 brazos y una familia con celdas en
# cero), se añade la comparación con `mada::madauni` y todo se escribe en `referencias/referencia-red-sroc-masiva.json`
# (`npm run validacion:masiva`).

lib <- file.path(Sys.getenv("HOME"), "Library", "R-anova")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
pdf(NULL)
suppressPackageStartupMessages({
  library(metafor)
  library(meta)
  library(netmeta)
  library(jsonlite)
  library(mada)
})

args <- commandArgs(trailingOnly = TRUE)
modo_masivo <- Sys.getenv("AVANZADA_MASIVA")      # "" = normal; "1" = masivo; "ceros" = solo las familias de celdas en cero (multibrazo con doble cero y eventos = n)
masiva <- nzchar(modo_masivo)
solo_ceros <- modo_masivo == "ceros"
semilla <- if (length(args) >= 1) as.integer(args[1]) else if (masiva) 20261010L else 20261006L
set.seed(semilla)
n_reg <- if (masiva) 0 else 60
n_sroc <- if (solo_ceros) 0 else if (masiva) 400 else 60
n_red <- if (solo_ceros) 120 else if (masiva) 300 else 40

source(file.path("validacion", "ayudas.R"))

nuevo_id <- local({ n <- 0; function(pref) { n <<- n + 1; sprintf("%s%03d", pref, n) } })
disenos <- c("ECA", "Cohortes", "Casos y Controles")

# ================================================================================================
# 1) META-REGRESIÓN de efectos mixtos (REML, pruebas z) — `metafor::rma(yi, vi, mods = …, method = "REML", test = "z")`
# ================================================================================================
# (Los identificadores de moderador van sin tildes —anio, n, diseno— porque Rscript escribe mal los acentos en algunas
# configuraciones regionales; `comparar-avanzada.mts` los traduce a los de aNova.)
# aNova construye la matriz de diseño así: intercepto; moderadores continuos CENTRADOS en su media; categóricos como
# variables ficticias con el primer nivel (orden de código de carácter) como referencia. Se arma igual aquí.
escenarios_reg <- list()
intentos <- 0
while (length(escenarios_reg) < n_reg && intentos < 400) {
  intentos <- intentos + 1
  k <- sample(c(6, 8, 12, 20, 40), 1)
  tipo <- sample(c("anio", "diseno", "anio+diseno", "anio+n"), 1)
  sei <- runif(k, 0.05, 0.5)
  anio <- sample(1995:2022, k, replace = TRUE)
  n <- round(runif(k, 20, 800))
  grupo <- sample(disenos, k, replace = TRUE)
  tau <- sample(c(0, 0.15, 0.4), 1)
  yi <- 0.2 + 0.03 * (anio - 2008) + ifelse(grupo == "ECA", 0, -0.2) + 0.0004 * (n - 400) + rnorm(k, 0, sqrt(sei^2 + tau^2))
  vi <- sei^2

  cols <- list(); nombres <- character(0); mods_aNova <- character(0)
  if (grepl("anio", tipo)) { cols[[length(cols) + 1]] <- anio - mean(anio); mods_aNova <- c(mods_aNova, "anio") }
  if (grepl("n$", tipo)) { cols[[length(cols) + 1]] <- n - mean(n); mods_aNova <- c(mods_aNova, "n") }
  if (grepl("diseno", tipo)) {
    niveles <- sort(unique(grupo), method = "radix")
    if (length(niveles) < 2) next
    for (nv in niveles[-1]) cols[[length(cols) + 1]] <- as.numeric(grupo == nv)
    mods_aNova <- c(mods_aNova, "diseno")
  }
  X <- do.call(cbind, cols)
  if (k <= ncol(X) + 1) next
  fit <- tryCatch(suppressWarnings(rma(yi, vi, mods = X, method = "REML", test = "z", control = ctrl)), error = function(e) NULL)
  if (is.null(fit)) next
  # τ² del modelo sin moderadores: el R² es (τ²₀ − τ²)/τ²₀ y solo tiene sentido si τ²₀ no es un cero numérico.
  fit0 <- tryCatch(suppressWarnings(rma(yi, vi, method = "REML", control = ctrl)), error = function(e) NULL)
  escenarios_reg[[length(escenarios_reg) + 1]] <- list(
    id = nuevo_id("reg"), k = k, tipo = tipo, moderadores = mods_aNova,
    yi = yi, sei = sei, anio = anio, n = n, grupo = grupo,
    ref = list(
      beta = as.numeric(fit$beta), se = as.numeric(fit$se), z = as.numeric(fit$zval), p = as.numeric(fit$pval),
      lb = as.numeric(fit$ci.lb), ub = as.numeric(fit$ci.ub), tau2 = num(fit$tau2),
      tau2_nulo = if (is.null(fit0)) NULL else num(fit0$tau2), qm = num(fit$QM), qmGl = as.integer(round(fit$QMdf[1])), qmP = num(fit$QMp), r2 = num(fit$R2)
    )
  )
}

# ================================================================================================
# 2) SROC / precisión diagnóstica
# ================================================================================================
# aNova corrige con +0,5 las CUATRO celdas de un estudio si alguna es cero (el mismo estudio para sensibilidad, especificidad
# y DOR). Se hace igual aquí y se pasa add = 0 a los cálculos, para comparar los métodos y no la regla de corrección.
#  · DOR y sensibilidad/especificidad combinadas: inverso-varianza sobre logit (`agrupar`, como metafor).
#  · Regresión de Moses-Littenberg: D = a + b·S ponderada por 1/var(D); el error estándar es el de efectos fijos
#    (`rma(method = "FE", mods = ~S)`), sin reescalar por la varianza residual.
#  · IC de Wilson por estudio sobre los conteos ORIGINALES (`prop.test(correct = FALSE)`).
escenarios_sroc <- list()
for (i in seq_len(n_sroc)) {
  k <- sample(c(3, 4, 5, 8, 15, 30), 1)
  nd <- round(runif(k, 15, 400)); nn <- round(runif(k, 15, 600))
  mus <- qlogis(runif(1, 0.6, 0.97)); mup <- qlogis(runif(1, 0.6, 0.97))
  umbral <- rnorm(k, 0, sample(c(0, 0.5, 1), 1)); signo <- sample(c(0, 0.5), 1)
  sens <- plogis(mus + umbral + rnorm(k, 0, 0.3)); spec <- plogis(mup - signo * umbral + rnorm(k, 0, 0.3))
  vp <- rbinom(k, nd, sens); vn <- rbinom(k, nn, spec); fn <- nd - vp; fp <- nn - vn
  cero <- vp == 0 | fn == 0 | fp == 0 | vn == 0
  a <- vp + 0.5 * cero; b <- fn + 0.5 * cero; c <- fp + 0.5 * cero; d <- vn + 0.5 * cero
  D <- log(a * d / (b * c)); seD <- sqrt(1 / a + 1 / b + 1 / c + 1 / d); S <- log(a / b) + log(c / d)
  wil <- function(x, m) { r <- suppressWarnings(prop.test(x, m, correct = FALSE)$conf.int); c(lb = r[1], ub = r[2]) }
  wil_sens <- t(mapply(function(x, m) wil(x, m), vp, vp + fn)); wil_spec <- t(mapply(function(x, m) wil(x, m), vn, fp + vn))
  reg <- NULL
  if (k >= 4) {
    f <- tryCatch(rma(D, seD^2, mods = ~ S, method = "FE"), error = function(e) NULL)
    if (!is.null(f)) reg <- list(a = as.numeric(f$beta[1]), b = as.numeric(f$beta[2]), seA = as.numeric(f$se[1]), seB = as.numeric(f$se[2]),
                                 zB = as.numeric(f$zval[2]), pB = as.numeric(f$pval[2]))
  }
  # mada::madauni como CAJA NEGRA (solo en modo masivo): ver `mada_dor` en ayudas.R
  mada_ref <- if (masiva) mada_dor(vp, fn, fp, vn) else NULL
  escenarios_sroc[[length(escenarios_sroc) + 1]] <- list(
    id = nuevo_id("sroc"), k = k, vp = vp, fp = fp, fn = fn, vn = vn,
    mada = mada_ref,
    ref = list(
      dor = agrupar(D, seD), sens = agrupar(log(a / b), sqrt(1 / a + 1 / b)), esp = agrupar(log(d / c), sqrt(1 / c + 1 / d)),
      wilson_sens = list(lb = as.numeric(wil_sens[, "lb"]), ub = as.numeric(wil_sens[, "ub"])),
      wilson_esp = list(lb = as.numeric(wil_spec[, "lb"]), ub = as.numeric(wil_spec[, "ub"])),
      regresion = reg
    )
  )
}

# ================================================================================================
# 3) META-ANÁLISIS EN RED (Rücker, `netmeta`): efectos de red, heterogeneidad, P-scores, descomposición de Q y
#    node-splitting (back-calculation y SIDDE), con modelo común (fijo) y aleatorio (τ² por DerSimonian-Laird).
# ================================================================================================
# Redes aleatorias conexas de 3 a 6 tratamientos y 8 a 25 estudios, con ~12 % de estudios de tres brazos. Dos medidas:
# diferencia de medias (brazos con media, DE, n) y OR (brazos con eventos, n; SIN celdas en cero, porque aNova corrige cada
# brazo por separado y `pairwise` corrige cada comparación: es otra convención, no un error).
conexa <- function(trts, estudios) {
  vistos <- trts[1]; cambio <- TRUE
  while (cambio) {
    cambio <- FALSE
    for (e in estudios) if (any(e %in% vistos) && !all(e %in% vistos)) { vistos <- union(vistos, e); cambio <- TRUE }
  }
  length(vistos) == length(trts)
}

extraer_modelo <- function(net, trts, sufijo) {
  TE <- net[[paste0("TE.", sufijo)]]; se <- net[[paste0("seTE.", sufijo)]]
  lb <- net[[paste0("lower.", sufijo)]]; ub <- net[[paste0("upper.", sufijo)]]
  n <- length(trts); pares <- list()
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    pares[[length(pares) + 1]] <- list(t1 = trts[i], t2 = trts[j], te = TE[trts[i], trts[j]], se = se[trts[i], trts[j]],
                                       lb = lb[trts[i], trts[j]], ub = ub[trts[i], trts[j]])
  }
  pares
}

extraer_split <- function(ns, sufijo) {
  etiquetas <- ns$comparison
  dir <- ns[[paste0("direct.", sufijo)]]; ind <- ns[[paste0("indirect.", sufijo)]]
  prop <- ns[[paste0("prop.", sufijo)]]; cmpr <- ns[[paste0("compare.", sufijo)]]
  out <- list()
  for (i in seq_along(etiquetas)) {
    partes <- strsplit(etiquetas[i], ":", fixed = TRUE)[[1]]
    out[[length(out) + 1]] <- list(
      t1 = partes[1], t2 = partes[2], prop = num(prop[i]),
      dir_te = num(dir$TE[i]), dir_lb = num(dir$lower[i]), dir_ub = num(dir$upper[i]),
      ind_te = num(ind$TE[i]), ind_lb = num(ind$lower[i]), ind_ub = num(ind$upper[i]),
      z = num(cmpr$statistic[i]), p = num(cmpr$p[i])
    )
  }
  out
}

escenarios_red <- list()
intentos <- 0
# Construye el escenario de red (datos + valores de R). Con la familia «m», `pairwise(allstudies = TRUE)` conserva los estudios cuyos brazos
# tienen TODOS cero eventos (contrastes con TE = 0 y varianza enorme); aNova los descarta por no informar ningún contraste (Handbook 10.4.4.2),
# así que a R se le entregan los mismos datos sin ellos, y a aNova se le entregan completos para comprobar que los descarta y lo declara.
construir_red <- function(d, medida, trts, fam, k, id = NULL) {
  d_r <- d
  if (medida == "OR" && fam == "m") {
    todo_cero <- tapply(d$event, d$study, function(e) all(e == 0))
    d_r <- d[!(d$study %in% names(todo_cero)[todo_cero]), ]
    if (nrow(d_r) == 0) return(NULL)
  }
  p <- tryCatch(suppressWarnings(
    if (medida == "DM") pairwise(treat = treat, n = n, mean = mean, sd = sd, studlab = study, data = d_r, sm = "MD")
    else pairwise(treat = treat, event = event, n = n, studlab = study, data = d_r, sm = "OR", allstudies = (fam == "m"))),
    error = function(e) NULL)
  if (is.null(p)) return(NULL)
  net <- tryCatch(suppressWarnings(netmeta(p, common = TRUE, random = TRUE, method.tau = "DL", reference.group = trts[1])), error = function(e) NULL)
  if (is.null(net)) return(NULL)
  nr <- netrank(net, small.values = "good")
  dd <- tryCatch(suppressWarnings(decomp.design(net)), error = function(e) NULL)
  qd <- if (is.null(dd)) NULL else dd$Q.decomp
  ns_bc <- tryCatch(suppressWarnings(netsplit(net, method = "Back-calculation")), error = function(e) NULL)
  ns_si <- tryCatch(suppressWarnings(netsplit(net, method = "SIDDE")), error = function(e) NULL)
  modelos <- list()
  for (suf in c("common", "random")) {
    modelos[[suf]] <- list(
      pares = extraer_modelo(net, trts, suf),
      pscore = as.numeric(nr[[paste0("ranking.", suf)]][trts]),
      split_bc = if (is.null(ns_bc)) NULL else extraer_split(ns_bc, suf),
      split_sidde = if (is.null(ns_si)) NULL else extraer_split(ns_si, suf)
    )
  }
  filas_json <- lapply(seq_len(nrow(d)), function(i) as.list(d[i, ]))
  list(
    id = if (is.null(id)) nuevo_id(if (fam == "") "red" else paste0("red", fam)) else id, medida = medida, k = k, tratamientos = trts, datos = filas_json,
    ref = list(
      tau2 = num(net$tau2), q = num(net$Q), gl = as.integer(round(net$df.Q)), pq = num(net$pval.Q), i2 = 100 * num(net$I2),
      m = as.integer(round(net$m)), n = as.integer(round(net$n)),
      decomp = if (is.null(qd)) NULL else list(total = as.numeric(qd["Total", ]), dentro = as.numeric(qd["Within designs", ]), entre = as.numeric(qd["Between designs", ])),
      common = modelos$common, random = modelos$random
    )
  )
}

# Recalcula, sobre el JSON de ceros ya generado, solo las redes «m» con algún estudio de todos-los-brazos-en-cero (mismos datos y mismo id).
if (nzchar(Sys.getenv("CEROS_RECALCULAR"))) {
  ruta <- file.path("validacion", "referencias", "referencia-red-ceros.json")
  ref <- fromJSON(ruta, simplifyVector = FALSE)
  cambiadas <- 0
  for (i in seq_along(ref$red)) {
    r <- ref$red[[i]]
    if (!startsWith(r$id, "redm")) next
    d <- do.call(rbind, lapply(r$datos, function(f) data.frame(study = f$study, treat = f$treat, n = f$n, event = f$event, stringsAsFactors = FALSE)))
    todo_cero <- tapply(d$event, d$study, function(e) all(e == 0))
    if (!any(todo_cero)) next
    nuevo <- construir_red(d, "OR", unlist(r$tratamientos), "m", r$k, r$id)
    if (is.null(nuevo)) stop("no se pudo recalcular ", r$id)
    ref$red[[i]] <- nuevo; cambiadas <- cambiadas + 1
    message("recalculada ", r$id)
  }
  write(toJSON(ref, digits = NA, auto_unbox = TRUE, na = "null", null = "null"), ruta)
  cat("Redes recalculadas:", cambiadas, "\n"); quit(save = "no")
}

while (length(escenarios_red) < n_red && intentos < (if (masiva) 4000 else 400)) {
  intentos <- intentos + 1
  medida <- if (solo_ceros) "OR" else if (length(escenarios_red) %% 2 == 0) "DM" else "OR"
  nt <- sample(3:6, 1); trts <- LETTERS[1:nt]
  k <- sample(8:25, 1)
  # Familias con celdas en cero (modo masivo): «z» = brazos con 0 eventos (corrección por estudio, doble cero de 2 brazos descartado);
  # «m» = estudios multibrazo con contrastes de doble cero (solo `pairwise(allstudies = TRUE)` los admite: aNova los conserva con su corrección);
  # «h» = tasas de eventos altas, con brazos donde eventos = n (celda de no-eventos en cero) y dobles «todos con evento» de 2 brazos.
  fam <- if (solo_ceros) c("m", "h")[1 + length(escenarios_red) %% 2] else if (masiva && medida == "OR" && length(escenarios_red) %% 4 == 1) "z" else ""
  con_ceros <- fam != ""
  disenos_red <- lapply(1:k, function(s) sample(trts, if (fam == "m" && nt >= 3 && runif(1) < 0.55) (if (nt >= 4 && runif(1) < 0.25) 4 else 3) else if (masiva && nt >= 4 && runif(1) < 0.06) 4 else if (nt >= 3 && runif(1) < 0.12) 3 else 2))
  if (!conexa(trts, disenos_red) || length(unique(unlist(disenos_red))) < nt) next
  tau <- sample(c(0, 0.15, 0.4), 1)
  efecto <- setNames(c(0, sort(rnorm(nt - 1, 0, 0.6))), trts)
  filas <- list(); ok <- TRUE
  for (s in seq_len(k)) {
    ts <- disenos_red[[s]]; na <- length(ts); base <- rnorm(1, 0, 0.5)
    if (medida == "DM") {
      n <- sample(25:200, na, replace = TRUE)
      filas[[s]] <- data.frame(study = paste0("S", s), treat = ts, n = n, mean = round(10 + base + efecto[ts] + rnorm(na, 0, tau), 3),
                               sd = round(runif(na, 1, 4), 3), stringsAsFactors = FALSE)
    } else {
      n <- sample(40:400, na, replace = TRUE)
      ev <- rbinom(na, n, plogis((if (fam == "m") -4.5 else if (fam == "h") 3.5 else if (con_ceros) -3 else -0.5) + base + efecto[ts] + rnorm(na, 0, tau)))
      if (!con_ceros && any(ev < 1 | ev > n - 1)) { ok <- FALSE; break }
      if (fam %in% c("z", "m") && any(ev > n - 1)) { ok <- FALSE; break }
      if (fam == "m" && na == 2 && all(ev == 0)) { ok <- FALSE; break }    # sin dobles ceros de 2 brazos (allstudies = TRUE los conservaría y aNova los descarta)
      filas[[s]] <- data.frame(study = paste0("S", s), treat = ts, n = n, event = ev, stringsAsFactors = FALSE)
    }
  }
  if (!ok) next
  d <- do.call(rbind, filas)
  esc <- construir_red(d, medida, trts, fam, k)
  if (is.null(esc)) next
  if (masiva) message(sprintf("red %d/%d (%s, %d tratamientos, %d estudios, intento %d)", length(escenarios_red) + 1, n_red, if (con_ceros) paste0("ceros ", fam) else medida, nt, k, intentos))
  escenarios_red[[length(escenarios_red) + 1]] <- esc
}

resultado <- list(
  meta = list(
    fecha = format(Sys.time(), "%Y-%m-%d"), semilla = semilla, R = R.version.string,
    paquetes = list(metafor = as.character(packageVersion("metafor")), meta = as.character(packageVersion("meta")),
                    netmeta = as.character(packageVersion("netmeta")), jsonlite = as.character(packageVersion("jsonlite"))),
    escenarios = list(reg = length(escenarios_reg), sroc = length(escenarios_sroc), red = length(escenarios_red))
  ),
  reg = escenarios_reg, sroc = escenarios_sroc, red = escenarios_red
)

destino <- file.path("validacion", "referencias")
dir.create(destino, recursive = TRUE, showWarnings = FALSE)
archivo <- if (solo_ceros) "referencia-red-ceros.json" else if (masiva) "referencia-red-sroc-masiva.json" else "referencia-avanzada.json"
write(toJSON(resultado, digits = NA, auto_unbox = TRUE, na = "null", null = "null"), file.path(destino, archivo))
cat("Referencia escrita:", file.path(destino, archivo), "\n")
print(resultado$meta$escenarios)
