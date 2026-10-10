# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Fernando Ramírez Nova
# Referencia de R para las conversiones de datos (Cochrane Handbook cap. 6 y 10.6).
#
# Casi todas son identidades algebraicas, así que se prueban contra DATOS SIMULADOS CON LA VERDAD CONOCIDA: se genera una muestra,
# R calcula con t.test su IC, su estadístico t y su p, y se comprueba que la conversión de aNova recupera la DE verdadera
# (la muestral). La combinación de grupos y la DE del cambio se comparan con el cálculo directo sobre los datos completos.
# Las aproximaciones de Wan et al. (2014) se comparan con `meta::metamean` usado como CAJA NEGRA (solo se guardan sus resultados).
#
#   Rscript validacion/generar-referencia-conversiones.R
#
# Escribe validacion/referencias/referencia-conversiones.json, que lee pruebas/test-conversiones.mts.

suppressMessages(library(meta, lib.loc = file.path(Sys.getenv("HOME"), "Library/R-anova")))
options(digits = 17, scipen = 999)
pdf(NULL)
set.seed(20261006)
fmt <- function(x) sprintf("%.17g", x)
niveles <- c(0.90, 0.95, 0.99)

# 1) IC de una media → DE
ic_media <- list()
for (i in 1:40) {
  n <- sample(c(4:30, 45, 59, 100, 250), 1)
  x <- rnorm(n, runif(1, -20, 60), runif(1, 0.5, 40))
  for (nv in niveles) {
    ci <- t.test(x, conf.level = nv)$conf.int
    se <- sd(x) / sqrt(n); z <- qnorm(1 - (1 - nv) / 2)
    ic_media[[length(ic_media) + 1]] <- list(n = n, nivel = nv, inf = ci[1], sup = ci[2], sd = sd(x),
                                             infZ = mean(x) - z * se, supZ = mean(x) + z * se, ee = se)
  }
}

# 2) Diferencia de medias: IC, t, p → EE y DE intragrupo (t.test con varianzas iguales, como supone el Handbook)
dif <- list()
for (i in 1:40) {
  n1 <- sample(c(4:40, 80, 150), 1); n2 <- sample(c(4:40, 80, 150), 1)
  x <- rnorm(n1, runif(1, 0, 30), runif(1, 1, 20)); y <- rnorm(n2, runif(1, 0, 30), runif(1, 1, 20))
  sp <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
  for (nv in niveles) {
    tt <- t.test(x, y, var.equal = TRUE, conf.level = nv)
    ee <- as.numeric(tt$stderr)
    dif[[length(dif) + 1]] <- list(n1 = n1, n2 = n2, nivel = nv, md = mean(x) - mean(y), inf = tt$conf.int[1], sup = tt$conf.int[2],
                                   t = as.numeric(tt$statistic), p = tt$p.value, ee = ee, sp = sp)
  }
}

# 3) Combinar dos grupos: contra el cálculo directo sobre los datos juntos
comb <- list()
for (i in 1:40) {
  n1 <- sample(2:120, 1); n2 <- sample(2:120, 1)
  x <- rnorm(n1, runif(1, -10, 50), runif(1, 0.5, 15)); y <- rnorm(n2, runif(1, -10, 50), runif(1, 0.5, 15))
  z <- c(x, y)
  comb[[length(comb) + 1]] <- list(n1 = n1, m1 = mean(x), sd1 = sd(x), n2 = n2, m2 = mean(y), sd2 = sd(y), n = length(z), m = mean(z), sd = sd(z))
}

# 4) DE del cambio desde la línea base: identidad exacta con la correlación de la muestra
cam <- list()
for (i in 1:40) {
  n <- sample(5:200, 1); rho <- runif(1, 0, 0.95)
  b <- rnorm(n, 50, runif(1, 3, 15)); f <- b * rho + rnorm(n, 0, runif(1, 3, 15)) + 5
  cam[[length(cam) + 1]] <- list(sdBase = sd(b), sdFinal = sd(f), sdCambio = sd(f - b), cor = cor(b, f))
}

# 5) Wan et al. (2014): meta::metamean como caja negra. Se guardan la media y la DE que devuelve para cada escenario.
wan <- list()
ns <- c(5, 9, 13, 21, 30, 50, 100, 201, 300, 500, 1000)
for (n in ns) for (rep in 1:4) {
  x <- rnorm(n, runif(1, 0, 40), runif(1, 2, 20)); s <- quantile(x, c(0, 0.25, 0.5, 0.75, 1), type = 7)
  a <- s[[1]]; q1 <- s[[2]]; m <- s[[3]]; q3 <- s[[4]]; b <- s[[5]]
  r1 <- metamean(n = n, median = m, min = a, max = b, method.mean = "Wan", method.sd = "Wan", sm = "MRAW", studlab = "x")
  r2 <- metamean(n = n, median = m, q1 = q1, q3 = q3, method.mean = "Wan", method.sd = "Wan", sm = "MRAW", studlab = "x")
  r3 <- metamean(n = n, median = m, q1 = q1, q3 = q3, min = a, max = b, method.mean = "Wan", method.sd = "Wan", sm = "MRAW", studlab = "x")
  wan[[length(wan) + 1]] <- list(n = n, a = a, q1 = q1, m = m, q3 = q3, b = b,
                                 media1 = r1$mean, sd1 = r1$sd, media2 = r2$mean, sd2 = r2$sd, media3 = r3$mean, sd3 = r3$sd,
                                 sdVerdadera = sd(x))
}

# 6) OR ↔ DME (Handbook 10.6) y EE desde IC/P de razones y de medidas absolutas
orr <- lapply(c(0.2, 0.5, 0.8, 1, 1.5, 2.7, 6), function(or) list(or = or, ee = 0.3, dme = log(or) * sqrt(3) / pi, eeDme = 0.3 * sqrt(3) / pi))
ee_ic <- lapply(c(0.9, 0.95, 0.99), function(nv) list(nivel = nv, inf = 0.62, sup = 1.91, eeAbs = (1.91 - 0.62) / (2 * qnorm(1 - (1 - nv) / 2)),
                                                       eeLog = (log(1.91) - log(0.62)) / (2 * qnorm(1 - (1 - nv) / 2))))
ee_p <- lapply(c(0.5, 0.2, 0.05, 0.01, 0.0001), function(p) list(p = p, est = 0.9, razon = 1.8, eeAbs = 0.9 / qnorm(1 - p / 2), eeLog = log(1.8) / qnorm(1 - p / 2)))

emitir <- function(lista) {
  paste0("[", paste(vapply(lista, function(e) {
    paste0("{", paste(vapply(names(e), function(k) paste0("\"", k, "\":", fmt(e[[k]])), ""), collapse = ","), "}")
  }, ""), collapse = ","), "]")
}
json <- paste0("{\"meta\":\"", as.character(packageVersion("meta")), "\",",
               "\"icMedia\":", emitir(ic_media), ",\"dif\":", emitir(dif), ",\"comb\":", emitir(comb), ",\"cambio\":", emitir(cam),
               ",\"wan\":", emitir(wan), ",\"orDme\":", emitir(orr), ",\"eeIc\":", emitir(ee_ic), ",\"eeP\":", emitir(ee_p), "}")
writeLines(json, "validacion/referencias/referencia-conversiones.json")
cat(sprintf("Escrito: icMedia=%d dif=%d comb=%d cambio=%d wan=%d\n", length(ic_media), length(dif), length(comb), length(cam), length(wan)))
