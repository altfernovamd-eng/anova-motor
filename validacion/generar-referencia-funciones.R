# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Fernando Ramírez Nova
# Referencia de R para las funciones especiales de probabilidad del motor
# (gamma incompleta, beta incompleta, ln Γ): cola superior de χ², cola de dos lados de la t de Student,
# cola de dos lados de la normal y ln Γ. Son funciones base de R (no hace falta ningún paquete).
#
#   Rscript validacion/generar-referencia-funciones.R
#
# Escribe validacion/referencias/referencia-funciones.json, que lee pruebas/test-funciones-especiales.mts.

options(digits = 17, scipen = 999)
pdf(NULL)

fmt <- function(x) sprintf("%.17g", x)

# --- χ², cola superior: gl enteros y semienteros, q desde la cola 1 hasta ~1e-300 ---
gls <- c(0.5, 1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 80, 100, 200, 500, 1000)
colas <- c(0.99, 0.9, 0.5, 0.2, 0.1, 0.05, 0.01, 1e-3, 1e-5, 1e-8, 1e-12, 1e-20, 1e-50, 1e-100, 1e-200, 1e-290)
chi <- list()
for (g in gls) {
  for (p in colas) {
    q <- qchisq(p, g, lower.tail = FALSE)
    if (is.finite(q) && q > 0) chi[[length(chi) + 1]] <- list(gl = g, q = q, p = pchisq(q, g, lower.tail = FALSE))
  }
  # puntos cerca de la frontera de cambio de método (x = a + 1, es decir q = gl + 2)
  for (q in c(g + 1.9, g + 2, g + 2.1, 0.001, 0.1)) {
    chi[[length(chi) + 1]] <- list(gl = g, q = q, p = pchisq(q, g, lower.tail = FALSE))
  }
}

# --- t de Student, cola de AMBOS lados: P(|T| > t) ---
gls_t <- c(0.5, 1, 2, 3, 4, 5, 8, 10, 20, 30, 60, 100, 500, 1000, 10000)
ts <- c(1e-6, 0.001, 0.05, 0.3, 0.7, 1, 1.5, 2, 2.5, 3, 4, 6, 10, 25, 100, 1000)
tt <- list()
for (g in gls_t) {
  for (t in ts) {
    tt[[length(tt) + 1]] <- list(gl = g, t = t, p = 2 * pt(-t, g))
  }
}

# --- cuantil 97,5 % de la t (lo que usan Hartung-Knapp y el intervalo de predicción) ---
gls_q <- c(1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100, 500, 1000)
cuant <- lapply(gls_q, function(g) list(gl = g, q = qt(0.975, g)))

# --- normal, dos lados ---
zs <- c(0, 1e-8, 0.01, 0.1, 0.5, 1, 1.5, 1.96, 2.5, 3, 4, 5, 6, 8, 10, 15, 20, 30)
nor <- lapply(zs, function(z) list(z = z, p = 2 * pnorm(-abs(z))))

# --- cuantiles bilaterales de la t y de la normal (IC y p reportados → error estándar) ---
alfas <- c(0.5, 0.2, 0.1, 0.05, 0.01, 0.001, 1e-6)
gls_c <- c(1, 2, 3, 5, 8, 10, 20, 30, 58, 100, 500, 5000)
cuantt <- list()
for (g in gls_c) for (a in alfas) cuantt[[length(cuantt) + 1]] <- list(gl = g, alfa = a, q = qt(1 - a/2, g))
cuantz <- lapply(alfas, function(a) list(alfa = a, q = qnorm(1 - a/2)))
ppz <- c(0.001, 0.01, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 0.99, 0.999, 0.9999)
cuantn <- lapply(ppz, function(p) list(p = p, q = qnorm(p)))

# --- ln Γ ---
xs <- c(0.01, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4.5, 7, 10, 14.9, 15, 20, 50, 100, 250, 500, 1000, 5000)
lg <- lapply(xs, function(x) list(x = x, v = lgamma(x)))

emitir <- function(lista, campos) {
  paste0("[", paste(vapply(lista, function(e) {
    paste0("{", paste(vapply(campos, function(k) paste0("\"", k, "\":", fmt(e[[k]])), ""), collapse = ","), "}")
  }, ""), collapse = ","), "]")
}

json <- paste0(
  "{\"version\":\"", R.version.string, "\",",
  "\"chi2\":", emitir(chi, c("gl", "q", "p")), ",",
  "\"t\":", emitir(tt, c("gl", "t", "p")), ",",
  "\"t975\":", emitir(cuant, c("gl", "q")), ",",
  "\"normal\":", emitir(nor, c("z", "p")), ",",
  "\"qt2\":", emitir(cuantt, c("gl", "alfa", "q")), ",",
  "\"qz2\":", emitir(cuantz, c("alfa", "q")), ",",
  "\"qn\":", emitir(cuantn, c("p", "q")), ",",
  "\"lgamma\":", emitir(lg, c("x", "v")),
  "}"
)
writeLines(json, "validacion/referencias/referencia-funciones.json")
cat("Escrito: validacion/referencias/referencia-funciones.json\n")
cat(sprintf("chi2=%d t=%d t975=%d normal=%d lgamma=%d\n", length(chi), length(tt), length(cuant), length(nor), length(lg)))
