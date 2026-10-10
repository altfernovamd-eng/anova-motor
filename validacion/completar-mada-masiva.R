# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Fernando Ramírez Nova
# Añade (o rehace) el bloque `mada` de cada conjunto de SROC de `referencias/referencia-red-sroc-masiva.json` sin regenerar las 300 redes
# (que tardan ~50 min). Usa la misma función `mada_dor` (ayudas.R) que el generador. Rscript validacion/completar-mada-masiva.R
lib <- file.path(Sys.getenv("HOME"), "Library", "R-anova")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({ library(mada); library(jsonlite) })
source(file.path("validacion", "ayudas.R"))
archivo <- file.path("validacion", "referencias", "referencia-red-sroc-masiva.json")
ref <- fromJSON(archivo, simplifyVector = FALSE)
for (i in seq_along(ref$sroc)) {
  e <- ref$sroc[[i]]
  ref$sroc[[i]]$mada <- mada_dor(unlist(e$vp), unlist(e$fn), unlist(e$fp), unlist(e$vn))
}
write(toJSON(ref, digits = NA, auto_unbox = TRUE, na = "null", null = "null"), archivo)
cat("mada añadido a", length(ref$sroc), "conjuntos de SROC\n")
