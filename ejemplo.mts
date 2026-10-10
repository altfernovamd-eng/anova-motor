// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
import { calcularMetaAnalisis } from "./src/statistics.ts";
import type { Estudio } from "./src/tipos.ts";

const estudio = (id: number, valor: string, icInferior: string, icSuperior: string): Estudio =>
  ({ id, autor: `Estudio ${id}`, año: "2020", diseño: "ECA", n: "200", medida: "OR", valor, icInferior, icSuperior });

const resultado = calcularMetaAnalisis(
  [estudio(1, "0.80", "0.50", "1.20"), estudio(2, "0.65", "0.40", "1.05"), estudio(3, "0.90", "0.60", "1.35")],
  "OR",
  "REML",
);
console.log(resultado?.aleatorio);      // OR combinado y su IC 95 %
console.log(resultado?.heterogeneidad); // Q, I², τ², ...
