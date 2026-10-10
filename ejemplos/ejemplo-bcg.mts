// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Ejemplo con datos reales: meta-análisis de 13 ensayos de la vacuna BCG contra la tuberculosis (conjunto clásico `dat.bcg`, paquete `metadat` de R).
 * Lee `bcg.csv` (una fila por estudio: eventos y no-eventos de tratados y controles), calcula el OR de cada estudio y los combina con efectos
 * aleatorios (REML). Para comprobarlo con R: `Rscript ejemplos/comparar-bcg.R`.
 *   npm run ejemplo:bcg
 */
import { readFileSync } from "node:fs";
import { calcularMetaAnalisis, efectoDesde2x2, formatoParaGuardar } from "../src/statistics.ts";
import type { Estudio } from "../src/tipos.ts";

const filas = readFileSync(new URL("./bcg.csv", import.meta.url), "utf8").trim().split("\n").slice(1).map((l) => l.replaceAll('"', "").split(","));

const estudios: Estudio[] = filas.map(([autor, año, tpos, tneg, cpos, cneg], i) => {
  const [eventosTrat, totalTrat] = [Number(tpos), Number(tpos) + Number(tneg)];
  const [eventosCtrl, totalCtrl] = [Number(cpos), Number(cpos) + Number(cneg)];
  const efecto = efectoDesde2x2(eventosTrat, totalTrat, eventosCtrl, totalCtrl, "OR")!;   // OR de cada estudio con su IC 95 %
  return {
    id: i + 1, autor, año, diseño: "ECA", n: String(totalTrat + totalCtrl), medida: "OR",
    // el motor recibe el efecto como texto (así sale de un formulario); formatoParaGuardar conserva 9 cifras significativas
    valor: formatoParaGuardar(efecto.valor), icInferior: formatoParaGuardar(efecto.icInferior), icSuperior: formatoParaGuardar(efecto.icSuperior),
  };
});

const r = calcularMetaAnalisis(estudios, "OR", "REML");
if (!r) throw new Error("No se pudo calcular el meta-análisis");
const f = (x: number, d = 4) => x.toFixed(d);
console.log(`Estudios: ${r.estudios.length}`);
console.log(`OR combinado (efectos aleatorios, REML): ${f(r.aleatorio.estimacion)}  IC 95 % [${f(r.aleatorio.icInferior)}, ${f(r.aleatorio.icSuperior)}]`);
console.log(`Intervalo de predicción: [${f(r.aleatorio.prediccionInferior ?? NaN)}, ${f(r.aleatorio.prediccionSuperior ?? NaN)}]`);
console.log(`Heterogeneidad: Q = ${f(r.heterogeneidad.q, 2)} (gl ${r.heterogeneidad.gl}), τ² = ${f(r.heterogeneidad.tau2)}, I² = ${f(r.heterogeneidad.i2, 1)} %`);
console.log("Nota: este I² es el clásico de Higgins (Q − gl)/Q; metafor muestra otro basado en τ² (ver docs/VALIDACION.md).");
