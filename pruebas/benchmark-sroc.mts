// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Prueba de regresión: valida `calcularSROC` (src/lib/meta-analysis/statistics.ts).
 *
 * Dataset: AuditC (14 estudios) — Kriston, Hölzel, Weiser, Berner & Härter (2008), incluido en el
 * paquete `mada` de R (`data("AuditC")`, el estándar de referencia para meta-análisis de precisión
 * diagnóstica). Obtenido parseando el .rda oficial del paquete (vía el paquete Python `rdata`),
 * TP/FP/FN/TN exactos, incluye 2 estudios con FN=0 (celda cero) — pone a prueba la corrección de
 * continuidad.
 *
 * **Por qué NO se compara literalmente contra `summary(madauni(AuditC))`** (que publica DOR=26.337,
 * IC[17.971,38.596], τ²=0.311): se investigó el código fuente de `mada` (R/madauni.R) y usa una
 * convención de corrección de continuidad DISTINTA a la del resto de este motor — aplica +0.5 a las
 * 4 celdas de TODOS los estudios siempre, no solo a los que tienen una celda en cero. aNova sigue,
 * aquí como en `efectoDesde2x2` y todo el motor, la recomendación del Cochrane Handbook (6.4.3.1):
 * corregir SOLO el estudio que de verdad tiene una celda cero — es la práctica moderna recomendada
 * (menos sesgo introducido en los estudios que no lo necesitan) y la que ya se valida en el resto de
 * `VALIDACION_ESTADISTICA.md`. Cambiar de convención solo para SROC rompería esa consistencia interna.
 * (Se confirmó, reproduciendo el algoritmo de mada con la corrección "a todos" en Python, que con ESA
 * convención sí se obtiene DOR=26.337/IC[17.971,38.596]/τ²=0.311 exacto — la fórmula de pooling en sí
 * es idéntica a mada; la única diferencia real es el alcance de la corrección de continuidad.)
 *
 * Por eso el DOR/τ²/IC combinados se validan aquí contra una **implementación independiente en
 * Python/numpy** que replica la misma convención (corrección solo en celdas cero) — confirma que el
 * código implementa correctamente las fórmulas publicadas (DerSimonian & Laird 1986; Glas et al. 2003
 * para el DOR), no que coincide con un paquete que corrige distinto.
 *
 * La regresión de Moses-Littenberg (Moses, Shapiro & Littenberg 1993; D~S por WLS) también se valida
 * contra Python/numpy (mínimos cuadrados ponderados por ecuaciones normales) sobre el mismo dataset.
 * El intervalo de Wilson se valida contra el ejemplo de texto clásico (Newcombe 1998): x=81, n=263 →
 * IC 95% [0.2553, 0.3662].
 *
 * Cómo correrlo:   npm run benchmark:sroc      (o)   node pruebas/benchmark-sroc.mts
 * Sale con código 1 si alguna cifra deja de coincidir.
 */
import { calcularSROC } from "../src/statistics.ts";
import type { Estudio } from "../src/tipos.ts";

// TP, FP, FN, TN — exactos, extraídos del AuditC.rda oficial de `mada` (data("AuditC")).
const AUDITC: [number, number, number, number][] = [
  [47, 101, 9, 738],
  [126, 272, 51, 1543],
  [19, 12, 10, 192],
  [36, 78, 3, 276],
  [130, 211, 19, 959],
  [84, 68, 2, 89],
  [68, 112, 0, 423],
  [752, 3226, 0, 2977],
  [59, 55, 5, 136],
  [142, 571, 50, 2788],
  [137, 107, 24, 358],
  [57, 103, 3, 437],
  [34, 21, 1, 56],
  [152, 88, 51, 264],
];

const estudios: Estudio[] = AUDITC.map(([tp, fp, fn, tn], i) => ({
  id: i + 1,
  autor: `Study ${i + 1}`,
  año: "",
  diseño: "Transversal",
  n: "",
  medida: "OR",
  valor: "",
  icInferior: "",
  icSuperior: "",
  datosDiagnosticos: { vp: String(tp), fp: String(fp), fn: String(fn), vn: String(tn) },
}));

const r = calcularSROC(estudios, "aleatorio", "DL");
if (!r) {
  console.error("FALLO: calcularSROC devolvió null");
  process.exit(1);
}

let fallos = 0;
let total = 0;
const chk = (nombre: string, anova: number, ref: number, dp: number) => {
  total++;
  const a = anova.toFixed(dp);
  const b = ref.toFixed(dp);
  const ok = a === b;
  if (!ok) fallos++;
  console.log(`${nombre.padEnd(32)} ${a.padEnd(12)} ${b.padEnd(12)} ${ok ? "✓" : "✗ DIFIERE"}`);
};

console.log(
  `\naNova (calcularSROC)  vs  Python/numpy independiente (misma convención: corrige solo celdas cero)  —  AuditC, k=${r.puntos.length}\n`,
);
console.log(`${"Cantidad".padEnd(32)} ${"aNova".padEnd(12)} ${"referencia".padEnd(12)} Estado`);
console.log("-".repeat(64));
chk("DOR combinado (aleatorio)", r.dor.estimacion, 27.943519, 3);
chk("DOR IC 95% inferior", r.dor.icInferior, 18.738169, 3);
chk("DOR IC 95% superior", r.dor.icSuperior, 41.671106, 3);
chk("τ² (DerSimonian-Laird)", r.dor.tau2, 0.338577, 3);
chk("Q de Cochran (pesos fijos)", r.dor.q, 54.583225, 3);
chk("I² (%)", r.dor.i2, 76.183159, 3);
console.log("-".repeat(64));

console.log(`\nRegresión de Moses-Littenberg  vs  WLS independiente (numpy, ecuaciones normales)\n`);
console.log(`${"Cantidad".padEnd(32)} ${"aNova".padEnd(12)} ${"numpy".padEnd(12)} Estado`);
console.log("-".repeat(64));
if (!r.regresion) {
  console.error("FALLO: la regresión de Moses-Littenberg no se calculó (se esperaban ≥4 estudios)");
  process.exit(1);
}
chk("Intercepto (a)", r.regresion.a, 2.915913, 4);
chk("SE(a)", r.regresion.seA, 0.083171, 4);
chk("Pendiente (b)", r.regresion.b, 0.386329, 4);
chk("SE(b)", r.regresion.seB, 0.084632, 4);
chk("z(b)", r.regresion.zB, 4.5648, 3);
console.log("-".repeat(64));

console.log(`\nIntervalo de Wilson  vs  ejemplo de texto (Newcombe 1998): x=81, n=263\n`);
const wilsonEstudio: Estudio = {
  id: 999,
  autor: "Wilson",
  año: "",
  diseño: "Transversal",
  n: "",
  medida: "OR",
  valor: "",
  icInferior: "",
  icSuperior: "",
  // vp+fn = 263 (todos "enfermos"), vp = 81 → sensibilidad = 81/263; fp/vn cualquier valor no-cero para
  // que el punto sea válido (no afecta el cálculo de sensibilidad, que es independiente de fp/vn).
  datosDiagnosticos: { vp: "81", fp: "10", fn: "182", vn: "50" },
};
// calcularSROC exige k≥3 estudios; se arma junto a 2 estudios neutros solo para poder leer el punto.
const relleno: Estudio[] = [1, 2].map((n) => ({
  id: n,
  autor: `Relleno ${n}`,
  año: "",
  diseño: "Transversal",
  n: "",
  medida: "OR",
  valor: "",
  icInferior: "",
  icSuperior: "",
  datosDiagnosticos: { vp: "20", fp: "5", fn: "5", vn: "20" },
}));
const rw = calcularSROC([wilsonEstudio, ...relleno], "aleatorio", "DL");
if (!rw) {
  console.error("FALLO: calcularSROC devolvió null (caso Wilson)");
  process.exit(1);
}
const punto = rw.puntos.find((p) => p.id === 999)!;
console.log(`${"Cantidad".padEnd(32)} ${"aNova".padEnd(12)} ${"referencia".padEnd(12)} Estado`);
console.log("-".repeat(64));
chk("Wilson IC inferior", punto.sensIcInferior, 0.2553, 4);
chk("Wilson IC superior", punto.sensIcSuperior, 0.3662, 4);
console.log("-".repeat(64));

if (fallos === 0) {
  console.log(`\n✓ PASA: las ${total} cantidades coinciden con las referencias a los decimales indicados.\n`);
  process.exit(0);
} else {
  console.error(`\n✗ FALLA: ${fallos} de ${total} cantidad(es) ya no coinciden. Revisa calcularSROC.\n`);
  process.exit(1);
}
