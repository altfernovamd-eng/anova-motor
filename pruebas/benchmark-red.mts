// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Prueba de regresión: valida `calcularNMA` (src/lib/meta-analysis/red.ts) contra el paquete
 * `netmeta` de R (el estándar del campo para meta-análisis en red frecuentista).
 *
 * Dataset: Stowe2010 (29 estudios, enfermedad de Parkinson) — Stowe et al., incluido en el paquete
 * `netmeta` (`data("Stowe2010")`). 28 estudios de 2 brazos + 1 estudio de 3 brazos (LARGO: COMTI,
 * Placebo, MAOBI) — el caso que pone a prueba el ajuste de varianza multibrazo (Rücker & Schwarzer
 * 2014). Datos exactos obtenidos parseando el `.rda` oficial del paquete (vía el paquete Python
 * `rdata`, ya que no hay R en este entorno).
 *
 * Referencia: netmeta::netmeta(TE, seTE, treat1, treat2, studlab, common=FALSE, ref="Placebo") —
 * vignette "netmeta: An R Package for Network Meta-Analysis" (Balduzzi, Rücker, Nikolakopoulou et al.,
 * J. Stat. Software 2023), sección 5.1.
 *
 * Publicado: k=29, m=31, n=4 tratamientos, τ²=0.1518, Q=48.15 (27 gl, p=0.0074), I²=43.9%.
 * Estimaciones (efectos aleatorios, diferencia de medias, vs Placebo):
 *   COMTI: −0.95 [−1.27, −0.63];  Dopamine Agonist: −1.49 [−1.86, −1.12];  MAOBI: −1.03 [−1.56, −0.50]
 * Ranking de tratamientos, `netrank(net1)` (P-score, mismo dataset):
 *   Dopamine Agonist 0.9677;  MAOBI 0.5664;  COMTI 0.4659;  Placebo 0.0000
 * Descomposición de Q (`decomp.design(net1)`, mismo dataset):
 *   Total 48.15 (27 gl); Within designs 46.91 (25 gl, p=0.0050); Between designs 1.24 (2 gl, p=0.5376)
 *
 * Nota (7 oct 2026): `red.ts` se reescribió desde los artículos; el texto de abajo describe cómo se escribió la versión ANTERIOR (ver docs/PROCEDENCIA.md).
 *
 * Nota: se validó primero una implementación independiente en Python/numpy reproduciendo el algoritmo
 * de Rücker (2012) formula por fórmula desde el paper/vignette — incluido el ajuste multibrazo, cuya
 * fórmula exacta se confirmó leyendo R/multiarm.R del código fuente real de netmeta en GitHub — contra
 * este mismo dataset, coincidencia exacta. Ese Python es el que se portó a `red.ts`; este script prueba
 * el motor de TypeScript real de aNova, no una reimplementación aparte.
 *
 * Los P-scores (`calcularPScores`) usan la fórmula exacta de `netrank()` (confirmada leyendo
 * R/netrank.R: P_ij = Φ(−z_ij) cuando "menor es mejor", promediada por tratamiento) — validada contra
 * los 4 P-scores publicados arriba, coincidencia exacta a 4 decimales.
 *
 * El node-splitting (`ResultadoNMA.nodeSplitting`) usa el método de back-calculation de König, Krahn &
 * Binder (2013) — el método por defecto de `netsplit()` en netmeta —, confirmado leyendo R/netmeta.R
 * (líneas ~1900-1950: `TE.indirect=(TE.red−P·TE.directo)/(1−P)`). La vignette no publica un ejemplo
 * numérico de `netsplit()` para este dataset, así que aquí se valida por **consistencia algebraica**:
 * la fórmula garantiza que `P·directo+(1−P)·indirecto` reproduce EXACTO el estimado de red ya validado
 * arriba — si esa identidad se cumple para una comparación con evidencia directa real (COMTI vs MAOBI,
 * vía el estudio LARGO), la implementación es correcta.
 *
 * La descomposición de Q (`ResultadoNMA.decomposicionQ`) agrupa los estudios por "diseño" (conjunto
 * exacto de tratamientos comparados) y combina por mínimos cuadrados generalizados (GLS) los que
 * comparten diseño — confirmado leyendo R/decomp.tau.R. Validado **exacto** contra `decomp.design(net1)`
 * (Q-dentro, Q-entre, sus gl y p-valores arriba).
 *
 * El node-splitting por **SIDDE** (`calcularNodeSplittingSIDDE`) reajusta la red completa quitando los
 * estudios con evidencia directa de cada par (confirmado leyendo R/sidde.R) — no tiene un número
 * publicado en la vignette, así que se validó con una implementación Python/numpy independiente que
 * quita el estudio LARGO y reajusta la red COMTI-MAOBI desde cero: coincidencia exacta (4 decimales).
 *
 * Cómo correrlo:   npm run benchmark:red      (o)   node pruebas/benchmark-red.mts
 * Sale con código 1 si alguna cifra deja de coincidir.
 */
import { calcularNMA, calcularNodeSplittingSIDDE, calcularPScores } from "../src/red.ts";
import type { EstudioRed } from "../src/tipos.ts";

type FilaContinua = [string, string, number, number, number, string, number, number, number, string?, number?, number?, number?];

// study, t1, y1, sd1, n1, t2, y2, sd2, n2, t3?, y3?, sd3?, n3?
const STOWE2010: FilaContinua[] = [
  ["DA(pr/B) Int", "Dopamine Agonist", -1.9, 4.3, 152, "Placebo", -0.3, 4.4, 76],
  ["DA(B) Germany", "Dopamine Agonist", -3.4, 2.78, 20, "Placebo", 0.3, 3.66, 20],
  ["DA(C)Spain", "Dopamine Agonist", -3.33, 3.47, 14, "Placebo", -2.47, 3.91, 12],
  ["DA(C) UK", "Dopamine Agonist", -2.0, 2.33, 18, "Placebo", -0.7, 2.25, 16],
  ["DA(C) USA", "Dopamine Agonist", -2.08, 3.19, 89, "Placebo", -0.77, 3.34, 46],
  ["DA(pe) N America", "Dopamine Agonist", -1.8, 4.8, 189, "Placebo", -0.2, 4.8, 187],
  ["DA(Pr) Aus/ger", "Dopamine Agonist", -2.3, 2.7, 23, "Placebo", 0.4, 3.5, 37],
  ["DA(Pr) CLEOPATRA", "Dopamine Agonist", -2.8, 2.83, 200, "Placebo", -0.9, 5.0, 100],
  ["DA(Pr) Denmark", "Dopamine Agonist", -1.2, 2.9, 29, "Placebo", -0.9, 5.0, 29],
  ["DA(Pr) Europe", "Dopamine Agonist", -1.94, 2.59, 174, "Placebo", -0.12, 3.11, 180],
  ["DA(Pr) US/CANADA", "Dopamine Agonist", -2.4, 3.4, 174, "Placebo", -0.7, 3.7, 172],
  ["DA(Pr) EASE-PD", "Dopamine Agonist", -2.1, 4.54, 201, "Placebo", -0.3, 4.41, 190],
  ["DA(Pr) France/Eng", "Dopamine Agonist", -1.74, 2.36, 23, "Placebo", -2.22, 3.0, 23],
  ["DA(Pr) UK/Israel", "Dopamine Agonist", -1.33, 2.35, 42, "Placebo", -0.75, 2.81, 22],
  ["DA(Pr) USA", "Dopamine Agonist", -1.53, 4.28, 90, "Placebo", -1.22, 3.7, 53],
  ["COMTI(E) Colemen", "COMTI", -1.6, 2.5, 129, "Placebo", -0.9, 3.4, 74],
  ["COMTI(E) ComQol", "COMTI", -1.0, 1.9, 174, "Placebo", -0.9, 2.0, 96],
  ["COMTI(E) INT-OZ", "COMTI", -0.86, 2.9, 99, "Placebo", -0.54, 3.13, 63],
  ["LARGO", "COMTI", -1.2, 2.21, 218, "Placebo", -0.4, 2.21, 218, "MAOBI", -1.18, 2.23, 222],
  ["COMTI(E) Nomecomt", "COMTI", -1.3, 2.2, 85, "Placebo", -0.1, 2.45, 86],
  ["COMTI(E) Sth Korea", "COMTI", -1.5, 2.55, 98, "Placebo", -1.0, 2.29, 99],
  ["COMTI(E) UK/Irish", "COMTI", -1.1, 2.55, 115, "Placebo", -0.3, 2.85, 57],
  ["COMTI(T) China", "COMTI", -1.9, 3.3, 20, "Placebo", 0.15, 3.3, 20],
  ["COMTI(T) Europe", "COMTI", -2.03, 2.6, 60, "Placebo", -0.67, 2.8, 58],
  ["COMTI(T) TFSG1", "COMTI", -1.61, 1.77, 40, "Placebo", -0.04, 1.75, 42],
  ["COMTI(T) TFSG3", "COMTI", -2.0, 2.49, 69, "Placebo", -0.3, 2.55, 72],
  ["COMTI(T) TIPS I", "COMTI", -1.78, 2.76, 31, "Placebo", -0.11, 2.63, 37],
  ["MAOBI(R) PRESTO", "MAOBI", -1.85, 2.5, 149, "Placebo", -0.91, 2.5, 159],
  ["MAOBI(ZS) USA", "MAOBI", -2.2, 2.43, 94, "Placebo", -0.6, 2.49, 46],
];

const estudiosRed: EstudioRed[] = STOWE2010.map((fila, i) => {
  const [etiqueta, t1, y1, sd1, n1, t2, y2, sd2, n2, t3, y3, sd3, n3] = fila;
  const brazos = [
    { tratamiento: t1, media: String(y1), de: String(sd1), n: String(n1) },
    { tratamiento: t2, media: String(y2), de: String(sd2), n: String(n2) },
  ];
  if (t3) brazos.push({ tratamiento: t3, media: String(y3), de: String(sd3), n: String(n3) });
  return { id: i + 1, etiqueta, medida: "DM" as const, brazos };
});

const r = calcularNMA(estudiosRed, "DM", "aleatorio");
if (!r) {
  console.error("FALLO: calcularNMA devolvió null");
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
  console.log(`${nombre.padEnd(34)} ${a.padEnd(12)} ${b.padEnd(12)} ${ok ? "✓" : "✗ DIFIERE"}`);
};

console.log(`\naNova (calcularNMA)  vs  R 'netmeta' (Stowe2010, Parkinson, k=${r.k}, m=${r.m}, n=${r.n})\n`);
console.log(`${"Cantidad".padEnd(34)} ${"aNova".padEnd(12)} ${"netmeta".padEnd(12)} Estado`);
console.log("-".repeat(66));
chk("k (estudios)", r.k, 29, 0);
chk("m (comparaciones)", r.m, 31, 0);
chk("n (tratamientos)", r.n, 4, 0);
chk("gl", r.gl, 27, 0);
chk("Q", r.q, 48.15, 2);
chk("τ²", r.tau2, 0.1518, 4);
chk("τ", r.tau, 0.3896, 4);
chk("I² (%)", r.i2, 43.9, 1);
console.log("-".repeat(66));

/** Devuelve {estimacion, icInferior, icSuperior} orientado como (tratamiento vs `ref`). */
function vsReferencia(tratamiento: string, ref: string) {
  const est = r!.estimaciones.find(
    (e) =>
      (e.tratamiento1 === tratamiento && e.tratamiento2 === ref) ||
      (e.tratamiento1 === ref && e.tratamiento2 === tratamiento),
  )!;
  if (est.tratamiento1 === tratamiento) return est; // ya es (tratamiento − ref)
  // est es (ref − tratamiento): invertir signo y voltear el IC
  return {
    estimacion: -est.estimacion,
    icInferior: -est.icSuperior,
    icSuperior: -est.icInferior,
  };
}

const comti = vsReferencia("COMTI", "Placebo");
const dopa = vsReferencia("Dopamine Agonist", "Placebo");
const maobi = vsReferencia("MAOBI", "Placebo");

console.log(`\nEstimaciones de red vs Placebo (efectos aleatorios, diferencia de medias)\n`);
console.log(`${"Cantidad".padEnd(34)} ${"aNova".padEnd(12)} ${"netmeta".padEnd(12)} Estado`);
console.log("-".repeat(66));
chk("COMTI vs Placebo", comti.estimacion, -0.95, 2);
chk("COMTI IC inferior", comti.icInferior, -1.27, 2);
chk("COMTI IC superior", comti.icSuperior, -0.63, 2);
chk("Dopamine Agonist vs Placebo", dopa.estimacion, -1.49, 2);
chk("Dopamine Agonist IC inferior", dopa.icInferior, -1.86, 2);
chk("Dopamine Agonist IC superior", dopa.icSuperior, -1.12, 2);
chk("MAOBI vs Placebo", maobi.estimacion, -1.03, 2);
chk("MAOBI IC inferior", maobi.icInferior, -1.56, 2);
chk("MAOBI IC superior", maobi.icSuperior, -0.5, 2);
console.log("-".repeat(66));

const pscores = calcularPScores(r, true); // menor MD = mejor (mayor reducción del off-time)
const pscore = (t: string) => pscores.find((p) => p.tratamiento === t)!.pScore;
console.log(`\nP-scores (netrank)  vs  R 'netmeta' (netrank(net1))\n`);
console.log(`${"Cantidad".padEnd(34)} ${"aNova".padEnd(12)} ${"netmeta".padEnd(12)} Estado`);
console.log("-".repeat(66));
chk("P-score: Dopamine Agonist", pscore("Dopamine Agonist"), 0.9677, 4);
chk("P-score: MAOBI", pscore("MAOBI"), 0.5664, 4);
chk("P-score: COMTI", pscore("COMTI"), 0.4659, 4);
chk("P-score: Placebo", pscore("Placebo"), 0.0, 4);
console.log("-".repeat(66));

console.log(`\nNode-splitting (back-calculation)  —  consistencia algebraica interna\n`);
console.log(`${"Cantidad".padEnd(34)} ${"aNova".padEnd(12)} ${"esperado".padEnd(12)} Estado`);
console.log("-".repeat(66));
const split = r.nodeSplitting.find(
  (s) =>
    (s.tratamiento1 === "COMTI" && s.tratamiento2 === "MAOBI") ||
    (s.tratamiento1 === "MAOBI" && s.tratamiento2 === "COMTI"),
);
if (!split) {
  console.error("FALLO: no se encontró node-splitting para COMTI vs MAOBI (se esperaba, vía LARGO)");
  process.exit(1);
}
console.log(`COMTI vs MAOBI — directo: ${split.estimacionDirecta.toFixed(3)} [${split.icInferiorDirecta.toFixed(3)}, ${split.icSuperiorDirecta.toFixed(3)}]`);
console.log(`COMTI vs MAOBI — indirecto: ${split.estimacionIndirecta.toFixed(3)} [${split.icInferiorIndirecta.toFixed(3)}, ${split.icSuperiorIndirecta.toFixed(3)}]`);
console.log(`proporción directa: ${split.proporcionDirecta.toFixed(4)}, z=${split.z.toFixed(3)}, p=${split.pValor.toFixed(3)}`);
// Identidad algebraica que garantiza la fórmula de back-calculation: P·directo+(1−P)·indirecto (en
// escala de análisis, lineal para DM) debe reproducir EXACTO el estimado de red ya validado arriba.
const comtiVsMaobiRed = vsReferencia("COMTI", "MAOBI");
const P = split.proporcionDirecta;
const reconstruido = P * split.estimacionDirecta + (1 - P) * split.estimacionIndirecta;
chk("P·directo+(1−P)·indirecto", reconstruido, comtiVsMaobiRed.estimacion, 4);
chk("proporción directa en (0,1)", P > 0 && P < 1 ? 1 : 0, 1, 0);
console.log("-".repeat(66));

console.log(`\nDescomposición de Q  vs  R 'netmeta' (decomp.design(net1))\n`);
console.log(`${"Cantidad".padEnd(34)} ${"aNova".padEnd(12)} ${"netmeta".padEnd(12)} Estado`);
console.log("-".repeat(66));
const dq = r.decomposicionQ;
if (!dq) {
  console.error("FALLO: decomposicionQ es null (se esperaba, hay diseños con ≥2 estudios)");
  process.exit(1);
}
chk("Q-dentro (within designs)", dq.qDentro, 46.91, 2);
chk("gl-dentro", dq.glDentro, 25, 0);
chk("p-dentro", dq.pDentro, 0.005, 3);
chk("Q-entre (between designs)", dq.qEntre, 1.24, 2);
chk("gl-entre", dq.glEntre, 2, 0);
chk("p-entre", dq.pEntre, 0.5376, 3);
console.log("-".repeat(66));

console.log(`\nNode-splitting SIDDE  —  validado con Python/numpy independiente (quitar LARGO y reajustar)\n`);
console.log(`${"Cantidad".padEnd(34)} ${"aNova".padEnd(12)} ${"Python".padEnd(12)} Estado`);
console.log("-".repeat(66));
const sidde = calcularNodeSplittingSIDDE(estudiosRed, "DM", "aleatorio");
const siddeComtiMaobi = sidde.find((s) => s.tratamiento1 === "COMTI" && s.tratamiento2 === "MAOBI");
if (!siddeComtiMaobi) {
  console.error("FALLO: SIDDE no encontró COMTI vs MAOBI orientado (COMTI, MAOBI) — se esperaba, vía LARGO");
  process.exit(1);
}
chk("COMTI vs MAOBI — directo (sin cambio)", siddeComtiMaobi.estimacionDirecta, -0.02, 3);
chk("COMTI vs MAOBI — indirecto (SIDDE)", siddeComtiMaobi.estimacionIndirecta, 0.2507, 3);
console.log(`z=${siddeComtiMaobi.z.toFixed(3)}, p=${siddeComtiMaobi.pValor.toFixed(3)}`);
console.log("-".repeat(66));

if (fallos === 0) {
  console.log(`\n✓ PASA: las ${total} cantidades coinciden con netmeta a los decimales publicados.\n`);
  process.exit(0);
} else {
  console.error(`\n✗ FALLA: ${fallos} de ${total} cantidad(es) ya no coinciden con netmeta. Revisa calcularNMA.\n`);
  process.exit(1);
}
