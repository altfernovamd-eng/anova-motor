// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Funciones especiales de probabilidad del motor (ln Γ, gamma y beta incompletas) contra R.
 *
 *   npm run test:funciones-especiales
 *
 * Referencia: validacion/referencias/referencia-funciones.json, generada por
 * validacion/generar-referencia-funciones.R con funciones base de R (pchisq, pt, pnorm, qt, lgamma).
 * Cubre colas desde p ≈ 1 hasta p ≈ 1e-290 y los puntos donde cambia el método de cálculo.
 * Además comprueba las propiedades matemáticas que no dependen de R (simetrías, identidades, casos límite).
 */
import { readFileSync } from "node:fs";
import {
  betaIncompleta,
  chiCuadradoCola,
  gammaIncompletaP,
  gammaIncompletaQ,
  lnGamma,
  normalCuantil,
  normalCuantilDosColas,
  normalDosColas,
  tStudentCuantilDosColas,
  tStudentDosColas,
} from "../src/funciones-especiales.ts";
import { chiSquarePValue, colaNormalDosLados, studentTQuantile975 } from "../src/statistics.ts";

const ref = JSON.parse(
  readFileSync(new URL("../validacion/referencias/referencia-funciones.json", import.meta.url), "utf8"),
);

let fallos = 0;
let comprobaciones = 0;
const ok = (cond: boolean, msg: string) => {
  comprobaciones++;
  if (!cond) {
    fallos++;
    console.log("  ✗", msg);
  }
};
const rel = (a: number, b: number) => (b === 0 ? Math.abs(a) : Math.abs(a - b) / Math.abs(b));

interface Fila {
  [k: string]: number;
}
function contraR(nombre: string, filas: Fila[], f: (r: Fila) => number, campo: string, tol: number) {
  let peor = 0;
  let donde = "";
  for (const r of filas) {
    const e = rel(f(r), r[campo]);
    if (!(e <= tol)) {
      fallos++;
      console.log(`  ✗ ${nombre}: error relativo ${e.toExponential(2)} en ${JSON.stringify(r)}`);
    }
    if (e > peor) {
      peor = e;
      donde = JSON.stringify(r);
    }
  }
  comprobaciones += filas.length;
  console.log(`${nombre.padEnd(34)} ${String(filas.length).padStart(4)} puntos · peor error relativo ${peor.toExponential(2)} (tolerancia ${tol.toExponential(0)})`);
  return donde;
}

console.log("— Contra R (funciones base) —");
contraR("χ² cola superior", ref.chi2, (r) => chiCuadradoCola(r.q, r.gl), "p", 1e-10);
contraR("t de Student, dos colas", ref.t, (r) => tStudentDosColas(r.t, r.gl), "p", 1e-10);
contraR("normal, dos colas", ref.normal, (r) => normalDosColas(r.z), "p", 1e-10);
contraR("ln Γ", ref.lgamma, (r) => lnGamma(r.x), "v", 1e-12);
contraR("cuantil bilateral de la t", ref.qt2, (r) => tStudentCuantilDosColas(r.alfa, r.gl), "q", 1e-10);
contraR("cuantil bilateral de la normal", ref.qz2, (r) => normalCuantilDosColas(r.alfa), "q", 1e-10);
contraR("cuantil de la normal Φ⁻¹(p)", ref.qn, (r) => normalCuantil(r.p), "q", 1e-12);

console.log("\n— Lo que usa el motor (statistics.ts) —");
contraR("chiSquarePValue", ref.chi2, (r) => chiSquarePValue(r.q, r.gl), "p", 1e-10);
contraR("colaNormalDosLados", ref.normal, (r) => colaNormalDosLados(r.z), "p", 1e-10);
contraR("studentTQuantile975", ref.t975, (r) => studentTQuantile975(r.gl), "q", 1e-12);

console.log("\n— Propiedades matemáticas (no dependen de R) —");
// P + Q = 1
for (const [a, x] of [[0.5, 0.3], [2, 1], [5, 9], [30, 25], [100, 130], [1000, 900]] as const) {
  ok(Math.abs(gammaIncompletaP(a, x) + gammaIncompletaQ(a, x) - 1) < 1e-12, `P + Q = 1 en a=${a}, x=${x}`);
}
// Para a = 1, Q(1, x) = e^(−x)
for (const x of [0.1, 1, 5, 20, 100]) ok(rel(gammaIncompletaQ(1, x), Math.exp(-x)) < 1e-13, `Q(1, ${x}) = e^−${x}`);
// Simetría de la beta: I_x(a, b) + I_(1−x)(b, a) = 1
for (const [x, a, b] of [[0.2, 2, 3], [0.7, 0.5, 4], [0.5, 10, 10], [0.01, 50, 0.5], [0.99, 3, 30]] as const) {
  ok(Math.abs(betaIncompleta(x, a, b) + betaIncompleta(1 - x, b, a) - 1) < 1e-12, `simetría de la beta en x=${x}, a=${a}, b=${b}`);
}
// I_x(1, 1) = x (uniforme) e I_x(a, 1) = x^a
for (const x of [0.1, 0.5, 0.9]) {
  ok(Math.abs(betaIncompleta(x, 1, 1) - x) < 1e-14, `I_${x}(1,1) = ${x}`);
  ok(rel(betaIncompleta(x, 3, 1), x ** 3) < 1e-12, `I_${x}(3,1) = ${x}³`);
}
// Valores exactos conocidos: χ²(2) tiene cola e^(−q/2); t con 1 gl tiene P(|T|>t) = 1 − (2/π)·atan(t)
for (const q of [0.5, 3, 10, 40]) ok(rel(chiCuadradoCola(q, 2), Math.exp(-q / 2)) < 1e-13, `χ² con 2 gl en q=${q}`);
for (const t of [0.3, 1, 4, 20]) ok(Math.abs(tStudentDosColas(t, 1) - (1 - (2 / Math.PI) * Math.atan(t))) < 1e-13, `t con 1 gl en t=${t}`);
// Γ(n) = (n−1)!, Γ(½) = √π
ok(Math.abs(lnGamma(10) - Math.log(362880)) < 1e-13, "ln Γ(10) = ln 9!");
ok(Math.abs(lnGamma(0.5) - 0.5 * Math.log(Math.PI)) < 1e-14, "ln Γ(½) = ½ ln π");
// Casos límite: sin números no válidos
ok(Number.isNaN(chiCuadradoCola(3, 0)), "χ² con 0 gl → NaN");
ok(chiCuadradoCola(0, 5) === 1 && chiCuadradoCola(-1, 5) === 1, "q ≤ 0 → p = 1");
ok(tStudentDosColas(0, 10) === 1, "t = 0 → p = 1");
ok(Number.isNaN(tStudentDosColas(2, -1)), "t con gl negativos → NaN");
ok(normalDosColas(0) === 1, "z = 0 → p = 1");
ok(normalDosColas(40) === 0 || normalDosColas(40) < 1e-300, "z = 40 → ≈ 0 sin romper");
ok(Number.isFinite(chiCuadradoCola(1e6, 5)) && chiCuadradoCola(1e6, 5) >= 0, "q enorme no rompe");
// Monotonía: el p baja al subir q
let anterior = 1;
let monotono = true;
for (let q = 0; q < 60; q += 0.37) {
  const p = chiCuadradoCola(q, 7);
  if (p > anterior + 1e-15) monotono = false;
  anterior = p;
}
ok(monotono, "el p de χ² baja al subir q");

console.log(`\n${comprobaciones - fallos}/${comprobaciones} comprobaciones correctas`);
if (fallos > 0) {
  console.log(`✗ ${fallos} fallos`);
  process.exit(1);
}
console.log("✓ Funciones especiales de acuerdo con R en todo el rango probado");
