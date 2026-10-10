// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * κ de Cohen y κ de Fleiss contra R (`irr::kappa2` y `irr::kappam.fleiss`) en conjuntos sintéticos con semilla fija.
 *
 *   npm run test:kappa-r
 *
 * Referencia: validacion/referencias/referencia-kappa.json, generada por validacion/generar-referencia-kappa.R
 * (el paquete `irr` se usa como caja negra: solo se guardan sus resultados). Además comprueba valores calculables a mano y
 * propiedades que no dependen de R.
 */
import { readFileSync } from "node:fs";
import { cohenKappa, fleissKappa } from "../src/kappa.ts";

const ref = JSON.parse(readFileSync(new URL("../validacion/referencias/referencia-kappa.json", import.meta.url), "utf8"));

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

const peor = { kappa: 0, z: 0, p: 0 };
const n = { cohen: 0, fleiss: 0 };
for (const c of ref.casos as { tipo: "cohen" | "fleiss"; datos: number[][]; kappa: number; z: number; p: number; pExacto: number }[]) {
  const filas = c.datos.map((f) => f.map((v) => `c${v}`));
  const r =
    c.tipo === "cohen" ? cohenKappa(filas.map((f) => [f[0], f[1]] as [string, string])) : fleissKappa(filas);
  n[c.tipo]++;
  if (!r) {
    ok(false, `${c.tipo}: devolvió null en un caso con referencia (${JSON.stringify(c.datos).slice(0, 60)}…)`);
    continue;
  }
  const dk = Math.abs(r.kappa - c.kappa);
  // z: tolerancia absoluta pequeña además de la relativa (con κ = 0, irr devuelve ruido numérico del orden de 1e-15).
  const dz = Math.abs(r.z - c.z) / Math.max(1, Math.abs(c.z));
  // p: contra la cola exacta de R, 2·pnorm(−|z|); contra el p que imprime irr solo donde éste no pierde cifras.
  const dp = rel(r.pValor, c.pExacto);
  if (c.p >= 1e-6) ok(rel(r.pValor, c.p) < 1e-9, `${c.tipo}: p difiere del que imprime irr (aNova ${r.pValor}, irr ${c.p})`);
  peor.kappa = Math.max(peor.kappa, dk);
  peor.z = Math.max(peor.z, dz);
  peor.p = Math.max(peor.p, dp);
  ok(dk < 1e-12, `${c.tipo}: κ difiere de irr en ${dk.toExponential(2)}`);
  ok(dz < 1e-9, `${c.tipo}: z difiere de irr (error ${dz.toExponential(2)}; aNova ${r.z}, irr ${c.z})`);
  ok(dp < 1e-9, `${c.tipo}: p difiere de la cola exacta de R (error relativo ${dp.toExponential(2)}; aNova ${r.pValor}, exacto ${c.pExacto})`);
}
console.log(`irr ${ref.irr}: ${n.cohen} conjuntos de Cohen y ${n.fleiss} de Fleiss`);
console.log(`peor diferencia: κ ${peor.kappa.toExponential(2)} (absoluta) · z ${peor.z.toExponential(2)} · p ${peor.p.toExponential(2)} (relativas)`);

console.log("\n— Valores calculables a mano —");
// 2×2: 20 de acuerdo en (a,a), 15 en (b,b), 5 y 10 en discordancias → N=50, po=0,70, filas (25,25)… ver cuentas abajo
{
  const pares: [string, string][] = [
    ...Array(20).fill(["a", "a"]),
    ...Array(5).fill(["a", "b"]),
    ...Array(10).fill(["b", "a"]),
    ...Array(15).fill(["b", "b"]),
  ];
  const r = cohenKappa(pares)!;
  // Marginales: A: a=25, b=25 · B: a=30, b=20 → pe = (25·30 + 25·20)/2500 = 0,5 · po = 35/50 = 0,7 → κ = 0,4
  ok(Math.abs(r.acuerdoObservado - 0.7) < 1e-15, "Cohen 2×2: po = 0,7");
  ok(Math.abs(r.acuerdoEsperado - 0.5) < 1e-15, "Cohen 2×2: pe = 0,5");
  ok(Math.abs(r.kappa - 0.4) < 1e-15, "Cohen 2×2: κ = 0,4");
  // Varianza bajo H0 con la fórmula desarrollada: (pe + pe² − Σ pᵢ. p.ᵢ (pᵢ. + p.ᵢ)) / (N(1−pe)²)
  //   Σ = 0,5·0,6·(0,5+0,6) + 0,5·0,4·(0,5+0,4) = 0,33 + 0,18 = 0,51 → (0,5 + 0,25 − 0,51)/(50·0,25) = 0,24/12,5 = 0,0192
  ok(Math.abs(r.se - Math.sqrt(0.0192)) < 1e-14, "Cohen 2×2: SE bajo H0 = √0,0192");
}
// Fleiss: 4 sujetos, 3 evaluadores, 2 categorías
{
  const r = fleissKappa([["x", "x", "x"], ["x", "x", "y"], ["x", "y", "y"], ["y", "y", "y"]])!;
  // P1=1, P2=1/3, P3=1/3, P4=1 → po = 2/3 · p_x = 6/12 = 0,5 → pe = 0,5 → κ = (2/3 − 1/2)/(1/2) = 1/3
  ok(Math.abs(r.acuerdoObservado - 2 / 3) < 1e-15, "Fleiss: po = 2/3");
  ok(Math.abs(r.acuerdoEsperado - 0.5) < 1e-15, "Fleiss: pe = 0,5");
  ok(Math.abs(r.kappa - 1 / 3) < 1e-15, "Fleiss: κ = 1/3");
  // Var0 = 2(pe + pe² − 2Σp³) / (N n (n−1) (1−pe)²) = 2(0,5 + 0,25 − 0,5)/(4·3·2·0,25) = 0,5/6
  ok(Math.abs(r.se - Math.sqrt(0.5 / 6)) < 1e-14, "Fleiss: SE bajo H0 = √(0,5/6)");
}

console.log("\n— Propiedades —");
{
  const perfecto = cohenKappa([["a", "a"], ["b", "b"], ["c", "c"], ["a", "a"]])!;
  ok(Math.abs(perfecto.kappa - 1) < 1e-15, "acuerdo perfecto → κ = 1");
  ok(cohenKappa([["a", "a"], ["a", "a"]]) === null, "una sola categoría (pe = 1) → null, sin dividir entre 0");
  ok(cohenKappa([]) === null, "sin sujetos → null");
  ok(fleissKappa([["a", "b"], ["a"]]) === null, "evaluadores desiguales → null");
  ok(fleissKappa([["a"]]) === null, "un solo evaluador → null");
  const contrario = cohenKappa([["a", "b"], ["b", "a"], ["a", "b"], ["b", "a"]])!;
  ok(contrario.kappa < 0 && contrario.interpretacion.startsWith("Pobre"), "acuerdo peor que el azar → κ < 0, «Pobre»");
  const simetrico = [["a", "a"], ["a", "b"], ["b", "a"], ["b", "b"], ["a", "a"], ["b", "b"]] as [string, string][];
  const dirA = cohenKappa(simetrico)!;
  const dirB = cohenKappa(simetrico.map(([x, y]) => [y, x] as [string, string]))!;
  ok(Math.abs(dirA.kappa - dirB.kappa) < 1e-15 && Math.abs(dirA.se - dirB.se) < 1e-15, "intercambiar evaluadores no cambia κ ni el SE");
  ok(Math.abs(cohenKappa(simetrico.map(([x, y]) => [x.toUpperCase(), y.toUpperCase()] as [string, string]))!.kappa - dirA.kappa) < 1e-15, "renombrar las categorías no cambia κ");
}

console.log(`\n${comprobaciones - fallos}/${comprobaciones} comprobaciones correctas`);
if (fallos > 0) {
  console.log(`✗ ${fallos} fallos`);
  process.exit(1);
}
console.log("✓ κ de acuerdo con irr en todos los conjuntos y con los cálculos a mano");
