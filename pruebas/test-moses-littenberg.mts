// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Pruebas de la regresión de Moses-Littenberg y de las piezas por estudio de la SROC univariada (moses-littenberg.ts), reescrita el 7 de
 * octubre de 2026 desde los artículos. La comparación con R (`mada`) está en `npm run validacion` y `npm run benchmark:sroc`.
 *   npm run test:moses
 */
import { curvaMosesLittenberg, intervaloWilson, puntoDiagnostico, regresionMosesLittenberg } from "../src/moses-littenberg.ts";

let fallos = 0;
const ok = (c: boolean, m: string) => { console.log(`${c ? "✓" : "✗"} ${m}`); if (!c) fallos++; };
const cerca = (a: number, b: number, tol: number) => Math.abs(a - b) <= tol * Math.max(1, Math.abs(a), Math.abs(b));
const logit = (p: number) => Math.log(p / (1 - p));

// ---------- Wilson (Newcombe 1998: 81/263 → [0,2553; 0,3662]) ----------
{
  const w = intervaloWilson(81, 263);
  ok(cerca(w.icInferior, 0.2553, 2e-4) && cerca(w.icSuperior, 0.3662, 2e-4) && cerca(w.estimacion, 81 / 263, 1e-12), `Wilson 81/263 = [${w.icInferior.toFixed(4)}; ${w.icSuperior.toFixed(4)}] (Newcombe 1998)`);
  const cero = intervaloWilson(0, 20), todo = intervaloWilson(20, 20);
  ok(cero.icInferior === 0 && cero.icSuperior > 0 && todo.icSuperior === 1 && todo.icInferior < 1, "Wilson con 0 o n eventos: el límite del lado de la frontera queda en 0 o 1 y el otro no colapsa");
}

// ---------- Punto por estudio ----------
{
  const p = puntoDiagnostico(45, 5, 5, 45)!; // VP=45 FP=5 FN=5 VN=45: sens 0,9; FPR 0,1
  ok(cerca(p.dLn, Math.log((45 * 45) / (5 * 5)), 1e-12) && cerca(p.dLn, logit(0.9) - logit(0.1), 1e-12), "D = ln(DOR) = logit(sens) − logit(FPR)");
  ok(cerca(p.s, logit(0.9) + logit(0.1), 1e-12) && cerca(p.s, 0, 1e-12), "S = logit(sens) + logit(FPR) (0 cuando sens = esp)");
  ok(cerca(p.seD, Math.sqrt(1 / 45 + 1 / 5 + 1 / 5 + 1 / 45), 1e-12) && !p.correccion, "EE de D = √(1/a + 1/b + 1/c + 1/d); sin ceros no se corrige");
  const z = puntoDiagnostico(0, 5, 5, 45)!;
  ok(z.correccion && z.a === 0.5 && z.b === 5.5 && z.c === 5.5 && z.d === 45.5, "una celda en cero: +0,5 en las cuatro celdas");
  ok(puntoDiagnostico(0, 0, 0, 10) === null && puntoDiagnostico(5, 0, 3, 0) === null && puntoDiagnostico(-1, 2, 3, 4) === null && puntoDiagnostico(NaN, 2, 3, 4) === null, "sin enfermos, sin sanos, negativos o NaN → null");
}

// ---------- Regresión: contra la forma centrada, calculada aparte ----------
const tablas: [number, number, number, number][] = [[45, 5, 5, 45], [40, 8, 10, 42], [90, 10, 12, 88], [30, 4, 6, 60], [70, 9, 8, 80], [55, 14, 5, 70], [25, 2, 11, 52], [62, 18, 4, 40]];
const pts = tablas.map(([vp, fp, fn, vn]) => puntoDiagnostico(vp, fp, fn, vn)!);
{
  const r = regresionMosesLittenberg(pts)!;
  const w = pts.map((p) => 1 / (p.seD * p.seD)), sw = w.reduce((s, x) => s + x, 0);
  const sBar = pts.reduce((s, p, i) => s + w[i] * p.s, 0) / sw, dBar = pts.reduce((s, p, i) => s + w[i] * p.dLn, 0) / sw;
  const sxx = pts.reduce((s, p, i) => s + w[i] * (p.s - sBar) ** 2, 0);
  const b = pts.reduce((s, p, i) => s + w[i] * (p.s - sBar) * (p.dLn - dBar), 0) / sxx;
  const a = dBar - b * sBar;
  ok(cerca(r.b, b, 1e-10) && cerca(r.a, a, 1e-10), `WLS: α = ${r.a.toFixed(5)}, β = ${r.b.toFixed(5)} coinciden con la forma centrada`);
  ok(cerca(r.seB, Math.sqrt(1 / sxx), 1e-10) && cerca(r.seA, Math.sqrt(1 / sw + (sBar * sBar) / sxx), 1e-10) && cerca(r.zB, b / Math.sqrt(1 / sxx), 1e-10), "EE de α y β y z de Wald con pesos conocidos");
  // los residuos ponderados son ortogonales a 1 y a S (condiciones normales)
  const res = pts.map((p) => p.dLn - (r.a + r.b * p.s));
  ok(cerca(res.reduce((s, e, i) => s + w[i] * e, 0), 0, 1e-9) && cerca(res.reduce((s, e, i) => s + w[i] * e * pts[i].s, 0), 0, 1e-9), "condiciones normales: Σ w·e = 0 y Σ w·e·S = 0");
  // invariancia al orden
  const inv = regresionMosesLittenberg([...pts].reverse())!;
  ok(cerca(inv.a, r.a, 1e-12) && cerca(inv.b, r.b, 1e-12), "el orden de los estudios no importa");
  ok(regresionMosesLittenberg(pts.slice(0, 2)) === null, "menos de 3 puntos → null");
  const iguales = [puntoDiagnostico(45, 5, 5, 45)!, puntoDiagnostico(90, 10, 10, 90)!, puntoDiagnostico(18, 2, 2, 18)!];
  ok(regresionMosesLittenberg(iguales) === null, "todos con el mismo S → null (la pendiente no es identificable)");
}

// ---------- Curva ----------
{
  const c = curvaMosesLittenberg({ a: 2, b: 0 })!;
  ok(c.length === 61 && cerca(c[0].fpr, 0.001, 1e-12) && cerca(c[60].fpr, 0.999, 1e-12), "61 puntos entre FPR 0,001 y 0,999");
  // con β = 0, logit(TPR) = α + logit(FPR): DOR constante e^α
  ok(c.every((q) => cerca(logit(q.tpr) - logit(q.fpr), 2, 1e-9)), "β = 0: la diferencia de logits es constante (DOR constante)");
  const c2 = curvaMosesLittenberg({ a: 1.5, b: 0.3 })!;
  // D = L_T − L_F y S = L_T + L_F cumplen D = α + βS en cada punto
  ok(c2.every((q) => { const D = logit(q.tpr) - logit(q.fpr), S = logit(q.tpr) + logit(q.fpr); return cerca(D, 1.5 + 0.3 * S, 1e-9); }), "cada punto de la curva satisface D = α + β·S");
  ok(curvaMosesLittenberg({ a: 1, b: 1 }) === null && curvaMosesLittenberg({ a: 1, b: -1.2 }) === null, "|β| ≥ 1 → no hay curva");
  ok(c2.every((q, i) => i === 0 || (q.tpr >= c2[i - 1].tpr && q.fpr > c2[i - 1].fpr)) && c2.every((q) => q.tpr > 0 && q.tpr < 1), "la curva crece con el FPR y queda dentro de (0, 1)");
}

if (fallos) { console.error(`\n${fallos} prueba(s) fallaron.`); process.exit(1); }
console.log("\nTodas las pruebas pasaron.");
