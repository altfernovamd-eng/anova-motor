// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Pruebas del meta-análisis en red (red.ts), reescrito el 7 de octubre de 2026 desde las fórmulas publicadas. La comparación numérica con
 * netmeta está en `npm run benchmark:red` (Stowe2010) y en `npm run validacion` (40 redes); aquí van las propiedades que no dependen de R:
 * una red en triángulo resuelta a mano, consistencia, identidades de la descomposición de Q, invariancias y los casos que se descartan.
 *   npm run test:red
 */
import { calcularNMA, calcularNodeSplittingSIDDE, calcularPScores } from "../src/red.ts";
import type { BrazoRed, EstudioRed } from "../src/tipos.ts";

let fallos = 0;
const ok = (c: boolean, m: string) => { console.log(`${c ? "✓" : "✗"} ${m}`); if (!c) fallos++; };
const cerca = (a: number, b: number, tol: number) => Math.abs(a - b) <= tol * Math.max(1, Math.abs(a), Math.abs(b));

const dm = (t: string, media: number, de: number, n: number): BrazoRed => ({ tratamiento: t, media: String(media), de: String(de), n: String(n) });
const est = (id: number, brazos: BrazoRed[], etiqueta = `E${id}`): EstudioRed => ({ id, etiqueta, medida: "DM", brazos });
// Varianza de un brazo = de²/n
const v = (de: number, n: number) => (de * de) / n;

// ---------- Triángulo A-B-C con un estudio por lado, resuelto a mano (efectos fijos) ----------
// Medias: A=0, B=1, C=3 en estudios distintos con un desvío conveniente → y_AB = 0−1, y_BC = 1−3, y_AC = 0−3 (tratamiento1 − tratamiento2)
const tri = [
  est(1, [dm("A", 0.0, 2, 50), dm("B", 1.0, 2, 50)]),
  est(2, [dm("B", 1.2, 2, 40), dm("C", 3.0, 2, 40)]),
  est(3, [dm("A", 0.1, 2, 60), dm("C", 2.5, 2, 60)]),
];
{
  const r = calcularNMA(tri, "DM", "fijo")!;
  const vAB = v(2, 50) * 2, vBC = v(2, 40) * 2, vAC = v(2, 60) * 2;
  const yAB = 0 - 1, yBC = 1.2 - 3, yAC = 0.1 - 2.5;
  // Estimado de red de A−B: combina el directo (peso 1/vAB) con el indirecto yAC − yBC (varianza vAC + vBC)
  const wD = 1 / vAB, wI = 1 / (vAC + vBC);
  const ab = (wD * yAB + wI * (yAC - yBC)) / (wD + wI);
  const seAB = Math.sqrt(1 / (wD + wI));
  const par = r.estimaciones.find((e) => e.tratamiento1 === "A" && e.tratamiento2 === "B")!;
  ok(cerca(par.diff, ab, 1e-10) && cerca(par.se, seAB, 1e-10), `triángulo: A−B de red = combinación directo + indirecto calculada a mano (${par.diff.toFixed(6)})`);
  ok(r.n === 3 && r.k === 3 && r.m === 3 && r.gl === 1 && r.tau2 === 0, "triángulo: 3 tratamientos, 3 estudios, 1 gl de inconsistencia, τ² = 0 en efectos fijos");
  // Q del triángulo: la inconsistencia (yAB + yBC − yAC) pesada por la varianza de la suma de los tres
  const incons = yAB + yBC - yAC;
  const qMano = (incons * incons) / (vAB + vBC + vAC);
  ok(cerca(r.q, qMano, 1e-10), `triángulo: Q = (inconsistencia)² / Σ varianzas = ${qMano.toFixed(6)}`);
  const ac = r.estimaciones.find((e) => e.tratamiento1 === "A" && e.tratamiento2 === "C")!;
  const bc = r.estimaciones.find((e) => e.tratamiento1 === "B" && e.tratamiento2 === "C")!;
  ok(cerca(par.diff + bc.diff, ac.diff, 1e-12), "consistencia: (A−B) + (B−C) = (A−C) en las estimaciones de red");
  const dq = r.decomposicionQ!;
  ok(cerca(dq.qDentro + dq.qEntre, r.q, 1e-12) && dq.glDentro === 0 && dq.glEntre === 1 && dq.disenos.length === 0, "Q = Q dentro + Q entre; sin diseños repetidos todo es inconsistencia entre diseños");
  const ps = calcularPScores(r);
  ok(cerca(ps.reduce((s, p) => s + p.pScore, 0), 1.5, 1e-12) && ps.every((p) => p.pScore >= 0 && p.pScore <= 1), "los P-scores suman n/2 = 1,5 y están en [0, 1]");
  ok(ps.find((p) => p.tratamiento === "A")!.pScore > ps.find((p) => p.tratamiento === "C")!.pScore, "menor es mejor: A (media más baja) queda por encima de C");
  ok(calcularPScores(r, false).find((p) => p.tratamiento === "C")!.pScore > calcularPScores(r, false).find((p) => p.tratamiento === "A")!.pScore, "mayor es mejor: se invierte el orden");
  // node-splitting: el par A−B tiene 1 estudio directo; P = Var(red)/Var(directa)
  const ns = r.nodeSplitting.find((e) => e.tratamiento1 === "A" && e.tratamiento2 === "B")!;
  ok(cerca(ns.proporcionDirecta, wD / (wD + wI), 1e-10), `node-splitting: proporción directa = peso directo / peso total (${ns.proporcionDirecta.toFixed(4)})`);
  ok(cerca(ns.diffDirecto, yAB, 1e-12) && cerca(ns.estimacionIndirecta, yAC - yBC, 1e-9), "node-splitting: el directo es el del estudio y el indirecto es A−C menos B−C");
  const sidde = calcularNodeSplittingSIDDE(tri, "DM", "fijo");
  ok(sidde.length === 0 || sidde.every((x) => Number.isFinite(x.z)), "SIDDE en un triángulo: sin estudios restantes no hay red que reajustar (no inventa)");
}

// ---------- Red más grande con un estudio de tres brazos ----------
const red4 = [
  est(1, [dm("A", 0.0, 1.8, 80), dm("B", 0.6, 1.8, 80)]),
  est(2, [dm("A", 0.1, 2.0, 90), dm("B", 0.8, 2.0, 90)]),
  est(3, [dm("B", 0.5, 1.9, 70), dm("C", 1.2, 1.9, 70)]),
  est(4, [dm("A", 0.0, 2.1, 100), dm("B", 0.7, 2.1, 100), dm("C", 1.4, 2.1, 100)]),
  est(5, [dm("C", 1.0, 2.0, 60), dm("D", 1.9, 2.0, 60)]),
  est(6, [dm("B", 0.4, 1.7, 75), dm("D", 1.6, 1.7, 75)]),
  est(7, [dm("A", 0.2, 2.2, 85), dm("D", 2.0, 2.2, 85)]),
];
{
  const f = calcularNMA(red4, "DM", "fijo")!, a = calcularNMA(red4, "DM", "aleatorio")!;
  ok(f.n === 4 && f.k === 7 && f.m === 6 + 3 && f.gl === (6 + 2) - 3, `4 tratamientos, 7 estudios, 9 contrastes (el de tres brazos aporta 3) y gl = Σ(kₛ − 1) − (n − 1) = ${f.gl}`);
  const dq = f.decomposicionQ!;
  ok(cerca(dq.qDentro + dq.qEntre, f.q, 1e-10) && dq.glDentro + dq.glEntre === f.gl, "Q total = Q dentro + Q entre y los gl suman");
  ok(dq.disenos.length === 1 && dq.disenos[0].tratamientos.join() === "A,B" && dq.disenos[0].nEstudios === 2 && dq.disenos[0].gl === 1, "único diseño con 2 estudios: A-B; su gl es 1");
  // para un diseño de dos brazos, la Q dentro es la Q de Cochran del par
  const y1 = 0.0 - 0.6, y2 = 0.1 - 0.8, w1 = 1 / (2 * v(1.8, 80)), w2 = 1 / (2 * v(2.0, 90));
  const ym = (w1 * y1 + w2 * y2) / (w1 + w2);
  ok(cerca(dq.disenos[0].q, w1 * (y1 - ym) ** 2 + w2 * (y2 - ym) ** 2, 1e-10), "Q del diseño A-B = Q de Cochran de sus dos estudios");
  // la evidencia directa de A−B junta los contrastes A−B de los estudios 1, 2 y 4 (este de tres brazos) con su varianza PROPIA (vₐ + v_b),
  // no con el peso reducido del estudio multibrazo
  {
    const ys = [0.0 - 0.6, 0.1 - 0.8, 0.0 - 0.7], vs = [2 * v(1.8, 80), 2 * v(2.0, 90), 2 * v(2.1, 100)];
    const ws = vs.map((x) => 1 / x), sw = ws.reduce((s, x) => s + x, 0);
    const directo = ws.reduce((s, w, i) => s + w * ys[i], 0) / sw;
    const nsAB = f.nodeSplitting.find((e) => e.tratamiento1 === "A" && e.tratamiento2 === "B")!;
    ok(cerca(nsAB.diffDirecto, directo, 1e-10) && cerca(nsAB.seDirecto, Math.sqrt(1 / sw), 1e-10) && nsAB.nEstudiosDirectos === 3, `node-splitting: el directo A−B pondera los 3 contrastes con su varianza propia (${directo.toFixed(5)})`);
  }
  ok(a.tau2 >= 0 && a.modelo === "aleatorio" && f.tau2 === a.tau2, "el τ² se estima igual en los dos modelos (el fijo no lo usa)");
  const seF = f.estimaciones.find((e) => e.tratamiento1 === "A" && e.tratamiento2 === "D")!.se;
  const seA = a.estimaciones.find((e) => e.tratamiento1 === "A" && e.tratamiento2 === "D")!.se;
  ok(a.tau2 > 0 ? seA > seF : cerca(seA, seF, 1e-12), "efectos aleatorios: el EE crece cuando τ² > 0 y es igual cuando τ² = 0");
  // el estudio de tres brazos aislado reproduce sus propios contrastes
  const solo = calcularNMA([red4[3], est(8, [dm("A", 0, 1, 10), dm("B", 1, 1, 10)]), est(9, [dm("B", 1, 1, 10), dm("C", 2, 1, 10)])], "DM", "fijo");
  ok(solo !== null && solo.n === 3, "una red con un estudio de tres brazos y otros dos conectados se ajusta");
  // invariancias
  const barajada = calcularNMA([...red4].reverse().map((e) => ({ ...e, brazos: [...e.brazos].reverse() })), "DM", "aleatorio")!;
  const dif = a.estimaciones.map((e, i) => Math.abs(e.diff - barajada.estimaciones[i].diff));
  ok(Math.max(...dif) < 1e-10 && cerca(a.tau2, barajada.tau2, 1e-10) && cerca(a.q, barajada.q, 1e-10), "el orden de estudios y de brazos no cambia ningún resultado");
  const renombrada = calcularNMA(red4.map((e) => ({ ...e, brazos: e.brazos.map((b) => ({ ...b, tratamiento: { A: "Z", B: "Y", C: "X", D: "W" }[b.tratamiento]! })) })), "DM", "aleatorio")!;
  const ab = a.estimaciones.find((e) => e.tratamiento1 === "A" && e.tratamiento2 === "B")!;
  const zy = renombrada.estimaciones.find((e) => e.tratamiento1 === "Y" && e.tratamiento2 === "Z")!;
  ok(cerca(ab.diff, -zy.diff, 1e-10) && cerca(ab.se, zy.se, 1e-10), "renombrar los tratamientos (cambia el orden alfabético) solo invierte el signo de los pares");
  const sidde = calcularNodeSplittingSIDDE(red4, "DM", "aleatorio");
  ok(sidde.every((x) => Number.isFinite(x.z) && x.pValor >= 0 && x.pValor <= 1 && Number.isNaN(x.proporcionDirecta)), "SIDDE: z y p válidos, sin proporción directa");
}

// ---------- Estudios con la misma etiqueta y estudios que se descartan ----------
{
  const dup = red4.map((e) => ({ ...e, etiqueta: "Mismo nombre" }));
  const a = calcularNMA(red4, "DM", "aleatorio")!, b = calcularNMA(dup, "DM", "aleatorio")!;
  ok(b.k === 7 && cerca(a.tau2, b.tau2, 1e-12) && cerca(a.q, b.q, 1e-12) && b.gl === a.gl, "estudios con la misma etiqueta no se confunden (antes se contaban como uno solo)");
  const conRepetido = [...red4, est(20, [dm("A", 0, 1, 10), dm("A", 1, 1, 10), dm("B", 1, 1, 10)], "Repetido"), est(21, [dm("A", 0, 1, 10), { tratamiento: "B", media: "x", de: "1", n: "10" }], "Sin datos")];
  const c = calcularNMA(conRepetido, "DM", "aleatorio")!;
  ok(c.k === 7 && c.descartados.length === 2 && c.descartados.some((d) => d.etiqueta === "Repetido" && /mismo tratamiento/.test(d.motivo)) && c.descartados.some((d) => d.etiqueta === "Sin datos" && /menos de dos brazos/.test(d.motivo)), "un estudio con el mismo tratamiento en dos brazos o con un solo brazo válido se descarta y se declara");
  ok(cerca(c.tau2, a.tau2, 1e-12), "…y no cambia el resultado de los demás");
}

// ---------- OR con ceros: la regla de `meta::pairwise` + `netmeta` (valores de R, red base de 5 estudios de 2 brazos) ----------
{
  const orb = (t: string, e: number, n: number): BrazoRed => ({ tratamiento: t, eventos: String(e), total: String(n) });
  const base: EstudioRed[] = [
    { id: 1, etiqueta: "S1", medida: "OR", brazos: [orb("A", 10, 80), orb("B", 14, 80)] },
    { id: 2, etiqueta: "S2", medida: "OR", brazos: [orb("B", 12, 70), orb("C", 9, 70)] },
    { id: 3, etiqueta: "S3", medida: "OR", brazos: [orb("C", 20, 90), orb("D", 11, 90)] },
    { id: 4, etiqueta: "S4", medida: "OR", brazos: [orb("A", 8, 60), orb("D", 15, 60)] },
    { id: 5, etiqueta: "S5", medida: "OR", brazos: [orb("A", 13, 100), orb("C", 9, 100)] },
  ];
  const con = (brazos: BrazoRed[]) => calcularNMA([...base, { id: 9, etiqueta: "T", medida: "OR", brazos }], "OR", "fijo")!;
  const ab = (r: ReturnType<typeof con>) => r.estimaciones.find((e) => e.tratamiento1 === "A" && e.tratamiento2 === "B")!;
  // Cada caso: (brazos del estudio T, k, m y ln OR de A:B con efecto común en netmeta 3.7.0)
  const casos: [string, BrazoRed[], number, number, number][] = [
    ["3 brazos, un brazo con 0 eventos: se corrigen TODOS sus contrastes (también el que no tiene ceros)", [orb("A", 0, 60), orb("B", 6, 60), orb("C", 9, 60)], 6, 8, -0.4301594834],
    ["3 brazos, un contraste de doble cero (A y B con 0): se conserva con su corrección (allstudies = TRUE)", [orb("A", 0, 60), orb("B", 0, 60), orb("C", 7, 60)], 6, 8, -0.41135107],
    ["3 brazos, doble cero entre B y C", [orb("A", 5, 60), orb("B", 0, 60), orb("C", 0, 60)], 6, 8, -0.26200967],
    ["2 brazos, todos con evento en los dos: no informa y se descarta", [orb("A", 30, 30), orb("B", 30, 30)], 5, 5, -0.46114961],
    ["2 brazos, 0 en un brazo y todos con evento en el otro: se conserva corregido", [orb("A", 0, 30), orb("B", 30, 30)], 6, 6, -0.70108396],
    ["2 brazos, un brazo con todos los eventos: se corrige", [orb("A", 30, 30), orb("B", 12, 30)], 6, 6, -0.18078681],
  ];
  for (const [nombre, brazos, k, m, esperado] of casos) {
    const r = con(brazos);
    ok(r.k === k && r.m === m && cerca(ab(r).diff, esperado, 1e-7), `netmeta: ${nombre} (k = ${r.k}, m = ${r.m}, ln OR A−B = ${ab(r).diff.toFixed(8)})`);
  }
  const doble = con([orb("A", 0, 60), orb("B", 0, 60)]);
  ok(doble.descartados.length === 1 && /dos brazos/.test(doble.descartados[0].motivo), "el doble cero de dos brazos se declara en `descartados`");
  const todosCero = con([orb("A", 0, 60), orb("B", 0, 60), orb("C", 0, 60)]);
  ok(todosCero.descartados.length === 1 && /todos los brazos/.test(todosCero.descartados[0].motivo) && todosCero.k === 5, "un estudio de tres brazos sin eventos en ninguno se descarta y se declara");
}

// ---------- Redes que no se calculan ----------
ok(calcularNMA(tri.slice(0, 1), "DM") === null, "un solo estudio → null");
ok(calcularNMA([est(1, [dm("A", 0, 1, 10), dm("B", 1, 1, 10)]), est(2, [dm("A", 0, 1, 10), dm("B", 1, 1, 10)])], "DM") === null, "menos de 3 tratamientos → null");
ok(calcularNMA([est(1, [dm("A", 0, 1, 10), dm("B", 1, 1, 10)]), est(2, [dm("C", 0, 1, 10), dm("D", 1, 1, 10)]), est(3, [dm("A", 0, 1, 10), dm("B", 1.1, 1, 10)])], "DM") === null, "red desconectada (A-B y C-D sin puente) → null");

// ---------- OR: ceros ----------
{
  const orb = (t: string, e: number, n: number): BrazoRed => ({ tratamiento: t, eventos: String(e), total: String(n) });
  const r = calcularNMA([
    { id: 1, etiqueta: "O1", medida: "OR", brazos: [orb("A", 0, 100), orb("B", 5, 100)] },
    { id: 2, etiqueta: "O2", medida: "OR", brazos: [orb("B", 4, 80), orb("C", 10, 80)] },
    { id: 3, etiqueta: "O3", medida: "OR", brazos: [orb("A", 2, 90), orb("C", 9, 90)] },
    { id: 4, etiqueta: "O4", medida: "OR", brazos: [orb("A", 3, 70), orb("B", 6, 70), orb("C", 11, 70)] },
  ], "OR", "fijo")!;
  ok(r.escala === "log" && r.estimaciones.every((e) => e.estimacion > 0 && e.icInferior < e.estimacion && e.estimacion < e.icSuperior), "OR: escala log, estimaciones positivas e IC que contienen la estimación (el brazo con 0 eventos lleva +0,5)");
}

if (fallos) { console.error(`\n${fallos} prueba(s) fallaron.`); process.exit(1); }
console.log("\nTodas las pruebas pasaron.");
