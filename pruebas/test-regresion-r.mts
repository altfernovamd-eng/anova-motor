// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
// Regresión de los errores que la validación diferencial contra R (oct 2026) encontró en el motor. Cada caso es un
// conjunto sintético de `validacion/referencias/referencia.json` con su valor de R (`metafor`/`meta`) copiado
// aquí, para que el test sea rápido y no dependa de leer el JSON. Ejecutar: npm run test:regresion-r
import {
  efectoDesde2x2,
  calcularMetaAnalisis,
  colaNormalDosLados,
  formatoParaGuardar,
  metaAnalisisMH,
  testEgger,
  testPeters,
  trimAndFill,
} from "../src/statistics.ts";
import { calcularNMA, calcularNodeSplittingSIDDE } from "../src/red.ts";
import type { Estudio } from "../src/tipos.ts";

let fallos = 0;
const ok = (c: boolean, m: string) => { console.log(`${c ? "✓" : "✗"} ${m}`); if (!c) fallos++; };
const cerca = (a: number, b: number, tolRel: number, tolAbs = 1e-12) => Math.abs(a - b) <= Math.max(tolAbs, tolRel * Math.max(Math.abs(a), Math.abs(b)));

const Z95 = 1.959963985;
let n = 0;
const est = (medida: Estudio["medida"], valor: number, lo: number, hi: number, extra: Partial<Estudio> = {}): Estudio =>
  ({ id: ++n, autor: `E${n}`, año: "", diseño: "ECA", n: "", medida, valor: String(valor), icInferior: String(lo), icSuperior: String(hi), ...extra }) as Estudio;

// ---- Conjunto genérico iv054 (20 estudios, escala lineal) ----
const yi = [0.257259079292431, 0.192182326596628, -0.185285637356688, 0.195110178097324, -0.419121636468224, -0.161740651138161, -0.837046585950171, -0.0908515864328519, -0.64085784247493, -0.519788024638282, 0.676708126653011, -0.133965349896889, -0.57135880709784, -0.004405827342237, -1.14211899023374, -1.32843327223507, -0.0468822847418429, -0.682223300133511, -0.15512283397462, -0.748067918176515];
const sei = [0.213530755450483, 0.505586704879533, 0.396812941762619, 0.348059448867571, 0.229450540314428, 0.176190611533821, 0.673025806003716, 0.368049828510266, 0.102674784488045, 0.334527855436318, 0.063931182725355, 0.411534186964855, 0.30263108148938, 0.583566818735562, 0.575267992459703, 0.497907673835289, 0.0855696592247114, 0.456276253541, 0.527474167325999, 0.233670008019544];
const estudiosIv = yi.map((y, i) => est("DM", y, y - Z95 * sei[i], y + Z95 * sei[i]));

{
  // Error 1: el p de Egger salía DUPLICADO (2 × una cola que ya era bilateral).
  const eg = testEgger(calcularMetaAnalisis(estudiosIv, "DM", "DL")!);
  ok(!!eg, "Egger: devuelve resultado con 20 estudios");
  ok(!!eg && cerca(eg.pValor, 0.0310847814453766, 1e-6), `Egger: p = R (metafor::regtest) (0.0310847814453766)`);
  ok(!!eg && cerca(eg.intercepto, -2.25561359941439, 1e-6) && cerca(eg.t, -2.3387487750756, 1e-6), "Egger: intercepto y t coinciden con R");
}
{
  // Error 4: trim-and-fill con efectos aleatorios usaba el centro de EFECTOS FIJOS (aquí 10 estudios imputados en vez de 2).
  const fe = calcularMetaAnalisis(estudiosIv, "DM", "DL")!;
  const tf = trimAndFill(fe, "aleatorio", "DL");
  ok(!!tf && tf.k0 === 2, "Trim-and-fill DL: k0 = R (2, no los 10 del centro de efectos fijos)");
  ok(!!tf && cerca(tf.estimacionAjustada, -0.197427632893434, 1e-6), "Trim-and-fill DL: efecto ajustado = R (metafor::trimfill)");
}

// ---- Conjunto 2×2 bin067 (25 estudios; varios con celdas en cero y dobles-ceros) ----
const ai = [1, 0, 1, 0, 2, 0, 2, 0, 2, 1, 1, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0, 1, 2, 1, 1], n1 = [96, 110, 325, 188, 220, 59, 237, 91, 362, 173, 198, 174, 49, 375, 151, 259, 383, 164, 118, 49, 69, 61, 311, 378, 241], ci = [3, 1, 5, 0, 0, 0, 3, 0, 1, 3, 2, 1, 0, 2, 1, 3, 2, 1, 3, 0, 0, 0, 3, 6, 2], n2 = [101, 118, 323, 168, 220, 38, 216, 79, 318, 188, 220, 186, 32, 365, 142, 275, 408, 170, 103, 63, 81, 71, 314, 381, 232];
const crudos = (i: number): Estudio["crudos"] => ({ tipo: "2x2", eventosTrat: String(ai[i]), totalTrat: String(n1[i]), eventosCtrl: String(ci[i]), totalCtrl: String(n2[i]) });
const orEst = ai.map((_, i) => {
  const e = efectoDesde2x2(ai[i], n1[i], ci[i], n2[i], "OR")!;
  return est("OR", e.valor, e.icInferior, e.icSuperior, { crudos: crudos(i) });
});
const rarEst = ai.map((_, i) => {
  const e = efectoDesde2x2(ai[i], n1[i], ci[i], n2[i], "RAR")!;
  return est("RAR", e.valor, e.icInferior, e.icSuperior, { crudos: crudos(i) });
});

{
  // Error 1 (Peters) + convención: el peso usa los conteos originales; el +0,5 solo corrige el ln(OR).
  const pt = testPeters(orEst, "OR");
  ok(!!pt, "Peters: devuelve resultado");
  ok(!!pt && cerca(pt.pValor, 0.614048268118289, 1e-4) && cerca(pt.t, 0.513737412066413, 1e-5), "Peters: t y p = R (meta::metabias), con celdas en cero");
}
{
  // Error 3: la diferencia de riesgos de Mantel-Haenszel excluía las tablas doble-cero; SE y Q según RevMan 5.
  const mh = metaAnalisisMH(rarEst, "RAR");
  ok(!!mh && mh.k === 25 && mh.excluidos === 0, "MH diferencia de riesgos: las 25 tablas entran (doble-cero incluidas)");
  ok(!!mh && cerca(mh.estimacion, -0.00436779990594252, 1e-9), "MH diferencia de riesgos: estimación = R (rma.mh)");
  ok(!!mh && cerca(mh.se, 0.00163565940811852, 1e-9), "MH diferencia de riesgos: SE = fórmula de RevMan 5 (Greenland-Robins)");
  ok(!!mh && cerca(mh.heterogeneidad.q, 13.3591068679277, 1e-8), "MH diferencia de riesgos: Q = RevMan 5 (doble-cero con +0,5)");
}
{
  // Error 2: el efecto se guardaba con 4 decimales y el combinado se apartaba de R hasta ~1,6 %.
  const guardados = ai.map((_, i) => {
    const e = efectoDesde2x2(ai[i], n1[i], ci[i], n2[i], "OR")!;
    const g = (x: number) => Number(formatoParaGuardar(x));
    return est("OR", g(e.valor), g(e.icInferior), g(e.icSuperior), { crudos: crudos(i) }); // la app guarda también los crudos
  });
  const a = calcularMetaAnalisis(guardados, "OR", "DL")!;
  ok(cerca(Math.log(a.aleatorio.estimacion), Math.log(calcularMetaAnalisis(orEst, "OR", "DL")!.aleatorio.estimacion), 1e-7), "Guardado a 9 cifras: el combinado DL no se mueve (antes hasta 1,6 %)");
  ok(formatoParaGuardar(0.123456789012) === "0.123456789" && formatoParaGuardar(123456.789012) === "123456.789", "formatoParaGuardar: 9 cifras significativas");
}


// ---- Node-splitting con un estudio de tres brazos (red139: 11 estudios, 4 tratamientos, S8 de tres brazos) ----
// Error 5: la evidencia directa de cada par se agrupaba con las varianzas YA ajustadas por multibrazo; `netmeta` usa las
// originales (y P = Var(red)/Var(directa)). Sin estudios multibrazo ambas coinciden, por eso solo se veía con ellos.
{
  const filas: [string, string, number, number, number][] = [["S1", "B", 119, 9.755, 1.3], ["S1", "A", 31, 9.515, 1.732], ["S2", "B", 50, 9.908, 3.936], ["S2", "C", 190, 10.762, 3.422], ["S3", "C", 200, 11.076, 1.292], ["S3", "B", 115, 9.826, 3.843], ["S4", "A", 171, 9.622, 3.259], ["S4", "D", 100, 10.17, 3.354], ["S5", "B", 175, 8.949, 3.009], ["S5", "A", 104, 9.856, 3.678], ["S6", "B", 188, 9.834, 1.375], ["S6", "D", 25, 11.368, 1.951], ["S7", "D", 98, 10.531, 3.529], ["S7", "B", 92, 10.023, 1.382], ["S8", "C", 68, 12.122, 3.443], ["S8", "B", 33, 10.128, 3.744], ["S8", "A", 36, 10.789, 3.967], ["S9", "D", 196, 10.325, 1.678], ["S9", "B", 136, 8.928, 1.981], ["S10", "B", 119, 10.357, 2.097], ["S10", "C", 111, 11.887, 1.667], ["S11", "C", 127, 10.502, 3.736], ["S11", "B", 193, 9.648, 3.098]];
  const por = new Map<string, { tratamiento: string; media: string; de: string; n: string }[]>();
  for (const [s, t, n, media, de] of filas) por.set(s, [...(por.get(s) ?? []), { tratamiento: t, media: String(media), de: String(de), n: String(n) }]);
  const red = [...por.entries()].map(([etiqueta, brazos], i) => ({ id: i + 1, etiqueta, medida: "DM" as const, brazos }));
  const r = calcularNMA(red, "DM", "aleatorio")!;
  const ns = r.nodeSplitting.find((s) => s.tratamiento1 === "A" && s.tratamiento2 === "D");
  ok(!!ns, "Node-splitting: existe la comparación A–D (netmeta la etiqueta D:A, con el signo opuesto)");
  ok(!!ns && cerca(ns.diffDirecto, -0.548, 1e-9), "Node-splitting (multibrazo): estimación directa = netmeta");
  ok(!!ns && cerca(ns.estimacionIndirecta, -0.963629322161634, 1e-9), "Node-splitting (multibrazo): estimación indirecta = netmeta");
  ok(!!ns && cerca(ns.proporcionDirecta, 0.377477072300471, 1e-9), "Node-splitting (multibrazo): proporción de evidencia directa = netmeta");
  ok(!!ns && cerca(ns.pValor, 0.515293642011421, 1e-8), "Node-splitting (multibrazo): p de la diferencia = netmeta");
  // B–C: los dos tratamientos están en el estudio de tres brazos S8, así que su evidencia directa mezcla contrastes ajustados.
  const bc = r.nodeSplitting.find((s) => s.tratamiento1 === "B" && s.tratamiento2 === "C");
  ok(!!bc && cerca(bc.diffDirecto, -1.27767659723935, 1e-9) && cerca(bc.estimacionIndirecta, -1.53783630137483, 1e-9) && cerca(bc.proporcionDirecta, 0.967147527399579, 1e-9),
    "Node-splitting B–C (par del estudio de tres brazos): directa, indirecta y proporción = netmeta");
  const si = calcularNodeSplittingSIDDE(red, "DM", "aleatorio").find((s) => s.tratamiento1 === "A" && s.tratamiento2 === "D");
  ok(!!si && cerca(si.estimacionIndirecta, -0.956602470077041, 1e-9) && cerca(si.pValor, 0.526714770176988, 1e-8), "Node-splitting SIDDE (multibrazo) = netmeta");
}

// ---- Cola normal exacta (valores de R: 2 * pnorm(-z)) ----
// Error 6: la cola normal usaba la aproximación de Abramowitz-Stegun (error absoluto hasta 1,5e-7); afectaba al p de los
// coeficientes de la meta-regresión, de la pendiente de Moses-Littenberg, del node-splitting y de κ.
{
  const casos: [number, number][] = [[0.1, 0.920344325445942], [0.5, 0.617075077451974], [1.96, 0.0499957902964409], [2.5, 0.0124193306515523], [4, 6.33424836662398e-05], [5, 5.73303143758388e-07]];
  ok(casos.every(([z, p]) => cerca(colaNormalDosLados(z), p, 1e-6, 1e-12)), "Cola normal de dos lados = R hasta 1e-6 relativo (antes ~1e-7 absoluto)");
}

console.log(fallos === 0 ? "\nTodo bien." : `\n${fallos} fallo(s).`);
process.exit(fallos === 0 ? 0 : 1);
