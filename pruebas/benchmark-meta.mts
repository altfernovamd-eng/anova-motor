// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Prueba de regresión: valida el motor estadístico de aNova (`src/lib/meta-analysis/statistics.ts`)
 * contra el paquete `meta` de R (el estándar de referencia mundial), reproduciendo su salida
 * publicada sobre un dataset canónico.
 *
 * Dataset:  ThirdWave (18 estudios) — Harrer et al., "Doing Meta-Analysis in R".
 * Referencia: metagen(TE, seTE, method.tau = "REML", hakn = TRUE) — cifras publicadas en el libro.
 *
 * Cómo correrlo:   npm run benchmark:meta      (o)   node pruebas/benchmark-meta.mts
 * Sale con código 1 si alguna cifra deja de coincidir (útil en CI / al tocar el motor).
 *
 * Por qué es válido: aNova agrupa por inverso-varianza sobre el efecto transformado; la matemática
 * de pooling es la misma sea cual sea la medida. Alimentamos TE/seTE como una diferencia de medias
 * (medida "DM", escala lineal) con IC = TE ± Z95·seTE, de modo que el SE se reconstruye EXACTO.
 */
import {
  ajustarMetaRegresion,
  calcularMetaAnalisis,
  confianzaHeterogeneidad,
  efectoDesde2x2,
  efectoDesdeMediasEstandarizada,
  metaAnalisisAcumulativo,
  metaAnalisisMH,
  metaAnalisisPeto,
  testBegg,
} from "../src/statistics.ts";

const Z95 = 1.959963985; // debe coincidir con statistics.ts para que el SE se reconstruya exacto

// --- Datos ThirdWave (TE = g de Hedges, seTE = su error estándar), precisión completa del .rda oficial ---
const THIRDWAVE: { autor: string; te: number; se: number }[] = [
  { autor: "Call et al.", te: 0.709136234756933, se: 0.260820178547769 },
  { autor: "Cavanagh et al.", te: 0.354864109077645, se: 0.196362410818884 },
  { autor: "DanitzOrsillo", te: 1.79117003272488, se: 0.345569205699001 },
  { autor: "de Vibe et al.", te: 0.182455216323761, se: 0.117787401662918 },
  { autor: "Frazier et al.", te: 0.42185093084487, se: 0.144812762024237 },
  { autor: "Frogeli et al.", te: 0.63, se: 0.196 },
  { autor: "Gallego et al.", te: 0.724883808783712, se: 0.224664137749212 },
  { autor: "Hazlett-Stevens & Oren", te: 0.528663812357321, se: 0.210460918848162 },
  { autor: "Hintz et al.", te: 0.284, se: 0.168 },
  { autor: "Kang et al.", te: 1.27506817139631, se: 0.337199688565405 },
  { autor: "Kuhlmann et al.", te: 0.10360817637506, se: 0.194727463484504 },
  { autor: "Lever Taylor et al.", te: 0.388390560376566, se: 0.230768893346391 },
  { autor: "Phang et al.", te: 0.540739823331358, se: 0.244313267301974 },
  { autor: "Rasanen et al.", te: 0.426159311552595, se: 0.257937944593587 },
  { autor: "Ratanasiripong", te: 0.515396915332498, se: 0.351273746997289 },
  { autor: "Shapiro et al.", te: 1.47972601084922, se: 0.315281680355345 },
  { autor: "Song & Lindquist", te: 0.612578239930467, se: 0.226683364977455 },
  { autor: "Warnecke et al.", te: 0.6, se: 0.249 },
];

// --- Salida publicada de `meta` (metagen, REML + Knapp-Hartung) para ThirdWave ---
const META = {
  efecto: 0.5771,
  hkInf: 0.3782,
  hkSup: 0.776,
  t: 6.12,
  tau2: 0.082,
  tau: 0.2863,
  i2: 62.6,
  q: 45.5,
  gl: 17,
  pQ: 0.0002,
};

const estudios = THIRDWAVE.map((d, i) => ({
  id: i + 1,
  autor: d.autor,
  año: "",
  diseño: "ECA",
  n: "",
  medida: "DM",
  valor: String(d.te),
  icInferior: String(d.te - Z95 * d.se),
  icSuperior: String(d.te + Z95 * d.se),
}));

const r = calcularMetaAnalisis(estudios as never, "DM", "REML");
if (!r) {
  console.error("FALLO: calcularMetaAnalisis devolvió null");
  process.exit(1);
}
const h = r.heterogeneidad;
const a = r.aleatorio;
const tCrit17 = 2.109815578; // t_{17, 0.975}
const seHK = (a.hkIcSuperior - a.hkIcInferior) / (2 * tCrit17);
const tStat = a.estimacion / seHK;

// Cada fila: valor de aNova, valor publicado por meta, y decimales a los que debe coincidir.
const filas: { nombre: string; anova: number; meta: number; dp: number }[] = [
  { nombre: "Efecto aleatorio (SMD)", anova: a.estimacion, meta: META.efecto, dp: 4 },
  { nombre: "IC 95% HK inferior", anova: a.hkIcInferior, meta: META.hkInf, dp: 4 },
  { nombre: "IC 95% HK superior", anova: a.hkIcSuperior, meta: META.hkSup, dp: 4 },
  { nombre: "Estadístico t (gl=17)", anova: tStat, meta: META.t, dp: 2 },
  { nombre: "tau^2 (REML)", anova: h.tau2, meta: META.tau2, dp: 4 },
  { nombre: "tau (REML)", anova: h.tau, meta: META.tau, dp: 4 },
  { nombre: "I^2 (%)", anova: h.i2, meta: META.i2, dp: 1 },
  { nombre: "Q de Cochran", anova: h.q, meta: META.q, dp: 2 },
  { nombre: "gl", anova: h.gl, meta: META.gl, dp: 0 },
  { nombre: "p-valor de Q", anova: h.pValor, meta: META.pQ, dp: 4 },
];

// IC Q-profile de τ² y τ (metagen los reporta por defecto con este método).
const cih = confianzaHeterogeneidad(r)!;
filas.push(
  { nombre: "IC τ² inferior (Q-profile)", anova: cih.tau2Inferior, meta: 0.0295, dp: 4 },
  { nombre: "IC τ² superior (Q-profile)", anova: cih.tau2Superior, meta: 0.3533, dp: 4 },
  { nombre: "IC τ inferior (Q-profile)", anova: cih.tauInferior, meta: 0.1717, dp: 4 },
  { nombre: "IC τ superior (Q-profile)", anova: cih.tauSuperior, meta: 0.5944, dp: 4 },
);

let fallos = 0;
console.log(`\naNova  vs  R 'meta' (metagen, REML + Knapp-Hartung)  —  ThirdWave, k=${r.estudios.length}\n`);
console.log(`${"Cantidad".padEnd(26)} ${"aNova".padEnd(12)} ${"meta".padEnd(12)} Estado`);
console.log("-".repeat(62));
for (const f of filas) {
  const anovaR = f.anova.toFixed(f.dp);
  const metaR = f.meta.toFixed(f.dp);
  const ok = anovaR === metaR;
  if (!ok) fallos++;
  console.log(`${f.nombre.padEnd(26)} ${anovaR.padEnd(12)} ${metaR.padEnd(12)} ${ok ? "✓" : "✗ DIFIERE"}`);
}
console.log("-".repeat(62));
let totalChecks = filas.length;

// ===========================================================================
//  Benchmark 2 — Mantel-Haenszel (datos binarios 2×2) vs metafor::rma.mh
//  Dataset: dat.bcg (13 estudios, vacuna BCG; sin celdas cero) — paquete metadat.
// ===========================================================================
// tpos, tneg, cpos, cneg
const BCG: number[][] = [
  [4, 119, 11, 128], [6, 300, 29, 274], [3, 228, 11, 209], [62, 13536, 248, 12619],
  [33, 5036, 47, 5761], [180, 1361, 372, 1079], [8, 2537, 10, 619], [505, 87886, 499, 87892],
  [29, 7470, 45, 7232], [17, 1699, 65, 1600], [186, 50448, 141, 27197], [5, 2493, 3, 2338],
  [27, 16886, 29, 17825],
];

// Salida publicada de metafor::rma.mh (por defecto) sobre dat.bcg — verbatim del sitio metafor.
const RMA_MH = {
  RR: { logEst: -0.4537, se: 0.0393, z: -11.5338, logInf: -0.5308, logSup: -0.3766, est: 0.6353, inf: 0.5881, sup: 0.6862, q: 152.5676, i2: 92.13 },
  OR: { logEst: -0.4734, se: 0.041, z: -11.5444, logInf: -0.5538, logSup: -0.393, est: 0.6229, inf: 0.5748, sup: 0.675, q: 163.9426, i2: 92.68 },
};

function bcgEstudios(medida: string) {
  return BCG.map(([tp, tn, cp, cn], i) => ({
    id: i + 1, autor: "s" + i, año: "", diseño: "ECA", n: "", medida, valor: "", icInferior: "", icSuperior: "",
    crudos: { tipo: "2x2", eventosTrat: String(tp), totalTrat: String(tp + tn), eventosCtrl: String(cp), totalCtrl: String(cp + cn) },
  }));
}

for (const medida of ["RR", "OR"] as const) {
  const r = metaAnalisisMH(bcgEstudios(medida) as never, medida as never)!;
  const ref = RMA_MH[medida];
  const logEst = Math.log(r.estimacion);
  const filasMH: { nombre: string; anova: number; meta: number; dp: number }[] = [
    { nombre: "log-estimación", anova: logEst, meta: ref.logEst, dp: 4 },
    { nombre: "SE (log)", anova: r.se, meta: ref.se, dp: 4 },
    { nombre: "z", anova: logEst / r.se, meta: ref.z, dp: 4 },
    { nombre: "IC log inferior", anova: logEst - Z95 * r.se, meta: ref.logInf, dp: 4 },
    { nombre: "IC log superior", anova: logEst + Z95 * r.se, meta: ref.logSup, dp: 4 },
    { nombre: `${medida} (escala natural)`, anova: r.estimacion, meta: ref.est, dp: 4 },
    { nombre: `IC ${medida} inferior`, anova: r.icInferior, meta: ref.inf, dp: 4 },
    { nombre: `IC ${medida} superior`, anova: r.icSuperior, meta: ref.sup, dp: 4 },
    { nombre: "Q de heterogeneidad", anova: r.heterogeneidad.q, meta: ref.q, dp: 4 },
    { nombre: "I^2 (%)", anova: r.heterogeneidad.i2, meta: ref.i2, dp: 2 },
  ];
  console.log(`\naNova  vs  metafor::rma.mh (Mantel-Haenszel, ${medida})  —  dat.bcg, k=${r.k}\n`);
  console.log(`${"Cantidad".padEnd(26)} ${"aNova".padEnd(12)} ${"metafor".padEnd(12)} Estado`);
  console.log("-".repeat(62));
  for (const f of filasMH) {
    const anovaR = f.anova.toFixed(f.dp);
    const metaR = f.meta.toFixed(f.dp);
    const ok = anovaR === metaR;
    if (!ok) fallos++;
    totalChecks++;
    console.log(`${f.nombre.padEnd(26)} ${anovaR.padEnd(12)} ${metaR.padEnd(12)} ${ok ? "✓" : "✗ DIFIERE"}`);
  }
  console.log("-".repeat(62));
}

// --- Benchmark 3 — Peto (OR de un paso) vs metafor::rma.peto, sobre dat.bcg ---
{
  const r = metaAnalisisPeto(bcgEstudios("OR") as never)!;
  const logEst = Math.log(r.estimacion);
  const ref = { logEst: -0.4744, se: 0.0407, z: -11.6689, logInf: -0.5541, logSup: -0.3948, est: 0.6222, inf: 0.5746, sup: 0.6738, q: 167.7302, i2: 92.85 };
  const filasP: { nombre: string; anova: number; meta: number; dp: number }[] = [
    { nombre: "log-estimación", anova: logEst, meta: ref.logEst, dp: 4 },
    { nombre: "SE (log)", anova: r.se, meta: ref.se, dp: 4 },
    { nombre: "z", anova: logEst / r.se, meta: ref.z, dp: 4 },
    { nombre: "IC log inferior", anova: logEst - Z95 * r.se, meta: ref.logInf, dp: 4 },
    { nombre: "IC log superior", anova: logEst + Z95 * r.se, meta: ref.logSup, dp: 4 },
    { nombre: "OR (escala natural)", anova: r.estimacion, meta: ref.est, dp: 4 },
    { nombre: "IC OR inferior", anova: r.icInferior, meta: ref.inf, dp: 4 },
    { nombre: "IC OR superior", anova: r.icSuperior, meta: ref.sup, dp: 4 },
    { nombre: "Q de heterogeneidad", anova: r.heterogeneidad.q, meta: ref.q, dp: 4 },
    { nombre: "I^2 (%)", anova: r.heterogeneidad.i2, meta: ref.i2, dp: 2 },
  ];
  console.log(`\naNova  vs  metafor::rma.peto (método de Peto, OR)  —  dat.bcg, k=${r.k}\n`);
  console.log(`${"Cantidad".padEnd(26)} ${"aNova".padEnd(12)} ${"metafor".padEnd(12)} Estado`);
  console.log("-".repeat(62));
  for (const f of filasP) {
    const anovaR = f.anova.toFixed(f.dp);
    const metaR = f.meta.toFixed(f.dp);
    const ok = anovaR === metaR;
    if (!ok) fallos++;
    totalChecks++;
    console.log(`${f.nombre.padEnd(26)} ${anovaR.padEnd(12)} ${metaR.padEnd(12)} ${ok ? "✓" : "✗ DIFIERE"}`);
  }
  console.log("-".repeat(62));
}

// --- Benchmark 4 — g de Hedges (SMD) vs metafor::escalc(measure="SMD"), sobre dat.normand1999 ---
// Cada fila: [n1, m1, sd1, n2, m2, sd2]  →  [g, v] de referencia (fórmula exacta de escalc, ec. Hedges 1982).
{
  const NORMAND: number[][] = [
    [155, 55, 47, 156, 75, 64], [31, 27, 7, 32, 29, 4], [75, 64, 17, 71, 119, 29],
    [18, 66, 20, 18, 137, 48], [8, 14, 8, 13, 18, 11], [57, 19, 7, 52, 18, 4],
    [34, 52, 45, 33, 41, 34], [110, 21, 16, 183, 31, 27], [60, 30, 27, 52, 23, 20],
  ];
  const REF: number[][] = [
    [-0.3552, 0.0131], [-0.3479, 0.0645], [-2.3176, 0.0458], [-1.888, 0.1606], [-0.384, 0.2054],
    [0.1721, 0.0369], [0.2721, 0.0603], [-0.4246, 0.0149], [0.2896, 0.0363],
  ];
  console.log(`\naNova  vs  metafor::escalc(measure="SMD") — g de Hedges  —  dat.normand1999, k=${NORMAND.length}\n`);
  console.log(`${"Estudio".padEnd(10)} ${"g aNova".padEnd(10)} ${"g meta".padEnd(10)} ${"v aNova".padEnd(10)} ${"v meta".padEnd(10)} Estado`);
  console.log("-".repeat(62));
  NORMAND.forEach(([n1, m1, s1, n2, m2, s2], i) => {
    const e = efectoDesdeMediasEstandarizada(m1, s1, n1, m2, s2, n2)!;
    const v = e.se * e.se;
    const gOk = e.valor.toFixed(4) === REF[i][0].toFixed(4);
    const vOk = v.toFixed(4) === REF[i][1].toFixed(4);
    if (!gOk) fallos++;
    if (!vOk) fallos++;
    totalChecks += 2;
    console.log(
      `${String(i + 1).padEnd(10)} ${e.valor.toFixed(4).padEnd(10)} ${REF[i][0].toFixed(4).padEnd(10)} ${v.toFixed(4).padEnd(10)} ${REF[i][1].toFixed(4).padEnd(10)} ${gOk && vOk ? "✓" : "✗ DIFIERE"}`,
    );
  });
  console.log("-".repeat(62));
}

// --- Benchmark 5 — Meta-análisis acumulativo vs metafor::cumul (RR, aleatorio REML, orden por año), dat.bcg ---
// Se valida la EVOLUCIÓN del efecto (log-estimación + IC), que es la sustancia del acumulativo.
{
  // [autor, año, a, tneg, c, cneg]
  const BCG_YR: [string, string, number, number, number, number][] = [
    ["Aronson", "1948", 4, 119, 11, 128], ["Ferguson & Simes", "1949", 6, 300, 29, 274],
    ["Rosenthal", "1960", 3, 228, 11, 209], ["Hart & Sutherland", "1977", 62, 13536, 248, 12619],
    ["Frimodt-Moller", "1973", 33, 5036, 47, 5761], ["Stein & Aronson", "1953", 180, 1361, 372, 1079],
    ["Vandiviere", "1973", 8, 2537, 10, 619], ["TPT Madras", "1980", 505, 87886, 499, 87892],
    ["Coetzee & Berjak", "1968", 29, 7470, 45, 7232], ["Rosenthal", "1961", 17, 1699, 65, 1600],
    ["Comstock", "1974", 186, 50448, 141, 27197], ["Comstock & Webster", "1969", 5, 2493, 3, 2338],
    ["Comstock", "1976", 27, 16886, 29, 17825],
  ];
  // metafor cumul (log-escala): [logEst, ci.lb, ci.ub] por paso acumulado (orden por año)
  const REFC: [number, number, number][] = [
    [-0.889, -2.008, 0.229], [-1.325, -2.009, -0.641], [-0.972, -1.453, -0.491], [-1.001, -1.437, -0.565],
    [-1.101, -1.493, -0.710], [-0.973, -1.341, -0.606], [-0.901, -1.275, -0.527], [-0.789, -1.175, -0.403],
    [-0.866, -1.250, -0.481], [-0.785, -1.138, -0.432], [-0.708, -1.058, -0.359], [-0.795, -1.147, -0.442],
    [-0.715, -1.067, -0.362],
  ];
  const est = BCG_YR.map(([autor, año, a, tn, c, cn], i) => {
    const ef = efectoDesde2x2(a, a + tn, c, c + cn, "RR")!;
    return { id: i + 1, autor, año, diseño: "ECA", n: "", medida: "RR", valor: String(ef.valor), icInferior: String(ef.icInferior), icSuperior: String(ef.icSuperior) };
  });
  const filas = metaAnalisisAcumulativo(est as never, "RR", "aleatorio", "REML")!;
  console.log(`\naNova  vs  metafor::cumul (acumulativo, RR aleatorio REML, orden por año)  —  dat.bcg, k=${filas.length}\n`);
  console.log(`${"paso".padEnd(6)} ${"log-est".padEnd(18)} ${"IC (log)".padEnd(26)} Estado`);
  console.log("-".repeat(62));
  filas.forEach((f, i) => {
    const le = Math.log(f.estimacion), li = Math.log(f.icInferior), ls = Math.log(f.icSuperior);
    const ok = le.toFixed(3) === REFC[i][0].toFixed(3) && li.toFixed(3) === REFC[i][1].toFixed(3) && ls.toFixed(3) === REFC[i][2].toFixed(3);
    if (!ok) fallos++;
    totalChecks += 3;
    console.log(`${String(f.k).padEnd(6)} ${`${le.toFixed(3)} (${REFC[i][0].toFixed(3)})`.padEnd(18)} ${`[${li.toFixed(3)};${ls.toFixed(3)}] vs [${REFC[i][1].toFixed(3)};${REFC[i][2].toFixed(3)}]`.padEnd(26)} ${ok ? "✓" : "✗ DIFIERE"}`);
  });
  console.log("-".repeat(62));
}

// --- Benchmark 6 — Test de Begg (correlación de rangos) vs metafor::ranktest, sobre dat.bcg (RR) ---
{
  const est = BCG.map(([tp, tn, cp, cn], i) => {
    const ef = efectoDesde2x2(tp, tp + tn, cp, cp + cn, "RR")!;
    return { id: i + 1, autor: "s" + i, año: "", diseño: "ECA", n: "", medida: "RR", valor: String(ef.valor), icInferior: String(ef.icInferior), icSuperior: String(ef.icSuperior) };
  });
  const r = calcularMetaAnalisis(est as never, "RR", "DL")!;
  const b = testBegg(r)!;
  console.log(`\naNova  vs  metafor::ranktest (test de Begg, RR)  —  dat.bcg, k=${b.k}\n`);
  console.log(`${"Cantidad".padEnd(26)} ${"aNova".padEnd(12)} ${"metafor".padEnd(12)} Estado`);
  console.log("-".repeat(62));
  for (const [nombre, anova, meta, dp] of [
    ["Kendall's tau", b.tau, 0.0256, 4],
    ["p-valor (exacta)", b.pValor, 0.9524, 4],
  ] as [string, number, number, number][]) {
    const a = anova.toFixed(dp);
    const m = meta.toFixed(dp);
    const ok = a === m;
    if (!ok) fallos++;
    totalChecks++;
    console.log(`${nombre.padEnd(26)} ${a.padEnd(12)} ${m.padEnd(12)} ${ok ? "✓" : "✗ DIFIERE"}`);
  }
  console.log("-".repeat(62));
}

// --- Benchmark 7 — Meta-regresión (mixta, REML) vs metafor::rma(mods=…), dat.bcg (RR) ---
// Referencia REML generada de forma independiente (numpy), anclada al τ² canónico 0.3132 de dat.bcg.
{
  const ABLAT = [44, 55, 42, 52, 13, 44, 19, 13, 27, 42, 18, 33, 33];
  const ALLOC = ["random", "random", "random", "random", "alternate", "alternate", "random", "random", "random", "systematic", "systematic", "systematic", "systematic"];
  const pts = BCG.map(([tp, tn, cp, cn]) => {
    const a = tp, n1 = tp + tn, c = cp, n2 = cp + cn;
    return { y: Math.log((a / n1) / (c / n2)), se: Math.sqrt(1 / a - 1 / n1 + 1 / c - 1 / n2) };
  });
  const mediaAb = ABLAT.reduce((s, x) => s + x, 0) / ABLAT.length;

  const chk = (nombre: string, anova: number, meta: number, dp: number) => {
    const a = anova.toFixed(dp);
    const m = meta.toFixed(dp);
    const ok = a === m;
    if (!ok) fallos++;
    totalChecks++;
    console.log(`${nombre.padEnd(30)} ${a.padEnd(11)} ${m.padEnd(11)} ${ok ? "✓" : "✗ DIFIERE"}`);
  };

  console.log(`\naNova  vs  metafor::rma (meta-regresión REML)  —  dat.bcg (RR)\n`);
  console.log(`${"Cantidad".padEnd(30)} ${"aNova".padEnd(11)} ${"metafor".padEnd(11)} Estado`);
  console.log("-".repeat(64));
  // continuo: ablat (centrado)
  const ra = ajustarMetaRegresion(pts, BCG.map((_, i) => [1, ABLAT[i] - mediaAb]), ["Intercepto", "ablat"])!;
  chk("ablat · τ² residual", ra.tau2Residual, 0.0763, 4);
  chk("ablat · QM", ra.qm, 16.3582, 4);
  chk("ablat · R²", ra.r2, 75.63, 2);
  chk("ablat · pendiente", ra.coeficientes[1].estimacion, -0.0291, 4);
  chk("ablat · SE(pendiente)", ra.coeficientes[1].se, 0.0072, 4);
  chk("ablat · intercepto", ra.coeficientes[0].estimacion, -0.7223, 4);
  // categórico: alloc (ref = alternate)
  const rc = ajustarMetaRegresion(
    pts,
    BCG.map((_, i) => [1, ALLOC[i] === "random" ? 1 : 0, ALLOC[i] === "systematic" ? 1 : 0]),
    ["Intercepto", "random", "systematic"],
  )!;
  chk("alloc · τ² residual", rc.tau2Residual, 0.3615, 4);
  chk("alloc · QM (df=2)", rc.qm, 1.7675, 4);
  chk("alloc · p de QM", rc.qmPValor, 0.4132, 4);
  chk("alloc · coef. random", rc.coeficientes[1].estimacion, -0.4478, 4);
  chk("alloc · coef. systematic", rc.coeficientes[2].estimacion, 0.089, 4);
  console.log("-".repeat(64));
}

if (fallos === 0) {
  console.log(`\n✓ PASA: las ${totalChecks} cantidades coinciden con R (meta / metafor) a los decimales publicados.\n`);
  process.exit(0);
} else {
  console.error(`\n✗ FALLA: ${fallos} de ${totalChecks} cantidad(es) ya no coinciden con R. Revisa el motor estadístico.\n`);
  process.exit(1);
}
