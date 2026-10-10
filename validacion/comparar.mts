// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Validación DIFERENCIAL del motor estadístico de aNova contra R (`metafor` y `meta`).
 *
 * Toma los cientos de conjuntos de datos sintéticos de `referencias/referencia.json` (generados por
 * `generar-referencia.R` con semilla fija), los pasa por el motor REAL de aNova y compara, número a número, con lo que
 * calculó R. No necesita R para correr: R solo se usa para regenerar la referencia.
 *
 * Uso:  npm run validacion         (o)   npx tsx validacion/comparar.mts [--detalle]
 * Sale con código 1 si alguna cantidad se aparta de R más que su tolerancia (y no es una diferencia documentada).
 */
import { existsSync, readFileSync } from "node:fs";
import {
  analisisSensibilidad,
  baujat,
  calcularMetaAnalisis,
  confianzaHeterogeneidad,
  diagnosticosInfluencia,
  galbraith,
  analisisSubgrupos,
  efectoDesde2x2,
  efectoDesdeCorrelacion,
  efectoDesdeMedias,
  efectoDesdeMediasEstandarizada,
  efectoDesdeOEV,
  efectoDesdeProporcion,
  efectoDesdeRazonDeMedias,
  efectoDesdeTasas,
  formatoParaGuardar,
  metaAnalisisAcumulativo,
  metaAnalisisMH,
  metaAnalisisPeto,
  testBegg,
  testEgger,
  testPeters,
  transformar,
  trimAndFill,
  type EscalaAnalisis,
  type MetodoTau2,
} from "../src/statistics.ts";
import type { Estudio, Medida } from "../src/tipos.ts";
import { compararAvanzada } from "./comparar-avanzada.mts";
import { compararBivariado } from "./comparar-bivariado.mts";
import { compararMultinivel } from "./comparar-multinivel.mts";
import { compararGLMM } from "./comparar-glmm.mts";

const Z95 = 1.959963985; // igual que statistics.ts: el SE se reconstruye exacto desde el IC
const DETALLE = process.argv.includes("--detalle");
// eslint-disable-next-line @typescript-eslint/no-explicit-any
type J = any;
const ref: J = JSON.parse(readFileSync(new URL("./referencias/referencia.json", import.meta.url), "utf8"));
// Con VALIDACION_AVANZADA=masiva (`npm run validacion:masiva`) la parte de SROC y red usa la referencia masiva (400 SROC y 300 redes, con
// 4 brazos y celdas en cero, más `mada::madauni`) y, si existe, la de celdas en cero (`referencia-red-ceros.json`: estudios multibrazo con
// doble cero y eventos = n), en lugar de la normal; el resto de la validación corre igual.
const MASIVA = process.env.VALIDACION_AVANZADA === "masiva";
const refAv: J = JSON.parse(readFileSync(new URL(`./referencias/${MASIVA ? "referencia-red-sroc-masiva" : "referencia-avanzada"}.json`, import.meta.url), "utf8"));
if (MASIVA && existsSync(new URL("./referencias/referencia-red-ceros.json", import.meta.url))) {
  const ceros: J = JSON.parse(readFileSync(new URL("./referencias/referencia-red-ceros.json", import.meta.url), "utf8"));
  refAv.red = [...(refAv.red ?? []), ...(ceros.red ?? [])];
  refAv.meta.escenarios.red = refAv.red.length;
}

// ---------------------------------------------------------------------------------------------
// Registro de comparaciones
// ---------------------------------------------------------------------------------------------
interface Grupo {
  n: number;
  maxAbs: number;
  maxRel: number;
  peor: string;
  fallos: string[];
  tolRel: number;
}
const grupos = new Map<string, Grupo>();
const sinPar: string[] = [];
const convencionRAR = { n: 0, seRel: 0, qDistintas: 0 };

/** Compara un número de aNova con el de R. Tolerancia: relativa o absoluta (la que sea más holgada). */
function cmp(grupo: string, id: string, a: number | null | undefined, r: number | null | undefined, tolRel = 1e-6, tolAbs = 1e-8) {
  const g = grupos.get(grupo) ?? { n: 0, maxAbs: 0, maxRel: 0, peor: "", fallos: [], tolRel };
  grupos.set(grupo, g);
  const aNull = a === null || a === undefined || Number.isNaN(a);
  const rNull = r === null || r === undefined || Number.isNaN(r);
  if (aNull && rNull) return;
  if (aNull !== rNull) {
    g.n++;
    g.fallos.push(`${id}: aNova=${aNull ? "—" : a} R=${rNull ? "—" : r}`);
    return;
  }
  const x = a as number;
  const y = r as number;
  g.n++;
  const abs = Math.abs(x - y);
  const rel = abs / Math.max(Math.abs(x), Math.abs(y), 1e-300);
  if (abs > g.maxAbs) g.maxAbs = abs;
  if (rel > g.maxRel && abs > tolAbs) {
    g.maxRel = rel;
    g.peor = `${id} (aNova=${x}, R=${y})`;
  }
  if (abs > tolAbs && rel > tolRel) g.fallos.push(`${id}: aNova=${x} R=${y} (abs ${abs.toExponential(2)}, rel ${rel.toExponential(2)})`);
}

const arr = <T = number,>(v: J): T[] => (v === null || v === undefined ? [] : Array.isArray(v) ? v : [v]);

// ---------------------------------------------------------------------------------------------
// Construcción de estudios para el motor
// ---------------------------------------------------------------------------------------------
let contador = 0;
function estudioDesde(medida: Medida, valor: number, icInf: number, icSup: number, extra: Partial<Estudio> = {}): Estudio {
  contador++;
  return {
    id: contador, autor: `E${contador}`, año: "", diseño: "ECA", n: "", medida,
    valor: String(valor), icInferior: String(icInf), icSuperior: String(icSup), ...extra,
  } as Estudio;
}

/** Estudios a partir de efectos en escala de análisis (yi, sei): el IC se arma con Z95 para reconstruir el SE exacto. */
function estudiosDesdeYSe(medida: Medida, escala: EscalaAnalisis, yi: number[], sei: number[], anio?: number[], grupo?: string[]): Estudio[] {
  const des = (t: number) => (escala === "log" ? Math.exp(t) : escala === "logit" ? 1 / (1 + Math.exp(-t)) : t);
  return yi.map((y, i) =>
    estudioDesde(medida, des(y), des(y - Z95 * sei[i]), des(y + Z95 * sei[i]), {
      ...(anio ? { año: String(anio[i]) } : {}),
      ...(grupo ? { diseño: grupo[i] as Estudio["diseño"] } : {}),
    }),
  );
}

const ESCALA: Record<string, EscalaAnalisis> = { DM: "lineal", DME: "lineal", RAR: "lineal", OR: "log", RR: "log", Proporción: "logit", RTI: "log", DTI: "lineal", RM: "log", r: "fisher", HR: "log" };

// ---------------------------------------------------------------------------------------------
// Pooling: calcularMetaAnalisis contra metafor/meta
// ---------------------------------------------------------------------------------------------
function compararPool(prefijo: string, id: string, estudios: Estudio[], medida: Medida, pool: J) {
  const escala = ESCALA[medida];
  const t = (v: number | null | undefined) => (v === null || v === undefined ? null : transformar(v, escala));
  for (const metodo of ["DL", "PM", "REML"] as MetodoTau2[]) {
    const R = pool[metodo];
    const a = calcularMetaAnalisis(estudios, medida, metodo);
    if (!R || !a) {
      if (!!R !== !!a) sinPar.push(`${prefijo} ${id} ${metodo}: aNova=${a ? "resultado" : "null"} R=${R ? "resultado" : "sin resultado"}`);
      continue;
    }
    const g = `${prefijo} · τ² ${metodo}`;
    cmp(`${g} · τ²`, id, a.heterogeneidad.tau2, R.tau2, 1e-6, 1e-9);
    cmp(`${g} · efecto aleatorio`, id, t(a.aleatorio.estimacion), R.est);
    cmp(`${g} · IC aleatorio (z)`, id, t(a.aleatorio.icInferior), R.lb);
    cmp(`${g} · IC aleatorio (z)`, id, t(a.aleatorio.icSuperior), R.ub);
    cmp(`${g} · IC Knapp-Hartung`, id, t(a.aleatorio.hkIcInferior), R.hk_lb);
    cmp(`${g} · IC Knapp-Hartung`, id, t(a.aleatorio.hkIcSuperior), R.hk_ub);
    cmp(`${g} · intervalo de predicción`, id, t(a.aleatorio.prediccionInferior), R.pred_lb);
    cmp(`${g} · intervalo de predicción`, id, t(a.aleatorio.prediccionSuperior), R.pred_ub);
  }
  const a = calcularMetaAnalisis(estudios, medida, "DL");
  if (a) {
    cmp(`${prefijo} · efecto fijo`, id, t(a.fijo.estimacion), pool.fijo.est);
    cmp(`${prefijo} · IC fijo`, id, t(a.fijo.icInferior), pool.fijo.lb);
    cmp(`${prefijo} · IC fijo`, id, t(a.fijo.icSuperior), pool.fijo.ub);
    cmp(`${prefijo} · Q`, id, a.heterogeneidad.q, pool.q);
    cmp(`${prefijo} · p de Q`, id, a.heterogeneidad.pValor, pool.pq);
    cmp(`${prefijo} · I² (Higgins)`, id, a.heterogeneidad.i2, pool.i2, 1e-6, 1e-6);
    const ic = confianzaHeterogeneidad(calcularMetaAnalisis(estudios, medida, "REML")!);
    if (ic && pool.tau2_ic) {
      cmp(`${prefijo} · IC de τ² (Q-profile)`, id, ic.tau2Inferior, pool.tau2_ic[0], 1e-6, 1e-9);
      cmp(`${prefijo} · IC de τ² (Q-profile)`, id, ic.tau2Superior, pool.tau2_ic[1], 1e-6, 1e-9);
      cmp(`${prefijo} · IC de I² (Q-profile)`, id, ic.i2Inferior, pool.i2_ic?.[0], 1e-6, 1e-6);
      cmp(`${prefijo} · IC de I² (Q-profile)`, id, ic.i2Superior, pool.i2_ic?.[1], 1e-6, 1e-6);
    }
  }
}

// ---------------------------------------------------------------------------------------------
// Sesgo, influencia, acumulado, subgrupos: contra metafor/meta
// ---------------------------------------------------------------------------------------------
function compararAdicionales(prefijo: string, id: string, estudios: Estudio[], medida: Medida, R: J) {
  const escala = ESCALA[medida];
  const t = (v: number) => transformar(v, escala);
  const fe = calcularMetaAnalisis(estudios, medida, "DL");
  if (!fe) return;

  const eg = testEgger(fe);
  if (R.egger && eg) {
    cmp(`${prefijo} · Egger`, id + " intercepto", eg.intercepto, R.egger.intercepto, 1e-6, 1e-8);
    cmp(`${prefijo} · Egger`, id + " t", eg.t, R.egger.t, 1e-6, 1e-8);
    cmp(`${prefijo} · Egger`, id + " p", eg.pValor, R.egger.p, 1e-5, 1e-8);
  } else if (!!R.egger !== !!eg) sinPar.push(`${prefijo} ${id} Egger: aNova=${eg ? "resultado" : "null"} R=${R.egger ? "resultado" : "sin resultado"}`);

  const bg = testBegg(fe);
  if (R.begg && bg) {
    cmp(`${prefijo} · Begg`, id + " τ", bg.tau, R.begg.tau, 1e-6, 1e-8);
    if (bg.exacta) cmp(`${prefijo} · Begg (p exacta)`, id, bg.pValor, R.begg.p, 1e-5, 1e-8);
    // Con empates R usa la aproximación normal con la varianza de Kendall corregida, igual que aNova: se comparan.
    else if (bg.empates) cmp(`${prefijo} · Begg (p con empates)`, id, bg.pValor, R.begg.p, 1e-5, 1e-8);
  } else if (!!R.begg !== !!bg) sinPar.push(`${prefijo} ${id} Begg: aNova=${bg ? "resultado" : "null"} R=${R.begg ? "resultado" : "sin resultado"}`);

  const bj = baujat(fe);
  if (R.baujat && bj) {
    bj.forEach((p, i) => {
      cmp(`${prefijo} · Baujat (x)`, `${id}#${i}`, p.x, arr(R.baujat.x)[i], 1e-6, 1e-9);
      cmp(`${prefijo} · Baujat (y)`, `${id}#${i}`, p.y, arr(R.baujat.y)[i], 1e-6, 1e-9);
    });
  }

  if (R.trimfill) {
    // R exige indicar el lado donde imputar; aNova lo decide por la suma de rangos. Se compara con el mismo lado.
    for (const [modeloR, modelo, metodo] of [["FE", "fijo", "DL"], ["DL", "aleatorio", "DL"], ["REML", "aleatorio", "REML"]] as const) {
      const tf = trimAndFill(fe, modelo, metodo as MetodoTau2);
      const L = (lado: string) => R.trimfill[`${modeloR}_${lado}`];
      if (!L("left") || !L("right")) continue;
      if (!tf) {
        if (L("left").k0 > 0 || L("right").k0 > 0) sinPar.push(`${prefijo} ${id} trim-fill ${modeloR}: aNova=null R k0 izq=${L("left").k0} der=${L("right").k0}`);
        continue;
      }
      const rr = L(tf.lado === "izquierda" ? "left" : "right");
      cmp(`${prefijo} · trim-and-fill k0 (${modeloR})`, id, tf.k0, rr.k0, 0, 0);
      if (tf.k0 === rr.k0) cmp(`${prefijo} · trim-and-fill efecto ajustado (${modeloR})`, id, t(tf.estimacionAjustada), rr.est_estricto ?? rr.est, 1e-6, 1e-8);
    }
  }

  // Diagnóstico de influencia (rstudent, DFFITS, Cook, razón de covarianzas, τ² y Q sin el estudio, hat, peso, DFBETAS) y Galbraith.
  for (const [modeloR, modelo, metodo] of [["FE", "fijo", "DL"], ["DL", "aleatorio", "DL"], ["REML", "aleatorio", "REML"]] as const) {
    const RI = R.influencia?.[modeloR];
    if (!RI) continue;
    const res = calcularMetaAnalisis(estudios, medida, metodo as MetodoTau2);
    const di = res ? diagnosticosInfluencia(res, modelo, metodo as MetodoTau2) : null;
    if (!di) {
      sinPar.push(`${prefijo} ${id} influencia ${modeloR}: aNova=null R=resultado`);
      continue;
    }
    const G = `${prefijo} · influencia (${modeloR})`;
    di.estudios.forEach((d, i) => {
      const f = (x: unknown) => arr(x)[i];
      cmp(`${G} · residuo estandarizado externamente`, `${id}#${i}`, d.rstudent, f(RI.rstudent), 1e-6, 1e-8);
      cmp(`${G} · DFFITS`, `${id}#${i}`, d.dffits, f(RI.dffits), 1e-6, 1e-8);
      cmp(`${G} · distancia de Cook`, `${id}#${i}`, d.cookD, f(RI.cook), 1e-6, 1e-8);
      cmp(`${G} · razón de covarianzas`, `${id}#${i}`, d.covR, f(RI.covr), 1e-6, 1e-8);
      cmp(`${G} · τ² sin el estudio`, `${id}#${i}`, d.tau2Sin, f(RI.tau2del), 1e-6, 1e-9);
      cmp(`${G} · Q sin el estudio`, `${id}#${i}`, d.qSin, f(RI.qedel), 1e-6, 1e-8);
      cmp(`${G} · hat`, `${id}#${i}`, d.hat, f(RI.hat), 1e-6, 1e-9);
      cmp(`${G} · peso (%)`, `${id}#${i}`, d.pesoPct, f(RI.weight), 1e-6, 1e-7);
      cmp(`${G} · DFBETAS`, `${id}#${i}`, d.dfbetas, f(RI.dfbs), 1e-6, 1e-8);
      cmp(`${G} · marca de estudio influyente`, `${id}#${i}`, d.influyente ? 1 : 0, ([] as boolean[]).concat(RI.marcado)[i] ? 1 : 0, 0, 0);
    });
    const gb = res ? galbraith(res, modelo) : null;
    if (gb && RI.radial_x) {
      gb.puntos.forEach((p, i) => {
        cmp(`${G.replace("influencia", "Galbraith")} · precisión (x)`, `${id}#${i}`, p.x, arr(RI.radial_x)[i], 1e-9, 1e-12);
        cmp(`${G.replace("influencia", "Galbraith")} · efecto estandarizado (z)`, `${id}#${i}`, p.z, arr(RI.radial_z)[i], 1e-9, 1e-12);
      });
    }
  }

  for (const [modeloR, modelo, metodo] of [["FE", "fijo", "DL"], ["DL", "aleatorio", "DL"], ["REML", "aleatorio", "REML"]] as const) {
    const lo = R.leave1out?.[modeloR];
    if (!lo) continue;
    const s = analisisSensibilidad(estudios, medida, modelo, metodo as MetodoTau2);
    if (!s) continue;
    s.forEach((f, i) => {
      cmp(`${prefijo} · dejar uno fuera (${modeloR})`, `${id}#${i}`, t(f.estimacion), arr(lo.est)[i], 1e-6, 1e-8);
      cmp(`${prefijo} · dejar uno fuera, IC (${modeloR})`, `${id}#${i}`, t(f.icInferior), arr(lo.lb)[i], 1e-6, 1e-8);
    });
  }

  for (const [modeloR, modelo, metodo] of [["FE", "fijo", "DL"], ["DL", "aleatorio", "DL"], ["REML", "aleatorio", "REML"]] as const) {
    const cu = R.acumulado?.[modeloR];
    if (!cu) continue;
    const a = metaAnalisisAcumulativo(estudios, medida, modelo, metodo as MetodoTau2);
    if (!a) continue;
    a.forEach((f, i) => {
      cmp(`${prefijo} · acumulado (${modeloR})`, `${id}#${i}`, t(f.estimacion), arr(cu.est)[i], 1e-6, 1e-8);
      cmp(`${prefijo} · acumulado, IC (${modeloR})`, `${id}#${i}`, t(f.icInferior), arr(cu.lb)[i], 1e-6, 1e-8);
    });
  }

  if (R.subgrupos) {
    for (const [metodoR, metodo] of [["DL", "DL"], ["REML", "REML"]] as const) {
      const sg = R.subgrupos[metodoR];
      if (!sg) continue;
      for (const [modeloR, modelo] of [["aleatorio", "aleatorio"], ["fijo", "fijo"]] as const) {
        const a = analisisSubgrupos(estudios, medida, "diseño", modelo, metodo as MetodoTau2);
        if (!a) continue;
        const niveles: string[] = arr(sg.niveles);
        const k = arr(sg.k);
        a.subgrupos.forEach((s) => {
          const i = niveles.indexOf(s.etiqueta);
          if (i < 0) { sinPar.push(`${prefijo} ${id} subgrupo ${s.etiqueta} sin par en R`); return; }
          cmp(`${prefijo} · subgrupos k (${modeloR})`, id, s.k, k[i], 0, 0);
          cmp(`${prefijo} · subgrupos efecto (${modeloR}, ${metodoR})`, `${id}/${s.etiqueta}`, t(s.estimacion), arr(sg[`est_${modeloR}`])[i], 1e-6, 1e-8);
          cmp(`${prefijo} · subgrupos IC (${modeloR}, ${metodoR})`, `${id}/${s.etiqueta}`, t(s.icInferior), arr(sg[`lb_${modeloR}`])[i], 1e-6, 1e-8);
        });
        cmp(`${prefijo} · subgrupos Q entre (${modeloR}, ${metodoR})`, id, a.qBet, sg[`qb_${modeloR}`], 1e-6, 1e-8);
        cmp(`${prefijo} · subgrupos p entre (${modeloR}, ${metodoR})`, id, a.pValor, sg[`p_${modeloR}`], 1e-5, 1e-8);
      }
    }
  }
}

// ---------------------------------------------------------------------------------------------
// 1) Efectos genéricos (pooling aislado)
// ---------------------------------------------------------------------------------------------
for (const e of ref.iv) {
  const estudios = estudiosDesdeYSe("DM", "lineal", arr(e.yi), arr(e.sei), arr(e.anio), arr(e.grupo));
  compararPool("Efectos genéricos", e.id, estudios, "DM", e.ref.pool);
  compararAdicionales("Efectos genéricos", e.id, estudios, "DM", e.ref);
}

// ---------------------------------------------------------------------------------------------
// 2) Tablas 2×2: del dato crudo al resultado
// ---------------------------------------------------------------------------------------------
const MEDIDA_BIN: Record<string, "OR" | "RR" | "RAR"> = { OR: "OR", RR: "RR", RD: "RAR" };
for (const e of ref.bin) {
  const ai = arr(e.ai), n1 = arr(e.n1), ci = arr(e.ci), n2 = arr(e.n2);
  for (const [claveR, medida] of Object.entries(MEDIDA_BIN)) {
    const R = e.ref[claveR];
    const efectos = ai.map((a, i) => efectoDesde2x2(a, n1[i], ci[i], n2[i], medida));
    const t = (v: number) => transformar(v, ESCALA[medida]);
    efectos.forEach((ef, i) => {
      if (!ef) return;
      cmp(`2×2 · efecto por estudio (${medida})`, `${e.id}#${i}`, t(ef.valor), arr(R.yi)[i], 1e-9, 1e-12);
      cmp(`2×2 · error estándar por estudio (${medida})`, `${e.id}#${i}`, ef.se, arr(R.sei)[i], 1e-9, 1e-12);
    });
    // Cadena completa: el efecto calculado por aNova (precisión completa) entra al motor
    const estCompleto = efectos.map((ef, i) =>
      estudioDesde(medida, ef!.valor, ef!.icInferior, ef!.icSuperior, {
        año: String(arr(e.anio)[i]), diseño: arr<string>(e.grupo)[i] as Estudio["diseño"],
        crudos: { tipo: "2x2", eventosTrat: String(ai[i]), totalTrat: String(n1[i]), eventosCtrl: String(ci[i]), totalCtrl: String(n2[i]) },
      }),
    );
    compararPool(`2×2 → ${medida}`, e.id, estCompleto, medida, R.pool);
    // Lo que realmente guarda la pantalla del Paso 5: el efecto con 9 cifras significativas (`formatoParaGuardar`, el mismo que usa ExtractionStep).
    // Hasta oct 2026 se guardaba con 4 decimales y el resultado combinado se apartaba de R hasta ~1,6 %.
    const fmtGuardar = (x: number) => Number(formatoParaGuardar(x));
    // La app guarda también los datos crudos 2×2 (con ellos reconoce las tablas sin eventos en ningún brazo y las excluye del OR/RR).
    const estUI = efectos.map((ef, i) =>
      estudioDesde(medida, fmtGuardar(ef!.valor), fmtGuardar(ef!.icInferior), fmtGuardar(ef!.icSuperior), {
        crudos: { tipo: "2x2", eventosTrat: String(ai[i]), totalTrat: String(n1[i]), eventosCtrl: String(ci[i]), totalCtrl: String(n2[i]) },
      }),
    );
    const aUI = calcularMetaAnalisis(estUI, medida, "DL");
    if (aUI) {
      cmp(`2×2 · tal como lo guarda la app · efecto aleatorio DL (${medida})`, e.id, t(aUI.aleatorio.estimacion), R.pool.DL.est, 1e-6, 1e-7);
      cmp(`2×2 · tal como lo guarda la app · IC aleatorio DL (${medida})`, e.id, t(aUI.aleatorio.icInferior), R.pool.DL.lb, 1e-6, 1e-7);
    }
    if (R.mh) {
      const mh = metaAnalisisMH(estCompleto, medida);
      if (mh) {
        cmp(`Mantel-Haenszel · efecto (${medida})`, e.id, t(mh.estimacion), R.mh.est, 1e-6, 1e-9);
        // Diferencia de riesgos (RAR): metafor usa otro estimador de varianza (Sato-Greenland-Robins 1989) y deja las tablas
        // doble-cero fuera del Q. aNova sigue a RevMan 5 / Cochrane Handbook (Greenland-Robins 1985, doble-cero con +0,5),
        // así que su SE, IC y Q se comparan con las fórmulas literales de RevMan transcritas en R (`mh_revman`); la
        // distancia a metafor se anota aparte como convención documentada, no como error.
        const RM = medida === "RAR" ? R.mh_revman : R.mh;
        cmp(`Mantel-Haenszel · error estándar (${medida})`, e.id, mh.se, RM.se, 1e-6, 1e-9);
        cmp(`Mantel-Haenszel · IC (${medida})`, e.id, t(mh.icInferior), RM.lb, 1e-6, 1e-9);
        cmp(`Mantel-Haenszel · IC (${medida})`, e.id, t(mh.icSuperior), RM.ub, 1e-6, 1e-9);
        cmp(`Mantel-Haenszel · estudios (incluidos + descartados) (${medida})`, e.id, mh.k + mh.excluidos, R.mh.k, 0, 0);
        cmp(`Mantel-Haenszel · Q (${medida})`, e.id, mh.heterogeneidad.q, RM.q, 1e-5, 1e-8);
        if (medida === "RAR") {
          convencionRAR.n++;
          convencionRAR.seRel = Math.max(convencionRAR.seRel, Math.abs(mh.se - R.mh.se) / R.mh.se);
          if (Math.abs(mh.heterogeneidad.q - R.mh.q) > 1e-5 * Math.max(1, R.mh.q)) convencionRAR.qDistintas++;
        }
      } else if (R.mh.k >= 2) sinPar.push(`MH ${medida} ${e.id}: aNova=null R k=${R.mh.k}`);
    }
    if (medida === "OR") {
      compararAdicionales("2×2 → OR", e.id, estCompleto, "OR", R);
      if (e.ref.peto) {
        const RP = e.ref.peto;
        const pe = metaAnalisisPeto(estCompleto);
        if (pe) {
          cmp("Peto · efecto", e.id, t(pe.estimacion), RP.est, 1e-6, 1e-9);
          cmp("Peto · error estándar", e.id, pe.se, RP.se, 1e-6, 1e-9);
          cmp("Peto · IC", e.id, t(pe.icInferior), RP.lb, 1e-6, 1e-9);
          cmp("Peto · IC", e.id, t(pe.icSuperior), RP.ub, 1e-6, 1e-9);
          cmp("Peto · estudios (incluidos + descartados)", e.id, pe.k + pe.excluidos, RP.k, 0, 0);
          cmp("Peto · Q", e.id, pe.heterogeneidad.q, RP.q, 1e-5, 1e-8);
        } else if (RP.k >= 2) sinPar.push(`Peto ${e.id}: aNova=null R k=${RP.k}`);
      }
      // Peters: como `meta::metabias`, el +0,5 corrige solo el ln(OR) y el peso usa los conteos originales (incluye conjuntos con ceros).
      if (e.ref.peters && e.ref.peters.p !== null && e.ref.peters.p !== undefined) {
        const RS = e.ref.peters;
        const pt = testPeters(estCompleto, "OR");
        if (pt) {
          cmp("Peters · pendiente", e.id, pt.pendiente, RS.pendiente, 1e-5, 1e-8);
          cmp("Peters · t", e.id, pt.t, RS.t, 1e-5, 1e-8);
          cmp("Peters · p", e.id, pt.pValor, RS.p, 1e-4, 1e-8);
        } else sinPar.push(`Peters ${e.id}: aNova=null R=resultado`);
      }
    }
  }
}

// ---------------------------------------------------------------------------------------------
// 3) Medias: DM y g de Hedges
// ---------------------------------------------------------------------------------------------
for (const e of ref.cont) {
  const m1 = arr(e.m1), sd1 = arr(e.sd1), n1 = arr(e.n1), m2 = arr(e.m2), sd2 = arr(e.sd2), n2 = arr(e.n2);
  for (const [claveR, medida, fn] of [["MD", "DM", efectoDesdeMedias], ["SMD", "DME", efectoDesdeMediasEstandarizada]] as const) {
    const R = e.ref[claveR];
    const efectos = m1.map((_, i) => fn(m1[i], sd1[i], n1[i], m2[i], sd2[i], n2[i]));
    efectos.forEach((ef, i) => {
      if (!ef) return;
      cmp(`Medias · efecto por estudio (${medida})`, `${e.id}#${i}`, ef.valor, arr(R.yi)[i], 1e-9, 1e-12);
      cmp(`Medias · error estándar por estudio (${medida})`, `${e.id}#${i}`, ef.se, arr(R.sei)[i], 1e-9, 1e-12);
    });
    const est = efectos.map((ef) => estudioDesde(medida, ef!.valor, ef!.icInferior, ef!.icSuperior));
    compararPool(`Medias → ${medida}`, e.id, est, medida, R.pool);
  }
}

// ---------------------------------------------------------------------------------------------
// 4) Proporciones (logit)
// ---------------------------------------------------------------------------------------------
for (const e of ref.prop) {
  const xi = arr(e.xi), ni = arr(e.ni);
  const R = e.ref.PLO;
  const efectos = xi.map((x, i) => efectoDesdeProporcion(x, ni[i]));
  efectos.forEach((ef, i) => {
    if (!ef) return;
    cmp("Proporciones · logit por estudio", `${e.id}#${i}`, transformar(ef.valor, "logit"), arr(R.yi)[i], 1e-9, 1e-12);
    cmp("Proporciones · error estándar por estudio", `${e.id}#${i}`, ef.se, arr(R.sei)[i], 1e-9, 1e-12);
  });
  const est = efectos.map((ef) => estudioDesde("Proporción", ef!.valor, ef!.icInferior, ef!.icSuperior));
  compararPool("Proporciones → logit", e.id, est, "Proporción", R.pool);
}

// ---------------------------------------------------------------------------------------------
// 4b) Medidas añadidas en octubre 2026: tasas (RTI, DTI), razón de medias, correlaciones (z de Fisher) y O−E/V (HR)
// ---------------------------------------------------------------------------------------------
const rmDiferenciaMetafor = { n: 0, max: 0 };
for (const e of ref.tasas) {
  const x1 = arr(e.x1), t1 = arr(e.t1), x2 = arr(e.x2), t2 = arr(e.t2), anio = arr(e.anio), grupo = arr(e.grupo) as unknown as string[];
  for (const [claveR, medida] of [["IRR", "RTI"], ["IRD", "DTI"]] as const) {
    const R = e.ref[claveR];
    const efectos = x1.map((_, i) => efectoDesdeTasas(x1[i], t1[i], x2[i], t2[i], medida));
    efectos.forEach((ef, i) => {
      if (!ef) return;
      cmp(`Tasas · efecto por estudio (${medida})`, `${e.id}#${i}`, transformar(ef.valor, ESCALA[medida]), arr(R.yi)[i], 1e-9, 1e-12);
      cmp(`Tasas · error estándar por estudio (${medida})`, `${e.id}#${i}`, ef.se, arr(R.sei)[i], 1e-9, 1e-12);
    });
    const f: number = R.factor ?? 1; // la diferencia de tasas se combina «por 1000 personas-tiempo» (ver generar-referencia.R)
    const est = efectos.map((ef, i) =>
      estudioDesde(medida, ef!.valor * f, ef!.icInferior * f, ef!.icSuperior * f, {
        año: String(anio[i]), diseño: grupo[i] as Estudio["diseño"],
        crudos: { tipo: "tasas", eventosTrat: String(x1[i]), tiempoTrat: String(t1[i]), eventosCtrl: String(x2[i]), tiempoCtrl: String(t2[i]) },
      }),
    );
    if (R.pool) {
      compararPool(`Tasas → ${medida}`, e.id, est, medida, R.pool);
      compararAdicionales(`Tasas → ${medida}`, e.id, est, medida, R);
    } else if (calcularMetaAnalisis(est, medida, "DL") !== null) {
      sinPar.push(`Tasas → ${medida} ${e.id}: aNova combina con menos de 2 estudios en R`);
    }
  }
}
for (const e of ref.rom) {
  const m1 = arr(e.m1), sd1 = arr(e.sd1), n1 = arr(e.n1), m2 = arr(e.m2), sd2 = arr(e.sd2), n2 = arr(e.n2), anio = arr(e.anio), grupo = arr(e.grupo) as unknown as string[];
  const R = e.ref.ROM;
  const efectos = m1.map((_, i) => efectoDesdeRazonDeMedias(m1[i], sd1[i], n1[i], m2[i], sd2[i], n2[i]));
  efectos.forEach((ef, i) => {
    if (!ef) return;
    cmp("Razón de medias · ln(RM) por estudio", `${e.id}#${i}`, Math.log(ef.valor), arr(R.yi)[i], 1e-9, 1e-12);
    cmp("Razón de medias · error estándar por estudio", `${e.id}#${i}`, ef.se, arr(R.sei)[i], 1e-9, 1e-12);
    // Convención documentada: metafor suma una corrección de sesgo de primer orden (Hedges 1999) al logaritmo; el Handbook no.
    rmDiferenciaMetafor.n++;
    rmDiferenciaMetafor.max = Math.max(rmDiferenciaMetafor.max, Math.abs(Math.log(ef.valor) - arr(e.yi_metafor)[i]));
  });
  const est = efectos.map((ef, i) =>
    estudioDesde("RM", ef!.valor, ef!.icInferior, ef!.icSuperior, {
      año: String(anio[i]), diseño: grupo[i] as Estudio["diseño"],
      crudos: { tipo: "medias", mediaTrat: String(m1[i]), deTrat: String(sd1[i]), nTrat: String(n1[i]), mediaCtrl: String(m2[i]), deCtrl: String(sd2[i]), nCtrl: String(n2[i]) },
    }),
  );
  compararPool("Razón de medias → RM", e.id, est, "RM", R.pool);
  compararAdicionales("Razón de medias → RM", e.id, est, "RM", R);
}
for (const e of ref.cor) {
  const ri = arr(e.ri), ni = arr(e.ni), anio = arr(e.anio), grupo = arr(e.grupo) as unknown as string[];
  const R = e.ref.ZCOR;
  const efectos = ri.map((r, i) => efectoDesdeCorrelacion(r, ni[i]));
  efectos.forEach((ef, i) => {
    if (!ef) return;
    cmp("Correlaciones · z de Fisher por estudio", `${e.id}#${i}`, transformar(ef.valor, "fisher"), arr(R.yi)[i], 1e-9, 1e-12);
    cmp("Correlaciones · error estándar por estudio", `${e.id}#${i}`, ef.se, arr(R.sei)[i], 1e-9, 1e-12);
  });
  const est = efectos.map((ef, i) =>
    estudioDesde("r", ef!.valor, ef!.icInferior, ef!.icSuperior, {
      año: String(anio[i]), diseño: grupo[i] as Estudio["diseño"], crudos: { tipo: "correlacion", r: String(ri[i]), n: String(ni[i]) },
    }),
  );
  compararPool("Correlaciones → z de Fisher", e.id, est, "r", R.pool);
  compararAdicionales("Correlaciones → z de Fisher", e.id, est, "r", R);
}
for (const e of ref.oev) {
  const oe = arr(e.oe), v = arr(e.v), anio = arr(e.anio), grupo = arr(e.grupo) as unknown as string[];
  const R = e.ref.HR;
  const efectos = oe.map((o, i) => efectoDesdeOEV(o, v[i]));
  efectos.forEach((ef, i) => {
    if (!ef) return;
    cmp("O−E y V · ln(HR) por estudio", `${e.id}#${i}`, Math.log(ef.valor), arr(R.yi)[i], 1e-9, 1e-12);
    cmp("O−E y V · error estándar por estudio", `${e.id}#${i}`, ef.se, arr(R.sei)[i], 1e-9, 1e-12);
  });
  const est = efectos.map((ef, i) =>
    estudioDesde("HR", ef!.valor, ef!.icInferior, ef!.icSuperior, {
      año: String(anio[i]), diseño: grupo[i] as Estudio["diseño"], crudos: { tipo: "oev", oe: String(oe[i]), v: String(v[i]) },
    }),
  );
  compararPool("O−E y V → HR", e.id, est, "HR", R.pool);
  compararAdicionales("O−E y V → HR", e.id, est, "HR", R);
  // El efecto fijo por varianza inversa de los ln(HR) es exactamente el método de Peto: Σ(O−E)/ΣV, con EE = 1/√ΣV.
  const fijo = calcularMetaAnalisis(est, "HR", "DL");
  if (fijo) {
    cmp("O−E y V · efecto fijo = Σ(O−E)/ΣV (Peto)", e.id, Math.log(fijo.fijo.estimacion), e.peto.est, 1e-9, 1e-12);
    cmp("O−E y V · EE fijo = 1/√ΣV (Peto)", e.id, (Math.log(fijo.fijo.icSuperior) - Math.log(fijo.fijo.icInferior)) / (2 * Z95), e.peto.se, 1e-8, 1e-12);
  }
}

// ---------------------------------------------------------------------------------------------
// 5) Meta-regresión, SROC y red (en comparar-avanzada.mts)
// ---------------------------------------------------------------------------------------------
compararAvanzada({ cmp, sinPar, arr, estudioDesde }, refAv);
const infoBivariado = compararBivariado({ cmp, sinPar, arr, estudioDesde });
const infoMultinivel = compararMultinivel({ cmp, sinPar, arr, estudioDesde });
const infoGLMM = compararGLMM({ cmp, sinPar, arr, estudioDesde });

// ---------------------------------------------------------------------------------------------
// Informe
// ---------------------------------------------------------------------------------------------
const lista = [...grupos.entries()].sort((a, b) => a[0].localeCompare(b[0], "es"));
let totalComparaciones = 0;
let gruposFallidos = 0;
const filas: string[] = [];
for (const [nombre, g] of lista) {
  totalComparaciones += g.n;
  const ok = g.fallos.length === 0;
  if (!ok) gruposFallidos++;
  filas.push(`${ok ? "✓" : "✗"} ${nombre.padEnd(66)} n=${String(g.n).padStart(6)}  dif. máx. abs ${g.maxAbs.toExponential(1).padStart(8)}  rel ${g.maxRel.toExponential(1).padStart(8)}${ok ? "" : `  ← ${g.fallos.length} fuera de tolerancia`}`);
}
console.log(`\nR: ${ref.meta.R} · metafor ${ref.meta.paquetes.metafor} · meta ${ref.meta.paquetes.meta} · semilla ${ref.meta.semilla}`);
console.log(`Escenarios: ${ref.meta.escenarios.iv} de efectos genéricos, ${ref.meta.escenarios.bin} de tablas 2×2, ${ref.meta.escenarios.cont} de medias, ${ref.meta.escenarios.prop} de proporciones`);
console.log(`Avanzadas (R ${refAv.meta.paquetes.netmeta ? "· netmeta " + refAv.meta.paquetes.netmeta : ""}): ${Object.entries(refAv.meta.escenarios).map(([k, v]) => `${v} de ${({ reg: "meta-regresión", sroc: "SROC", red: "red" } as Record<string, string>)[k] ?? k}`).join(", ")}\n`);
console.log(`Bivariado diagnóstico (R · mada ${infoBivariado.mada} · lme4 ${infoBivariado.lme4}): ${infoBivariado.conjuntos} conjuntos\n`);
console.log(`Multinivel con varianza robusta (R · metafor ${infoMultinivel.metafor} · clubSandwich ${infoMultinivel.clubSandwich}): ${infoMultinivel.conjuntos} conjuntos\n`);
console.log(`Modelos mixtos logísticos (R · metafor ${infoGLMM.metafor} · BiasedUrn ${infoGLMM.BiasedUrn}): ${infoGLMM.conjuntos} conjuntos\n`);
filas.forEach((f) => console.log(f));
if (rmDiferenciaMetafor.n) {
  console.log(
    `\nRazón de medias: aNova usa la fórmula del Handbook (Friedrich 2008, sin corrección de sesgo). metafor::escalc("ROM") añade una corrección de primer orden (Hedges 1999): diferencia máxima del logaritmo en estos ${rmDiferenciaMetafor.n} estudios = ${rmDiferenciaMetafor.max.toExponential(1)} (convención documentada, no error).`,
  );
}
if (convencionRAR.n) console.log(`\nMH diferencia de riesgos: aNova sigue a RevMan 5 (Greenland-Robins 1985; doble-cero con +0,5 en el Q). Frente a metafor (Sato-Greenland-Robins 1989; Q sin doble-cero) el SE se aparta hasta ${(convencionRAR.seRel * 100).toFixed(1)} % y el Q difiere en ${convencionRAR.qDistintas} de ${convencionRAR.n} conjuntos: convención documentada, no error.`);
console.log(`\n${totalComparaciones} comparaciones numéricas en ${lista.length} grupos; ${gruposFallidos} grupo(s) con diferencias fuera de tolerancia.`);
if (sinPar.length) {
  console.log(`\nCasos donde solo uno de los dos programas da resultado (${sinPar.length}):`);
  sinPar.slice(0, DETALLE ? 500 : 12).forEach((s) => console.log("  - " + s));
  if (!DETALLE && sinPar.length > 12) console.log(`  … y ${sinPar.length - 12} más (--detalle para verlos todos)`);
}
for (const [nombre, g] of lista) {
  if (g.fallos.length === 0) continue;
  console.log(`\n✗ ${nombre} — peor caso: ${g.peor}`);
  g.fallos.slice(0, DETALLE ? 500 : 4).forEach((f) => console.log("    " + f));
  if (!DETALLE && g.fallos.length > 4) console.log(`    … y ${g.fallos.length - 4} más`);
}
process.exit(gruposFallidos === 0 && sinPar.length === 0 ? 0 : 1);
