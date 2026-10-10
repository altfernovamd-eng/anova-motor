// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Parte «avanzada» de la validación diferencial contra R: meta-regresión, SROC (precisión diagnóstica) y meta-análisis en red.
 * La invoca `comparar.mts` con sus ayudantes (`cmp`, `sinPar`, …) para que todo salga en un solo informe.
 * Referencia: `referencias/referencia-avanzada.json`, generada por `generar-referencia-avanzada.R`.
 */
import {
  calcularSROC,
  metaRegresionMulti,
  MODERADORES_DISPONIBLES,
  type MetodoTau2,
} from "../src/statistics.ts";
import { calcularNMA, calcularNodeSplittingSIDDE, calcularPScores } from "../src/red.ts";
import type { BrazoRed, Estudio, EstudioRed, Medida } from "../src/tipos.ts";

const Z95 = 1.959963985;
// eslint-disable-next-line @typescript-eslint/no-explicit-any
type J = any;

export interface AyudasComparar {
  cmp: (grupo: string, id: string, a: number | null | undefined, r: number | null | undefined, tolRel?: number, tolAbs?: number) => void;
  sinPar: string[];
  arr: (v: J) => number[];
  estudioDesde: (medida: Medida, valor: number, icInf: number, icSup: number, extra?: Partial<Estudio>) => Estudio;
}

const expit = (x: number) => 1 / (1 + Math.exp(-x));

// ---------------------------------------------------------------------------------------------
// Meta-regresión (REML, pruebas z) — metafor::rma(mods = …)
// ---------------------------------------------------------------------------------------------
function compararMetaRegresion(h: AyudasComparar, escenarios: J[]) {
  const { cmp, sinPar, arr, estudioDesde } = h;
  for (const e of escenarios) {
    const yi = arr(e.yi), sei = arr(e.sei), anio = arr(e.anio), n = arr(e.n), grupo = ([] as string[]).concat(e.grupo);
    const estudios = yi.map((y, i) =>
      estudioDesde("DM", y, y - Z95 * sei[i], y + Z95 * sei[i], { año: String(anio[i]), n: String(n[i]), diseño: grupo[i] as Estudio["diseño"] }),
    );
    const ID_ANOVA: Record<string, string> = { anio: "año", n: "n", diseno: "diseño" };
    const mods = ([] as string[]).concat(e.moderadores).map((id) => MODERADORES_DISPONIBLES.find((m) => m.id === ID_ANOVA[id])!);
    const r = metaRegresionMulti(estudios, "DM", mods);
    if (!r) {
      sinPar.push(`Meta-regresión ${e.id} (${e.tipo}): aNova=null R=resultado`);
      continue;
    }
    const R = e.ref;
    const beta = arr(R.beta);
    if (r.coeficientes.length !== beta.length) {
      sinPar.push(`Meta-regresión ${e.id}: aNova tiene ${r.coeficientes.length} coeficientes, R ${beta.length}`);
      continue;
    }
    const G = "Meta-regresión (REML)";
    r.coeficientes.forEach((c, j) => {
      cmp(`${G} · coeficiente`, `${e.id}#${j}`, c.estimacion, beta[j], 1e-6, 1e-8);
      cmp(`${G} · error estándar`, `${e.id}#${j}`, c.se, arr(R.se)[j], 1e-6, 1e-9);
      cmp(`${G} · z`, `${e.id}#${j}`, c.z, arr(R.z)[j], 1e-6, 1e-7);
      cmp(`${G} · p del coeficiente`, `${e.id}#${j}`, c.pValor, arr(R.p)[j], 1e-6, 1e-10);
      cmp(`${G} · IC del coeficiente`, `${e.id}#${j}`, c.icInferior, arr(R.lb)[j], 1e-6, 1e-8);
      cmp(`${G} · IC del coeficiente`, `${e.id}#${j}`, c.icSuperior, arr(R.ub)[j], 1e-6, 1e-8);
    });
    cmp(`${G} · τ² residual`, e.id, r.tau2Residual, R.tau2 ?? 0, 1e-6, 1e-9);
    cmp(`${G} · QM (test ómnibus)`, e.id, r.qm, R.qm, 1e-6, 1e-7);
    cmp(`${G} · grados de libertad de QM`, e.id, r.qmGl, R.qmGl, 0, 0);
    cmp(`${G} · p de QM`, e.id, r.qmPValor, R.qmP, 1e-6, 1e-10);
    // R² = (τ²₀ − τ²)/τ²₀: con τ²₀ ≈ 0 en R (un cero numérico, p. ej. 1e-15) el cociente no tiene significado.
    if ((R.tau2_nulo ?? 0) > 1e-8) cmp(`${G} · R² (% de heterogeneidad explicada)`, e.id, r.r2, R.r2 ?? 0, 1e-5, 1e-5);
    else cmp(`${G} · R² con τ² nulo ≈ 0 (debe ser 0)`, e.id, r.r2, 0, 0, 1e-6);
  }
}

// ---------------------------------------------------------------------------------------------
// SROC / precisión diagnóstica
// ---------------------------------------------------------------------------------------------
function compararSROC(h: AyudasComparar, escenarios: J[]) {
  const { cmp, sinPar, arr } = h;
  for (const e of escenarios) {
    const vp = arr(e.vp), fp = arr(e.fp), fn = arr(e.fn), vn = arr(e.vn);
    const estudios = vp.map((_, i) => ({
      id: i + 1, autor: `D${i + 1}`, año: "", diseño: "Transversal", n: "", medida: "OR", valor: "", icInferior: "", icSuperior: "",
      datosDiagnosticos: { vp: String(vp[i]), fp: String(fp[i]), fn: String(fn[i]), vn: String(vn[i]) },
    })) as Estudio[];
    const R = e.ref;

    for (const [etiqueta, modelo, metodo] of [["fijo", "fijo", "DL"], ["DL", "aleatorio", "DL"], ["PM", "aleatorio", "PM"], ["REML", "aleatorio", "REML"]] as const) {
      const s = calcularSROC(estudios, modelo, metodo as MetodoTau2, { bivariado: false }); // el modelo bivariado se valida aparte (comparar-bivariado.mts)
      if (!s) {
        sinPar.push(`SROC ${e.id} ${etiqueta}: aNova=null R=resultado`);
        continue;
      }
      const G = `SROC · ${etiqueta === "fijo" ? "efectos fijos" : `efectos aleatorios ${etiqueta}`}`;
      const pool = (clave: "dor" | "sens" | "esp") => (etiqueta === "fijo" ? R[clave].fijo : R[clave][etiqueta]);
      const pd = pool("dor"), ps = pool("sens"), pe = pool("esp");
      if (!pd || !ps || !pe) continue;
      cmp(`${G} · DOR combinado`, e.id, s.dor.estimacion, Math.exp(pd.est), 1e-6, 1e-9);
      cmp(`${G} · IC del DOR`, e.id, s.dor.icInferior, Math.exp(pd.lb), 1e-6, 1e-9);
      cmp(`${G} · IC del DOR`, e.id, s.dor.icSuperior, Math.exp(pd.ub), 1e-6, 1e-9);
      cmp(`${G} · sensibilidad combinada`, e.id, s.sensibilidadPromedio.estimacion, expit(ps.est), 1e-6, 1e-9);
      cmp(`${G} · IC de la sensibilidad`, e.id, s.sensibilidadPromedio.icInferior, expit(ps.lb), 1e-6, 1e-9);
      cmp(`${G} · IC de la sensibilidad`, e.id, s.sensibilidadPromedio.icSuperior, expit(ps.ub), 1e-6, 1e-9);
      cmp(`${G} · especificidad combinada`, e.id, s.especificidadPromedio.estimacion, expit(pe.est), 1e-6, 1e-9);
      cmp(`${G} · IC de la especificidad`, e.id, s.especificidadPromedio.icInferior, expit(pe.lb), 1e-6, 1e-9);
      cmp(`${G} · IC de la especificidad`, e.id, s.especificidadPromedio.icSuperior, expit(pe.ub), 1e-6, 1e-9);
      if (etiqueta !== "fijo") cmp(`${G} · τ² del DOR`, e.id, s.dor.tau2, pd.tau2, 1e-6, 1e-9);
      // mada::madauni como caja negra (solo en la referencia masiva): DOR, IC, τ² y Q de DerSimonian-Laird con la corrección «single» (la de aNova)
      if (etiqueta === "DL" && e.mada) {
        const M = e.mada;
        cmp("SROC · mada::madauni (DSL) · DOR combinado", e.id, s.dor.estimacion, Math.exp(M.lnDOR), 1e-6, 1e-9);
        cmp("SROC · mada::madauni (DSL) · IC del DOR", e.id, s.dor.icInferior, Math.exp(M.lnDOR - Z95 * M.se), 1e-6, 1e-9);
        cmp("SROC · mada::madauni (DSL) · IC del DOR", e.id, s.dor.icSuperior, Math.exp(M.lnDOR + Z95 * M.se), 1e-6, 1e-9);
        cmp("SROC · mada::madauni (DSL) · τ² del DOR", e.id, s.dor.tau2, M.tau2, 1e-6, 1e-9);
        // La «Q de Cochran» que devuelve `madauni` con DSL usa los pesos de EFECTOS ALEATORIOS (1/(v + τ²)); la de aNova es la de Cochran con pesos
        // de efectos fijos (la que da el I²), igual que `metafor` y que la comprobación a mano. No es comparable: se verifica contra metafor/meta arriba.
        cmp("SROC · mada::madauni (DSL) · gl del DOR", e.id, s.dor.gl, M.gl, 1e-12, 1e-12);
      }
      if (etiqueta === "fijo") {
        cmp("SROC · Q del DOR", e.id, s.dor.q, R.dor.q, 1e-6, 1e-8);
        cmp("SROC · p de Q del DOR", e.id, s.dor.pValor, R.dor.pq, 1e-5, 1e-6);
        cmp("SROC · I² del DOR", e.id, s.dor.i2, R.dor.i2, 1e-6, 1e-6);
        // Regresión de Moses-Littenberg (≥4 estudios) y Wilson por estudio: no dependen del modelo; se comparan una vez.
        if (R.regresion) {
          if (!s.regresion) sinPar.push(`SROC ${e.id}: regresión de Moses-Littenberg aNova=null R=resultado`);
          else {
            cmp("SROC · Moses-Littenberg · intercepto a", e.id, s.regresion.a, R.regresion.a, 1e-6, 1e-8);
            cmp("SROC · Moses-Littenberg · pendiente b", e.id, s.regresion.b, R.regresion.b, 1e-6, 1e-8);
            cmp("SROC · Moses-Littenberg · EE de a", e.id, s.regresion.seA, R.regresion.seA, 1e-6, 1e-9);
            cmp("SROC · Moses-Littenberg · EE de b", e.id, s.regresion.seB, R.regresion.seB, 1e-6, 1e-9);
            cmp("SROC · Moses-Littenberg · z de la pendiente", e.id, s.regresion.zB, R.regresion.zB, 1e-6, 1e-7);
            cmp("SROC · Moses-Littenberg · p de la pendiente", e.id, s.regresion.pValorB, R.regresion.pB, 1e-6, 1e-10);
          }
        } else if (s.regresion) sinPar.push(`SROC ${e.id}: regresión de Moses-Littenberg aNova=resultado R=sin resultado`);
        s.puntos.forEach((p, i) => {
          cmp("SROC · IC de Wilson de la sensibilidad (por estudio)", `${e.id}#${i}`, p.sensIcInferior, arr(R.wilson_sens.lb)[i], 1e-6, 1e-8);
          cmp("SROC · IC de Wilson de la sensibilidad (por estudio)", `${e.id}#${i}`, p.sensIcSuperior, arr(R.wilson_sens.ub)[i], 1e-6, 1e-8);
          cmp("SROC · IC de Wilson de la especificidad (por estudio)", `${e.id}#${i}`, p.espIcInferior, arr(R.wilson_esp.lb)[i], 1e-6, 1e-8);
          cmp("SROC · IC de Wilson de la especificidad (por estudio)", `${e.id}#${i}`, p.espIcSuperior, arr(R.wilson_esp.ub)[i], 1e-6, 1e-8);
        });
      }
    }
  }
}

// ---------------------------------------------------------------------------------------------
// Meta-análisis en red — netmeta (Rücker)
// ---------------------------------------------------------------------------------------------
const lista = <T,>(v: T | T[] | null | undefined): T[] => (v === null || v === undefined ? [] : ([] as T[]).concat(v));

/** Estudios de la red de un escenario: una fila por brazo en la referencia, agrupadas por estudio. */
function estudiosDeRed(e: J): EstudioRed[] {
  const medida = e.medida as "DM" | "OR";
  const porEstudio = new Map<string, BrazoRed[]>();
  for (const f of lista<J>(e.datos)) {
    const brazos = porEstudio.get(f.study) ?? [];
    brazos.push(medida === "DM"
      ? { tratamiento: f.treat, media: String(f.mean), de: String(f.sd), n: String(f.n) }
      : { tratamiento: f.treat, eventos: String(f.event), total: String(f.n) });
    porEstudio.set(f.study, brazos);
  }
  return [...porEstudio.entries()].map(([etiqueta, brazos], i) => ({ id: i + 1, etiqueta, medida, brazos }));
}

/** Brazos de OR por estudio de un escenario (para contar los estudios de dos brazos con doble cero). */
function porEstudioOR(e: J): Map<string, { tratamiento: string; eventos: number; total: number }[]> {
  const m = new Map<string, { tratamiento: string; eventos: number; total: number }[]>();
  for (const f of lista<J>(e.datos)) m.set(f.study, [...(m.get(f.study) ?? []), { tratamiento: f.treat, eventos: f.event, total: f.n }]);
  return m;
}

function compararRed(h: AyudasComparar, escenarios: J[]) {
  const { cmp, sinPar } = h;
  for (const e of escenarios) {
    const medida = e.medida as "DM" | "OR";
    const trts: string[] = lista(e.tratamientos);
    const estudiosRed = estudiosDeRed(e);
    const R = e.ref;
    const t = (x: number) => (medida === "OR" ? Math.log(x) : x);

    for (const [clave, modelo] of [["common", "fijo"], ["random", "aleatorio"]] as const) {
      const G = `Red (${medida === "DM" ? "diferencia de medias" : "OR"}) · ${modelo === "fijo" ? "efecto común" : "efectos aleatorios"}`;
      const RM = R[clave];
      const r = calcularNMA(estudiosRed, medida, modelo);
      if (!r) {
        sinPar.push(`Red ${e.id} ${modelo}: aNova=null R=resultado`);
        continue;
      }
      // `netmeta` descarta del OR los estudios de dos brazos con doble cero (sin eventos, o con evento en todos, en los dos brazos);
      // aNova descarta además los multibrazo cuyos brazos están todos en cero (ningún contraste informa; a R se le entregan sin ellos): k es el de la red resultante
      const doblesCeros = [...porEstudioOR(e)].filter(([, brazos]) => brazos.every((b) => b.eventos === 0) || brazos.every((b) => b.eventos === b.total)).length;
      cmp(`${G} · estudios (k)`, e.id, r.k, e.k - (medida === "OR" ? doblesCeros : 0), 0, 0);
      cmp(`${G} · comparaciones (m)`, e.id, r.m, R.m, 0, 0);
      cmp(`${G} · tratamientos (n)`, e.id, r.n, R.n, 0, 0);
      cmp(`${G} · Q total`, e.id, r.q, R.q, 1e-6, 1e-8);
      cmp(`${G} · gl de Q`, e.id, r.gl, R.gl, 0, 0);
      cmp(`${G} · p de Q`, e.id, r.pValor, R.pq, 1e-6, 1e-10);
      cmp(`${G} · I²`, e.id, r.i2, R.i2, 1e-6, 1e-6);
      if (modelo === "aleatorio") cmp(`${G} · τ²`, e.id, r.tau2, R.tau2, 1e-6, 1e-9);

      // Efecto de red para cada par de tratamientos (orientación tratamiento1 − tratamiento2)
      const pares: J[] = lista(RM.pares);
      for (const est of r.estimaciones) {
        const directo = pares.find((p) => p.t1 === est.tratamiento1 && p.t2 === est.tratamiento2);
        const inverso = pares.find((p) => p.t1 === est.tratamiento2 && p.t2 === est.tratamiento1);
        const par = directo ?? inverso;
        if (!par) { sinPar.push(`Red ${e.id}: par ${est.tratamiento1}–${est.tratamiento2} sin referencia en R`); continue; }
        const signo = directo ? 1 : -1;
        const id = `${e.id} ${est.tratamiento1}-${est.tratamiento2}`;
        cmp(`${G} · efecto de red`, id, est.diff, signo * par.te, 1e-6, 1e-9);
        cmp(`${G} · error estándar de red`, id, est.se, par.se, 1e-6, 1e-9);
        const lo = directo ? par.lb : -par.ub, hi = directo ? par.ub : -par.lb;
        cmp(`${G} · IC de red`, id, t(est.icInferior), lo, 1e-6, 1e-9);
        cmp(`${G} · IC de red`, id, t(est.icSuperior), hi, 1e-6, 1e-9);
      }

      // P-scores (menor valor = mejor, como small.values = "good")
      const ps = calcularPScores(r, true);
      trts.forEach((tr, i) => {
        const x = ps.find((p) => p.tratamiento === tr)?.pScore;
        cmp(`${G} · P-score`, `${e.id} ${tr}`, x, lista<number>(RM.pscore)[i], 1e-6, 1e-8);
      });

      // Descomposición de Q en dentro de diseños / entre diseños (netmeta::decomp.design, modelo común)
      if (modelo === "fijo") {
        const dq = r.decomposicionQ;
        const Rd = R.decomp;
        const dentroR = Rd ? lista<number>(Rd.dentro) : null;
        if (dq && Rd) {
          const entreR = lista<number>(Rd.entre);
          cmp("Red · descomposición de Q · Q dentro de diseños", e.id, dq.qDentro, dentroR![0], 1e-6, 1e-8);
          cmp("Red · descomposición de Q · gl dentro", e.id, dq.glDentro, dentroR![1], 0, 0);
          if (dq.glDentro > 0) cmp("Red · descomposición de Q · p dentro", e.id, dq.pDentro, dentroR![2], 1e-5, 1e-9);
          cmp("Red · descomposición de Q · Q entre diseños", e.id, dq.qEntre, entreR[0], 1e-6, 1e-8);
          cmp("Red · descomposición de Q · gl entre", e.id, dq.glEntre, entreR[1], 0, 0);
          if (dq.glEntre > 0) cmp("Red · descomposición de Q · p entre", e.id, dq.pEntre, entreR[2], 1e-5, 1e-9);
        } else if (!!dq !== (!!dentroR && dentroR[1] > 0)) {
          sinPar.push(`Red ${e.id}: descomposición de Q aNova=${dq ? "resultado" : "null"} R df dentro=${dentroR ? dentroR[1] : "—"}`);
        }
      }

      // Node-splitting: back-calculation (en r.nodeSplitting) y SIDDE
      const comparar = (etiqueta: string, nsA: typeof r.nodeSplitting, nsR: J[], prop: boolean) => {
        const conSplit = nsR.filter((x) => x.dir_te !== null && x.ind_te !== null && x.dir_te !== undefined && x.ind_te !== undefined);
        for (const a of nsA) {
          const dir = conSplit.find((x) => x.t1 === a.tratamiento1 && x.t2 === a.tratamiento2);
          const inv = conSplit.find((x) => x.t1 === a.tratamiento2 && x.t2 === a.tratamiento1);
          const x = dir ?? inv;
          const id = `${e.id} ${a.tratamiento1}-${a.tratamiento2}`;
          if (!x) { sinPar.push(`Red ${etiqueta} ${modelo} ${id}: aNova tiene la comparación y R no`); continue; }
          const sg = dir ? 1 : -1;
          cmp(`${G} · ${etiqueta} · estimación directa`, id, a.diffDirecto, sg * x.dir_te, 1e-6, 1e-9);
          cmp(`${G} · ${etiqueta} · estimación indirecta`, id, t(a.estimacionIndirecta), sg * x.ind_te, 1e-6, 1e-8);
          cmp(`${G} · ${etiqueta} · |z| directa vs indirecta`, id, Math.abs(a.z), Math.abs(x.z), 1e-5, 1e-7);
          cmp(`${G} · ${etiqueta} · p de la diferencia`, id, a.pValor, x.p, 1e-6, 1e-10);
          if (prop) cmp(`${G} · ${etiqueta} · proporción de evidencia directa`, id, a.proporcionDirecta, x.prop, 1e-6, 1e-9);
        }
        for (const x of conSplit) {
          const hay = nsA.some((a) => (a.tratamiento1 === x.t1 && a.tratamiento2 === x.t2) || (a.tratamiento1 === x.t2 && a.tratamiento2 === x.t1));
          if (!hay) sinPar.push(`Red ${etiqueta} ${modelo} ${e.id} ${x.t1}-${x.t2}: R tiene la comparación y aNova no`);
        }
      };
      comparar("node-splitting (back-calculation)", r.nodeSplitting, lista<J>(RM.split_bc), true);
      comparar("node-splitting (SIDDE)", calcularNodeSplittingSIDDE(estudiosRed, medida, modelo), lista<J>(RM.split_sidde), false);
    }
  }
}

export function compararAvanzada(h: AyudasComparar, refAv: J) {
  compararMetaRegresion(h, refAv.reg ?? []);
  compararSROC(h, refAv.sroc ?? []);
  compararRed(h, refAv.red ?? []);
}
