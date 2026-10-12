// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * **Efectos dependientes: modelo multinivel y varianza robusta** (Cochrane Handbook 6.2.4 y 6.2.9 para el problema; Pustejovsky y Tipton 2022 para el
 * modelo CHE; Hedges, Tipton y Johnson 2010, Tipton 2015 y Pustejovsky y Tipton 2018 para la varianza robusta).
 *
 * Cuando un mismo estudio aporta varios efectos (varios desenlaces, seguimientos o comparaciones), esos efectos no son independientes y
 * un meta-análisis que los trata como tales sobreestima la precisión. Aquí se modela la dependencia de dos maneras a la vez:
 *
 *  1. **Modelo CHE** (correlated and hierarchical effects): yᵢⱼ = x'ᵢⱼβ + uᵢ + wᵢⱼ + eᵢⱼ, con uᵢ ~ N(0, σ²ᵤ) (estudio), wᵢⱼ ~ N(0, σ²ᵥ)
 *     (efecto dentro del estudio) y eᵢⱼ el error de muestreo, con varianza conocida vᵢⱼ y correlación supuesta ρ entre los de un mismo
 *     estudio (Cov(eᵢⱼ, eᵢₖ) = ρ·√(vᵢⱼ vᵢₖ)). Los dos componentes de varianza se estiman por REML; β por mínimos cuadrados
 *     generalizados con ese modelo de trabajo.
 *  2. **Varianza robusta (RVE)** sobre ese mismo ajuste: el error estándar de β no depende de que el modelo de trabajo (en particular la
 *     ρ supuesta) sea correcto. Se usa el estimador sándwich CR2 de Bell y McCaffrey con la corrección de Pustejovsky y Tipton (2018)
 *     y los grados de libertad de Satterthwaite para un coeficiente (Tipton 2015), con los que se hacen el contraste t y el intervalo.
 *
 * Con una sola observación por estudio el modelo se reduce a efectos aleatorios con REML (σ²ᵤ no es identificable: se fija en 0).
 *
 * Sin dependencias de la app (solo `funciones-especiales` y `numerica-compartida`). Lógica pura.
 *
 * Fuentes de las fórmulas: Pustejovsky y Tipton 2022 (CHE); Pustejovsky y Tipton 2018 (CR2: A_j = D_j' B_j^{-1/2} D_j); Tipton 2015
 * (Satterthwaite: df = (Σⱼ bⱼ'Φbⱼ)² / Σⱼₖ (bⱼ'Φbₖ)², bajo el modelo de trabajo Φ); Nelder y Mead 1965. Detalle de procedencia y de
 * la validación contra R en docs/PROCEDENCIA.md y VALIDACION_ESTADISTICA.md.
 */
import { normalCuantilDosColas, tStudentCuantilDosColas, tStudentDosColas, normalDosColas } from "./funciones-especiales.ts";
import { cholesky, inversaSPD, matriz, mult, multVec, nelderMead, raizInversaSimetrica, trans, type Mat } from "./numerica-compartida.ts";

// ---------------------------------------------------------------------------------------------
// Tipos
// ---------------------------------------------------------------------------------------------

/** Un efecto: valor y varianza de muestreo en la escala de análisis, y a qué estudio (conglomerado) pertenece. */
export interface EfectoDependiente {
  etiqueta?: string;
  y: number;
  /** Varianza de muestreo (el cuadrado del error estándar). */
  v: number;
  estudio: string;
  /** Valores de las covariables (moderadores) de este efecto, en el mismo orden que `nombresModeradores`. */
  x?: number[];
}

export interface OpcionesMultinivel {
  /** Correlación supuesta entre los errores de muestreo de los efectos de un mismo estudio. Por defecto 0,6. */
  rho?: number;
  nombresModeradores?: string[];
}

export interface CoeficienteMultinivel {
  nombre: string;
  estimacion: number;
  /** Error estándar del modelo de trabajo (supone que el modelo es correcto). */
  seModelo: number;
  /** Error estándar robusto CR2. */
  seRobusto: number;
  /** Grados de libertad de Satterthwaite del contraste robusto. */
  gl: number;
  t: number;
  pValor: number;
  /** IC 95 % robusto (t con los grados de libertad de Satterthwaite). */
  icInferior: number;
  icSuperior: number;
  /** IC 95 % del modelo de trabajo (normal), para comparar. */
  icModeloInferior: number;
  icModeloSuperior: number;
}

export interface ResultadoMultinivel {
  /** Número de efectos y de estudios (conglomerados). */
  k: number;
  m: number;
  /** Cuántos estudios aportan más de un efecto, y el máximo por estudio. */
  estudiosConVarios: number;
  maxPorEstudio: number;
  rho: number;
  /** Varianza entre estudios (σ²ᵤ) y entre efectos de un mismo estudio (σ²ᵥ). */
  sigma2Estudio: number;
  sigma2Efecto: number;
  coeficientes: CoeficienteMultinivel[];
  /** Log-verosimilitud restringida (REML), con sus constantes. */
  logVerosimilitudREML: number;
  convergio: boolean;
  /** Componentes de varianza que quedaron en cero. */
  limites: ("estudio" | "efecto")[];
  advertencias: string[];
}

// ---------------------------------------------------------------------------------------------
// Ajuste CHE por REML
// ---------------------------------------------------------------------------------------------

interface Cluster {
  idx: number[];
  y: number[];
  X: Mat;
  /** Covarianza de muestreo con la ρ supuesta. */
  V: Mat;
}

function construirClusters(datos: EfectoDependiente[], X: Mat, rho: number): Cluster[] {
  const mapa = new Map<string, number[]>();
  datos.forEach((d, i) => {
    const l = mapa.get(d.estudio) ?? [];
    l.push(i);
    mapa.set(d.estudio, l);
  });
  return [...mapa.values()].map((idx) => {
    const n = idx.length;
    const V = matriz(n, n);
    for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) V[a][b] = a === b ? datos[idx[a]].v : rho * Math.sqrt(datos[idx[a]].v * datos[idx[b]].v);
    return { idx, y: idx.map((i) => datos[i].y), X: idx.map((i) => X[i]), V };
  });
}

/** Σⱼ = σ²ᵤ·11' + σ²ᵥ·I + Vⱼ. */
function sigmaCluster(c: Cluster, s2u: number, s2w: number): Mat {
  return c.V.map((fila, a) => fila.map((v, b) => v + s2u + (a === b ? s2w : 0)));
}

interface Gls {
  beta: number[];
  M: Mat;
  /** Σ⁻¹ por conglomerado y log-determinantes. */
  W: Mat[];
  sumaLogDet: number;
  residuo2: number;
  logDetXWX: number;
}

function gls(clusters: Cluster[], s2u: number, s2w: number, p: number): Gls | null {
  const W: Mat[] = [];
  let sumaLogDet = 0;
  let A = matriz(p, p);
  const b = new Array<number>(p).fill(0);
  for (const c of clusters) {
    const r = inversaSPD(sigmaCluster(c, s2u, s2w));
    if (!r) return null;
    W.push(r.inv);
    sumaLogDet += r.logDet;
    const XtW = mult(trans(c.X), r.inv);
    const XtWX = mult(XtW, c.X);
    A = A.map((fila, i) => fila.map((v, j) => v + XtWX[i][j]));
    const Xy = multVec(XtW, c.y);
    for (let i = 0; i < p; i++) b[i] += Xy[i];
  }
  const inv = inversaSPD(A);
  if (!inv) return null;
  const beta = multVec(inv.inv, b);
  let residuo2 = 0;
  clusters.forEach((c, j) => {
    const r = c.y.map((yy, a) => yy - c.X[a].reduce((s, x, i) => s + x * beta[i], 0));
    const Wr = multVec(W[j], r);
    residuo2 += r.reduce((s, v, a) => s + v * Wr[a], 0);
  });
  return { beta, M: inv.inv, W, sumaLogDet, residuo2, logDetXWX: inv.logDet };
}

/** Objetivo REML (a minimizar, sin constantes): ½[Σ ln|Σⱼ| + ln|X'WX| + r'Wr]. */
function objetivoREML(clusters: Cluster[], s2u: number, s2w: number, p: number): number {
  const g = gls(clusters, s2u, s2w, p);
  return g ? 0.5 * (g.sumaLogDet + g.logDetXWX + g.residuo2) : 1e12;
}

const SIGMA2_MIN = 1e-9;

// ---------------------------------------------------------------------------------------------
// Varianza robusta CR2 + Satterthwaite
// ---------------------------------------------------------------------------------------------

interface Robusto {
  seRobusto: number[];
  gl: number[];
}

/**
 * CR2 con los grados de libertad de Satterthwaite para cada coeficiente, a partir del ajuste con el modelo de trabajo Φⱼ = Σⱼ (y pesos
 * Wⱼ = Φⱼ⁻¹). Con Hⱼ,c = Xⱼ M X'c Wc (proyección), eⱼ = (I − H)ⱼ· y, y Aⱼ = Dⱼ' Bⱼ^{-1/2} Dⱼ, con Dⱼ'Dⱼ = Φⱼ y
 * Bⱼ = Dⱼ (I − H)ⱼ· Φ (I − H)ⱼ·' Dⱼ'. Varianza: M (Σⱼ Xⱼ'Wⱼ Aⱼ eⱼ eⱼ' Aⱼ Wⱼ Xⱼ) M.
 */
function varianzaRobusta(clusters: Cluster[], Phi: Mat[], g: Gls, p: number): Robusto {
  const m = clusters.length;
  const { M, W, beta } = g;
  // Residuos y pre-cálculos por conglomerado
  const Dj: Mat[] = Phi.map((P) => {
    const L = cholesky(P) ?? P; // P = L L' → D = L'
    return trans(L);
  });
  const XtW: Mat[] = clusters.map((c, j) => mult(trans(c.X), W[j])); // p × n_j
  const Ares: Mat[] = []; // A_j
  const resid: number[][] = clusters.map((c) => c.y.map((yy, a) => yy - c.X[a].reduce((s, x, i) => s + x * beta[i], 0)));
  for (let j = 0; j < m; j++) {
    const nj = clusters[j].idx.length;
    // Σ_c (I−H)_{j,c} Φ_c (I−H)_{j,c}'  con (I−H)_{j,c} = δ_{jc} I − X_j M X_c' W_c
    let B = matriz(nj, nj);
    const XjM = mult(clusters[j].X, M); // n_j × p
    for (let c = 0; c < m; c++) {
      const T = mult(XjM, XtW[c]); // n_j × n_c  (= H_{j,c})
      const IH = T.map((fila, a) => fila.map((v, b) => (j === c ? (a === b ? 1 : 0) : 0) - v));
      const termino = mult(mult(IH, Phi[c]), trans(IH));
      B = B.map((fila, a) => fila.map((v, b) => v + termino[a][b]));
    }
    const Bd = mult(mult(Dj[j], B), trans(Dj[j]));
    // simetrizar por redondeo
    const Bs = Bd.map((fila, a) => fila.map((v, b) => 0.5 * (v + Bd[b][a])));
    Ares.push(mult(mult(trans(Dj[j]), raizInversaSimetrica(Bs)), Dj[j]));
  }
  // Varianza sándwich
  let meat = matriz(p, p);
  const q: number[][][] = []; // para Satterthwaite: q_j(coef) = A_j W_j X_j M e_coef
  for (let j = 0; j < m; j++) {
    const AWX = mult(mult(Ares[j], W[j]), clusters[j].X); // n_j × p
    const Ae = multVec(Ares[j], resid[j]);
    const We = multVec(W[j], Ae);
    const s = multVec(trans(clusters[j].X), We); // X_j' W_j A_j e_j  (p)
    meat = meat.map((fila, a) => fila.map((v, b) => v + s[a] * s[b]));
    q.push(mult(AWX, M)); // n_j × p : columnas = q_j para cada coeficiente
  }
  const V = mult(mult(M, meat), M);
  const seRobusto = V.map((fila, i) => Math.sqrt(Math.max(0, fila[i])));
  // Satterthwaite: b_j = [δ_{jc} q_j − W_c X_c u_j]_c  con u_j = M X_j' q_j
  const gl: number[] = [];
  for (let k = 0; k < p; k++) {
    const bj: number[][][] = []; // [j][c] vector de longitud n_c
    for (let j = 0; j < m; j++) {
      const qj = q[j].map((fila) => fila[k]); // n_j
      const u = multVec(M, multVec(trans(clusters[j].X), qj)); // p
      bj.push(clusters.map((cl, c) => {
        const WXu = multVec(mult(W[c], cl.X), u);
        return WXu.map((v, a) => (j === c ? qj[a] : 0) - v);
      }));
    }
    // G_{jl} = Σ_c b_{j,c}' Φ_c b_{l,c}
    const PhiB: number[][][] = bj.map((fila) => fila.map((v, c) => multVec(Phi[c], v)));
    let num = 0, den = 0;
    for (let j = 0; j < m; j++) {
      for (let l = j; l < m; l++) {
        let G = 0;
        for (let c = 0; c < m; c++) for (let a = 0; a < bj[j][c].length; a++) G += bj[j][c][a] * PhiB[l][c][a];
        if (j === l) { num += G; den += G * G; } else den += 2 * G * G;
      }
    }
    gl.push(den > 0 ? (num * num) / den : NaN);
  }
  return { seRobusto, gl };
}

// ---------------------------------------------------------------------------------------------
// Ajuste completo
// ---------------------------------------------------------------------------------------------

/**
 * Resultado de la parte de varianza robusta con los componentes de varianza dados (los que ajusta `ajustarMultinivel`, o los de otro
 * programa para comprobarla). Expuesta para validar el estimador sándwich de forma aislada.
 */
export function inferenciaRobusta(
  datos: EfectoDependiente[],
  sigma2Estudio: number,
  sigma2Efecto: number,
  opc: OpcionesMultinivel = {},
): { beta: number[]; seModelo: number[]; seRobusto: number[]; gl: number[] } | null {
  const rho = opc.rho ?? 0.6;
  const p = 1 + (opc.nombresModeradores?.length ?? 0);
  const X = datos.map((d) => [1, ...(d.x ?? [])]);
  const clusters = construirClusters(datos, X, rho);
  const g = gls(clusters, sigma2Estudio, sigma2Efecto, p);
  if (!g) return null;
  const Phi = clusters.map((c) => sigmaCluster(c, sigma2Estudio, sigma2Efecto));
  const rb = varianzaRobusta(clusters, Phi, g, p);
  return { beta: g.beta, seModelo: g.M.map((f, i) => Math.sqrt(f[i])), seRobusto: rb.seRobusto, gl: rb.gl };
}

/**
 * Ajusta el modelo CHE por REML y calcula la inferencia robusta. `null` con menos de 3 efectos, con menos de 2 estudios o si el ajuste no se
 * puede calcular. La escala es la de los `y` (log para OR/RR/HR…); la conversión a la escala natural la hace quien llama.
 */
export function ajustarMultinivel(datos: EfectoDependiente[], opc: OpcionesMultinivel = {}): ResultadoMultinivel | null {
  const rho = opc.rho ?? 0.6;
  if (!(rho >= 0 && rho < 1)) return null;
  const validos = datos.filter((d) => Number.isFinite(d.y) && Number.isFinite(d.v) && d.v > 0 && d.estudio !== "");
  const nombresMod = opc.nombresModeradores ?? [];
  const p = 1 + nombresMod.length;
  if (validos.length < 3 || validos.length <= p) return null;
  if (validos.some((d) => (d.x ?? []).length !== nombresMod.length || (d.x ?? []).some((v) => !Number.isFinite(v)))) return null;
  const X = validos.map((d) => [1, ...(d.x ?? [])]);
  const clusters = construirClusters(validos, X, rho);
  const m = clusters.length;
  if (m < 2) return null;
  const maxPorEstudio = Math.max(...clusters.map((c) => c.idx.length));
  const estudiosConVarios = clusters.filter((c) => c.idx.length > 1).length;

  // Punto de partida: DerSimonian-Laird sobre todos los efectos como si fueran independientes
  const w = validos.map((d) => 1 / d.v);
  const sw = w.reduce((a, b) => a + b, 0);
  const media = validos.reduce((s, d, i) => s + w[i] * d.y, 0) / sw;
  const q = validos.reduce((s, d, i) => s + w[i] * (d.y - media) ** 2, 0);
  const c0 = sw - w.reduce((s, x) => s + x * x, 0) / sw;
  const t20 = Math.max(0.01, (q - (validos.length - 1)) / c0);

  const identificable = estudiosConVarios > 0;
  const objetivo = (x: number[]) => {
    const v = objetivoREML(clusters, identificable ? x[0] * x[0] : 0, x[1] * x[1], p);
    return Number.isFinite(v) ? v : 1e12;
  };
  let mejor: { x: number[]; fx: number } | null = null;
  const inicios = [[Math.sqrt(t20 / 2), Math.sqrt(t20 / 2)], [Math.sqrt(t20), 0.01], [0.01, Math.sqrt(t20)], [0.3, 0.3]];
  for (const x0 of inicios) {
    let sol = nelderMead(objetivo, x0, { paso: 0.2, maxEval: 1500, reinicios: 3, tolX: 1e-9, tolF: 1e-12 });
    sol = nelderMead(objetivo, sol.x, { paso: 0.01, maxEval: 1500, reinicios: 4, tolX: 1e-12, tolF: 1e-14 });
    if (!mejor || sol.fx < mejor.fx - 1e-12) mejor = sol;
  }
  if (!mejor) return null;
  const reintento = nelderMead(objetivo, mejor.x, { paso: 0.002, maxEval: 600, reinicios: 0, tolX: 1e-12, tolF: 1e-14 });
  const estable = Math.abs(reintento.fx - mejor.fx) < 1e-9;
  if (reintento.fx < mejor.fx) mejor = reintento;
  let s2u = identificable ? mejor.x[0] ** 2 : 0;
  let s2w = mejor.x[1] ** 2;
  const limites: ("estudio" | "efecto")[] = [];
  if (s2u < SIGMA2_MIN) { s2u = 0; if (identificable) limites.push("estudio"); }
  if (s2w < SIGMA2_MIN) { s2w = 0; limites.push("efecto"); }

  const g = gls(clusters, s2u, s2w, p);
  if (!g) return null;
  const Phi = clusters.map((c) => sigmaCluster(c, s2u, s2w));
  const rb = varianzaRobusta(clusters, Phi, g, p);
  const z = normalCuantilDosColas(0.05);
  const N = validos.length;
  const nombres = ["Efecto combinado (intercepto)", ...nombresMod];
  const coeficientes: CoeficienteMultinivel[] = g.beta.map((b, i) => {
    const seM = Math.sqrt(g.M[i][i]);
    const se = rb.seRobusto[i];
    const gl = rb.gl[i];
    const t = b / se;
    const tcrit = Number.isFinite(gl) && gl > 0 ? tStudentCuantilDosColas(0.05, gl) : NaN;
    return {
      nombre: nombres[i],
      estimacion: b,
      seModelo: seM,
      seRobusto: se,
      gl,
      t,
      pValor: Number.isFinite(gl) && gl > 0 ? tStudentDosColas(t, gl) : normalDosColas(t),
      icInferior: b - tcrit * se,
      icSuperior: b + tcrit * se,
      icModeloInferior: b - z * seM,
      icModeloSuperior: b + z * seM,
    };
  });
  // Log-verosimilitud REML con constantes: −½(N−p)ln 2π + ½ ln|X'X| − ½[Σ ln|Σⱼ| + ln|X'WX| + r'Wr]
  const XtX = X.reduce((A, f) => A.map((fila, a) => fila.map((v, b) => v + f[a] * f[b])), matriz(p, p));
  const detXtX = inversaSPD(XtX)?.logDet ?? 0;
  const ll = -0.5 * (N - p) * Math.log(2 * Math.PI) + 0.5 * detXtX - 0.5 * (g.sumaLogDet + g.logDetXWX + g.residuo2);

  const advertencias: string[] = [];
  if (m < 10) advertencias.push(`Con solo ${m} estudios, los componentes de varianza y la varianza robusta se estiman con poca precisión.`);
  const glMin = Math.min(...rb.gl.filter(Number.isFinite));
  if (Number.isFinite(glMin) && glMin < 4) advertencias.push(`Los grados de libertad de Satterthwaite son ${glMin.toFixed(1)} (< 4): el intervalo y el valor p robustos no son fiables.`);
  if (!identificable) advertencias.push("Ningún estudio aporta más de un efecto: el modelo se reduce a efectos aleatorios y no hay dependencia que modelar.");
  if (!estable) advertencias.push("La optimización de los componentes de varianza no mostró convergencia estricta: contrasta con otro programa antes de publicar.");
  return {
    k: N,
    m,
    estudiosConVarios,
    maxPorEstudio,
    rho,
    sigma2Estudio: s2u,
    sigma2Efecto: s2w,
    coeficientes,
    logVerosimilitudREML: ll,
    convergio: estable,
    limites,
    advertencias,
  };
}

/** Sensibilidad a la correlación supuesta: el mismo ajuste con distintas ρ. */
export const RHO_SENSIBILIDAD: readonly number[] = [0, 0.2, 0.4, 0.6, 0.8];

/** Repite el ajuste multinivel con varias correlaciones supuestas entre efectos (por defecto 0; 0,2; 0,4; 0,6 y 0,8) para ver cuánto depende el resultado de esa suposición. */
export function sensibilidadRho(datos: EfectoDependiente[], opc: OpcionesMultinivel = {}, rhos: readonly number[] = RHO_SENSIBILIDAD): { rho: number; resultado: ResultadoMultinivel | null }[] {
  return rhos.map((rho) => ({ rho, resultado: ajustarMultinivel(datos, { ...opc, rho }) }));
}
