// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * **Modelo bivariado de precisión diagnóstica** (Reitsma y cols. 2005; Chu y Cole 2006), equivalente al HSROC de Rutter y Gatsonis cuando
 * no hay covariables (Harbord y cols. 2007): la sensibilidad y la especificidad de cada estudio, en escala logit, siguen una normal
 * bivariada entre estudios, con medias (μA, μB), desviaciones (σA, σB) y correlación ρ. Es lo que recomienda el Cochrane Handbook de
 * revisiones de precisión diagnóstica (cap. 10.5.2) para resumir una prueba, y de él salen el punto resumen, su región de confianza, la
 * región de predicción y la curva SROC.
 *
 * Dos formas de ajustarlo:
 *  · **«binomial»** (la principal): verosimilitud binomial exacta dentro de cada estudio e integración numérica sobre los dos efectos
 *    aleatorios (cuadratura de Gauss-Hermite adaptativa en 2 dimensiones), máxima verosimilitud. No necesita corregir celdas en cero.
 *    Es el modelo de `metandi` (Stata) y de SAS NLMIXED, el que el manual Cochrane recomienda frente a la aproximación normal.
 *  · **«normal»**: la aproximación de Reitsma (2005): las proporciones se pasan a logit con la varianza 1/a + 1/b (con +0,5 en las
 *    cuatro celdas de un estudio si alguna es cero) y se ajusta un modelo de efectos aleatorios con REML. Es `mada::reitsma`. El manual
 *    Cochrane advierte de que esta aproximación puede funcionar mal con pocos eventos; aquí se usa para comparar y como punto de partida.
 *
 * Sin dependencias de la app (solo `funciones-especiales`). Lógica pura.
 *
 * Fuentes de las fórmulas: Reitsma 2005; Chu y Cole 2006; Harbord 2007 (equivalencia con el HSROC); Cochrane Handbook DTA cap. 10.5.2;
 * Golub y Welsch 1969 (nodos y pesos de Gauss-Hermite por los valores propios de la matriz de Jacobi); Nelder y Mead 1965.
 * Detalle de procedencia y de la validación contra R en docs/PROCEDENCIA.md y VALIDACION_ESTADISTICA.md.
 */
import { normalCuantilDosColas } from "./funciones-especiales.ts";
import {
  hessianoNumerico,
  invertir,
  limitar,
  lnChoose,
  logit,
  logSigmoide,
  nelderMead,
  nodosGaussHermite,
  sigmoide,
  type NodosGH,
} from "./numerica-compartida.ts";


export { nodosGaussHermite };

// ---------------------------------------------------------------------------------------------
// Tipos
// ---------------------------------------------------------------------------------------------

/** Tabla 2×2 de un estudio: verdaderos y falsos positivos y negativos (conteos originales, sin corregir). */
export interface TablaDiagnostica {
  vp: number;
  fp: number;
  fn: number;
  vn: number;
}

export type MetodoBivariado = "binomial" | "normal";

export interface ParametrosBivariado {
  /** logit de la sensibilidad media (μA). */
  muSens: number;
  /** logit de la especificidad media (μB). */
  muEsp: number;
  /** Desviación estándar entre estudios del logit de la sensibilidad (σA). */
  sdSens: number;
  /** Desviación estándar entre estudios del logit de la especificidad (σB). */
  sdEsp: number;
  /** Correlación entre los logit de la sensibilidad y de la especificidad (ρ_AB). Suele ser negativa si hay efecto de umbral. */
  rho: number;
}

export interface IntervaloProporcion {
  estimacion: number;
  icInferior: number;
  icSuperior: number;
}

export interface PuntoROC {
  /** Tasa de falsos positivos (1 − especificidad). */
  fpr: number;
  /** Sensibilidad. */
  tpr: number;
}

export interface ResultadoBivariado {
  metodo: MetodoBivariado;
  /** «ML» (binomial) o «REML» (normal). */
  estimacion: "ML" | "REML";
  k: number;
  parametros: ParametrosBivariado;
  /** Matriz de covarianza de (μA, μB): la de las estimaciones de las medias, no la variación entre estudios. */
  vcovMu: [[number, number], [number, number]];
  sensibilidad: IntervaloProporcion;
  especificidad: IntervaloProporcion;
  /** Tasa de falsos positivos = 1 − especificidad (su IC es el de la especificidad invertido). */
  fpr: IntervaloProporcion;
  /** Log-verosimilitud en el óptimo (binomial: con los coeficientes binomiales; normal: la REML sin constantes). */
  logVerosimilitud: number;
  convergio: boolean;
  /** Parámetros que quedaron en el límite del rango permitido (varianza ≈ 0 o |ρ| ≈ 1). */
  limites: string[];
  advertencias: string[];
  /** Rango de FPR observado en los estudios (la curva solo se dibuja ahí: fuera es extrapolar). */
  rangoFprObservado: [number, number];
  /** Área bajo la curva SROC (generalización de Rutter-Gatsonis), integrada en (0, 1) y por tanto con extrapolación. */
  auc: number;
}

// ---------------------------------------------------------------------------------------------
// Utilidades numéricas
// ---------------------------------------------------------------------------------------------

/** Límites del espacio de búsqueda (en las parametrizaciones sin restricciones). */
const LIM = { logSd: [Math.log(1e-3), Math.log(30)] as [number, number], atanhRho: [-6.9, 6.9] as [number, number], mu: [-14, 14] as [number, number] };

// ---------------------------------------------------------------------------------------------
// Aproximación normal (Reitsma 2005), REML
// ---------------------------------------------------------------------------------------------

interface DatosLogit {
  /** [logit(sens), logit(esp)] por estudio. */
  y: [number, number][];
  /** Varianzas dentro de cada estudio [1/vp + 1/fn, 1/vn + 1/fp] (sin covarianza). */
  v: [number, number][];
}

/** Pasa las tablas a logit con +0,5 en las cuatro celdas de un estudio si alguna es cero (la misma regla del resto del SROC y de `mada` con `correction.control = "single"`). */
function datosLogit(tablas: TablaDiagnostica[]): DatosLogit {
  const y: [number, number][] = [];
  const v: [number, number][] = [];
  for (const t of tablas) {
    const c = t.vp === 0 || t.fp === 0 || t.fn === 0 || t.vn === 0 ? 0.5 : 0;
    const vp = t.vp + c, fn = t.fn + c, fp = t.fp + c, vn = t.vn + c;
    y.push([Math.log(vp / fn), Math.log(vn / fp)]);
    v.push([1 / vp + 1 / fn, 1 / vn + 1 / fp]);
  }
  return { y, v };
}

const covarianza = (sdA: number, sdB: number, rho: number): [number, number, number] => [sdA * sdA, rho * sdA * sdB, sdB * sdB];

/**
 * Objetivo REML (a minimizar; sin constantes): ½[Σ ln|Vᵢ| + ln|Σ Wᵢ| + Σ rᵢ'Wᵢrᵢ], con Vᵢ = diag(vᵢ) + Σ, Wᵢ = Vᵢ⁻¹ y rᵢ = yᵢ − μ̂,
 * donde μ̂ es la media de mínimos cuadrados generalizados. Devuelve también μ̂ y la covarianza de μ̂ = (Σ Wᵢ)⁻¹.
 */
function remlNormal(d: DatosLogit, sdA: number, sdB: number, rho: number) {
  const [s11, s12, s22] = covarianza(sdA, sdB, rho);
  let logDet = 0;
  let a11 = 0, a12 = 0, a22 = 0; // Σ Wᵢ
  let b1 = 0, b2 = 0; // Σ Wᵢ yᵢ
  const Ws: [number, number, number][] = [];
  for (let i = 0; i < d.y.length; i++) {
    const v11 = d.v[i][0] + s11;
    const v22 = d.v[i][1] + s22;
    const det = v11 * v22 - s12 * s12;
    if (!(det > 0)) return null;
    logDet += Math.log(det);
    const w11 = v22 / det, w12 = -s12 / det, w22 = v11 / det;
    Ws.push([w11, w12, w22]);
    a11 += w11; a12 += w12; a22 += w22;
    b1 += w11 * d.y[i][0] + w12 * d.y[i][1];
    b2 += w12 * d.y[i][0] + w22 * d.y[i][1];
  }
  const detA = a11 * a22 - a12 * a12;
  if (!(detA > 0)) return null;
  const mu: [number, number] = [(a22 * b1 - a12 * b2) / detA, (-a12 * b1 + a11 * b2) / detA];
  let q = 0;
  for (let i = 0; i < d.y.length; i++) {
    const r1 = d.y[i][0] - mu[0];
    const r2 = d.y[i][1] - mu[1];
    const [w11, w12, w22] = Ws[i];
    q += w11 * r1 * r1 + 2 * w12 * r1 * r2 + w22 * r2 * r2;
  }
  return { objetivo: 0.5 * (logDet + Math.log(detA) + q), mu, vcov: [[a22 / detA, -a12 / detA], [-a12 / detA, a11 / detA]] as [[number, number], [number, number]] };
}

/** Parámetros no restringidos → (σA, σB, ρ), limitados al rango permitido. */
function desempaquetarVar(t: number[]): { sdA: number; sdB: number; rho: number } {
  return { sdA: Math.exp(limitar(t[0], ...LIM.logSd)), sdB: Math.exp(limitar(t[1], ...LIM.logSd)), rho: Math.tanh(limitar(t[2], ...LIM.atanhRho)) };
}

function ajusteNormal(d: DatosLogit): { sdA: number; sdB: number; rho: number; mu: [number, number]; vcov: [[number, number], [number, number]]; objetivo: number; convergio: boolean; theta: number[] } | null {
  const k = d.y.length;
  const media = (xs: number[]) => xs.reduce((s, x) => s + x, 0) / xs.length;
  const varMuestral = (xs: number[]) => {
    const m = media(xs);
    return xs.reduce((s, x) => s + (x - m) ** 2, 0) / Math.max(1, xs.length - 1);
  };
  const tauA0 = Math.max(0.05, varMuestral(d.y.map((p) => p[0])) - media(d.v.map((p) => p[0])));
  const tauB0 = Math.max(0.05, varMuestral(d.y.map((p) => p[1])) - media(d.v.map((p) => p[1])));
  const f = (t: number[]) => {
    const { sdA, sdB, rho } = desempaquetarVar(t);
    const r = remlNormal(d, sdA, sdB, rho);
    return r ? r.objetivo : 1e12;
  };
  let mejor: { x: number[]; fx: number; convergio: boolean } | null = null;
  for (const z0 of [0, 0.8, -0.8]) {
    const sol = nelderMead(f, [Math.log(Math.sqrt(tauA0)), Math.log(Math.sqrt(tauB0)), z0], { paso: 0.7, maxEval: 3000, reinicios: 4 });
    if (!mejor || sol.fx < mejor.fx - 1e-12) mejor = { x: sol.x, fx: sol.fx, convergio: sol.convergio };
  }
  if (!mejor || k < 2) return null;
  const { sdA, sdB, rho } = desempaquetarVar(mejor.x);
  const r = remlNormal(d, sdA, sdB, rho);
  if (!r) return null;
  return { sdA, sdB, rho, mu: r.mu, vcov: r.vcov, objetivo: r.objetivo, convergio: mejor.convergio, theta: mejor.x };
}

// ---------------------------------------------------------------------------------------------
// Modelo binomial-normal (verosimilitud exacta, cuadratura adaptativa)
// ---------------------------------------------------------------------------------------------

/**
 * Logaritmo de ∫∫ Bin(vp; vp+fn, σ(a))·Bin(vn; vn+fp, σ(b))·N₂((a,b); μ, Σ) da db para un estudio, sin los coeficientes binomiales.
 * Cuadratura de Gauss-Hermite adaptativa: se centra en la moda posterior (Newton) y se escala con su curvatura.
 */
function logIntegralEstudio(t: TablaDiagnostica, mu: [number, number], sdA: number, sdB: number, rho: number, nodos: NodosGH): number {
  const [s11, s12, s22] = covarianza(sdA, sdB, rho);
  const det = s11 * s22 - s12 * s12;
  const p11 = s22 / det, p12 = -s12 / det, p22 = s11 / det; // Σ⁻¹
  const n1 = t.vp + t.fn; // enfermos
  const n0 = t.vn + t.fp; // sanos
  const logNorm = -Math.log(2 * Math.PI) - 0.5 * Math.log(det);
  const g = (a: number, b: number) => {
    const da = a - mu[0], db = b - mu[1];
    return t.vp * logSigmoide(a) + t.fn * logSigmoide(-a) + t.vn * logSigmoide(b) + t.fp * logSigmoide(-b) +
      logNorm - 0.5 * (p11 * da * da + 2 * p12 * da * db + p22 * db * db);
  };
  // Moda por Newton con retroceso (la función log-posterior es cóncava).
  let a = mu[0], b = mu[1];
  let gv = g(a, b);
  let h11 = 0, h12 = 0, h22 = 0;
  for (let it = 0; it < 60; it++) {
    const pa = sigmoide(a), pb = sigmoide(b);
    const da = a - mu[0], db = b - mu[1];
    const g1 = t.vp - n1 * pa - (p11 * da + p12 * db);
    const g2 = t.vn - n0 * pb - (p12 * da + p22 * db);
    h11 = -n1 * pa * (1 - pa) - p11;
    h12 = -p12;
    h22 = -n0 * pb * (1 - pb) - p22;
    const dh = h11 * h22 - h12 * h12;
    const sa = -(h22 * g1 - h12 * g2) / dh;
    const sb = -(-h12 * g1 + h11 * g2) / dh;
    let paso = 1;
    let na = a + sa, nb = b + sb, ng = g(na, nb);
    while (!(ng >= gv - 1e-14) && paso > 1e-8) {
      paso /= 2;
      na = a + paso * sa;
      nb = b + paso * sb;
      ng = g(na, nb);
    }
    const cambio = Math.abs(na - a) + Math.abs(nb - b);
    a = na; b = nb; gv = ng;
    if (cambio < 1e-11) break;
  }
  const pa = sigmoide(a), pb = sigmoide(b);
  h11 = -n1 * pa * (1 - pa) - p11;
  h12 = -p12;
  h22 = -n0 * pb * (1 - pb) - p22;
  // Covarianza de la aproximación gaussiana: (−H)⁻¹ y su factor de Cholesky inferior.
  const m11 = -h11, m12 = -h12, m22 = -h22;
  const dm = m11 * m22 - m12 * m12;
  const c11 = m22 / dm, c12 = -m12 / dm, c22 = m11 / dm;
  const l11 = Math.sqrt(c11);
  const l21 = c12 / l11;
  const l22 = Math.sqrt(c22 - l21 * l21);
  const q = nodos.z.length;
  let maxTermino = -Infinity;
  const terminos: number[] = [];
  for (let j = 0; j < q; j++) {
    for (let l = 0; l < q; l++) {
      const z1 = nodos.z[j], z2 = nodos.z[l];
      const aa = a + l11 * z1;
      const bb = b + l21 * z1 + l22 * z2;
      const tt = g(aa, bb) + (z1 * z1 + z2 * z2) / 2 + nodos.logw[j] + nodos.logw[l];
      terminos.push(tt);
      if (tt > maxTermino) maxTermino = tt;
    }
  }
  let suma = 0;
  for (const tt of terminos) suma += Math.exp(tt - maxTermino);
  return Math.log(2 * Math.PI) + Math.log(l11 * l22) + maxTermino + Math.log(suma);
}

/** Log-verosimilitud marginal del modelo binomial-normal (con los coeficientes binomiales). */
export function logVerosimilitudBinomial(
  tablas: TablaDiagnostica[],
  p: ParametrosBivariado,
  nodosPorDimension = 15,
): number {
  const nodos = nodosGaussHermite(nodosPorDimension);
  let total = 0;
  for (const t of tablas) {
    total += logIntegralEstudio(t, [p.muSens, p.muEsp], p.sdSens, p.sdEsp, p.rho, nodos);
    total += lnChoose(t.vp + t.fn, t.vp) + lnChoose(t.vn + t.fp, t.vn);
  }
  return total;
}

const empaquetar = (p: ParametrosBivariado): number[] => [p.muSens, p.muEsp, Math.log(p.sdSens), Math.log(p.sdEsp), Math.atanh(p.rho)];
const desempaquetar = (x: number[]): ParametrosBivariado => ({
  muSens: limitar(x[0], ...LIM.mu),
  muEsp: limitar(x[1], ...LIM.mu),
  sdSens: Math.exp(limitar(x[2], ...LIM.logSd)),
  sdEsp: Math.exp(limitar(x[3], ...LIM.logSd)),
  rho: Math.tanh(limitar(x[4], ...LIM.atanhRho)),
});

// ---------------------------------------------------------------------------------------------
// Salidas comunes: curva, AUC, regiones
// ---------------------------------------------------------------------------------------------

/**
 * Curva SROC implícita en el modelo: en el espacio logit-logit es la recta que pasa por el punto resumen con pendiente σA/σB
 * (la «generalización de Rutter-Gatsonis», coincide con el HSROC cuando hay equivalencia con el bivariado: Harbord 2007):
 *   logit(sensibilidad) = μA + (σA/σB)·(logit(FPR) − logit(FPR̄)),  con logit(FPR̄) = −μB.
 */
export function sensibilidadEnCurva(p: ParametrosBivariado, fpr: number): number {
  return sigmoide(p.muSens + (p.sdSens / p.sdEsp) * (logit(fpr) + p.muEsp));
}

/** Curva SROC en `n` valores de FPR entre `desde` y `hasta`, equiespaciados en la escala logit (más puntos donde la curva se dobla). */
export function curvaSROC(p: ParametrosBivariado, desde: number, hasta: number, n = 80): PuntoROC[] {
  const l0 = logit(limitar(desde, 1e-6, 1 - 1e-6));
  const l1 = logit(limitar(hasta, 1e-6, 1 - 1e-6));
  return Array.from({ length: n }, (_, i) => {
    const fpr = sigmoide(l0 + ((l1 - l0) * i) / (n - 1));
    return { fpr, tpr: sensibilidadEnCurva(p, fpr) };
  });
}

/**
 * Área bajo la curva SROC, ∫₀¹ TPR(FPR) dFPR, por la regla del trapecio en la escala logit (FPR = σ(u), dFPR = σ(u)(1−σ(u)) du) entre
 * u = −40 y 40 con 8000 intervalos: error ≈ 1e-10. Incluye la extrapolación fuera del rango de FPR observado.
 */
export function areaBajoCurva(p: ParametrosBivariado): number {
  const n = 8000;
  const u0 = -40, u1 = 40;
  const h = (u1 - u0) / n;
  const integrando = (u: number) => {
    const fpr = sigmoide(u);
    return sigmoide(p.muSens + (p.sdSens / p.sdEsp) * (u + p.muEsp)) * fpr * (1 - fpr);
  };
  let s = 0.5 * (integrando(u0) + integrando(u1));
  for (let i = 1; i < n; i++) s += integrando(u0 + i * h);
  return s * h;
}

/**
 * Elipse de nivel `nivel` en el plano logit (centrada en `centro`, con matriz de covarianza `cov`) pasada al espacio ROC
 * (FPR, sensibilidad): d' cov⁻¹ d = χ²₂(nivel) = −2 ln(1 − nivel). Los puntos están en el orden del ángulo.
 * Coordenadas logit: (logit sens, logit esp); FPR = σ(−logit esp).
 */
export function elipseROC(centro: [number, number], cov: [[number, number], [number, number]], nivel = 0.95, n = 100): PuntoROC[] {
  const q = -2 * Math.log(1 - nivel);
  const l11 = Math.sqrt(cov[0][0]);
  const l21 = cov[1][0] / l11;
  const l22 = Math.sqrt(Math.max(0, cov[1][1] - l21 * l21));
  return Array.from({ length: n }, (_, i) => {
    const ang = (2 * Math.PI * i) / n;
    const c = Math.cos(ang) * Math.sqrt(q);
    const s = Math.sin(ang) * Math.sqrt(q);
    const sens = centro[0] + l11 * c;
    const esp = centro[1] + l21 * c + l22 * s;
    return { fpr: sigmoide(-esp), tpr: sigmoide(sens) };
  });
}

/** Región de confianza (para el punto resumen) y de predicción (para un estudio nuevo) en el espacio ROC. */
export function regionesROC(r: ResultadoBivariado, nivel = 0.95): { confianza: PuntoROC[]; prediccion: PuntoROC[] } {
  const { muSens, muEsp, sdSens, sdEsp, rho } = r.parametros;
  const [t11, t12, t22] = covarianza(sdSens, sdEsp, rho);
  const pred: [[number, number], [number, number]] = [
    [r.vcovMu[0][0] + t11, r.vcovMu[0][1] + t12],
    [r.vcovMu[1][0] + t12, r.vcovMu[1][1] + t22],
  ];
  return { confianza: elipseROC([muSens, muEsp], r.vcovMu, nivel), prediccion: elipseROC([muSens, muEsp], pred, nivel) };
}

// ---------------------------------------------------------------------------------------------
// Ajuste completo
// ---------------------------------------------------------------------------------------------

export interface OpcionesAjuste {
  metodo?: MetodoBivariado;
  /** Puntos de cuadratura por dimensión durante la optimización (binomial). */
  nodosOptimizacion?: number;
  /** Puntos de cuadratura por dimensión para el valor final. */
  nodosFinales?: number;
}

const Z95_EXACTO = () => normalCuantilDosColas(0.05);

function resumen(
  metodo: MetodoBivariado,
  tablas: TablaDiagnostica[],
  p: ParametrosBivariado,
  vcovMu: [[number, number], [number, number]],
  logVer: number,
  convergio: boolean,
  limites: string[],
  advertencias: string[],
): ResultadoBivariado {
  const z = Z95_EXACTO();
  const se1 = Math.sqrt(Math.max(0, vcovMu[0][0]));
  const se2 = Math.sqrt(Math.max(0, vcovMu[1][1]));
  const sens: IntervaloProporcion = { estimacion: sigmoide(p.muSens), icInferior: sigmoide(p.muSens - z * se1), icSuperior: sigmoide(p.muSens + z * se1) };
  const esp: IntervaloProporcion = { estimacion: sigmoide(p.muEsp), icInferior: sigmoide(p.muEsp - z * se2), icSuperior: sigmoide(p.muEsp + z * se2) };
  const fprs = tablas.map((t) => (t.fp + (t.vp === 0 || t.fp === 0 || t.fn === 0 || t.vn === 0 ? 0.5 : 0)) / (t.fp + t.vn + (t.vp === 0 || t.fp === 0 || t.fn === 0 || t.vn === 0 ? 1 : 0)));
  return {
    metodo,
    estimacion: metodo === "binomial" ? "ML" : "REML",
    k: tablas.length,
    parametros: p,
    vcovMu,
    sensibilidad: sens,
    especificidad: esp,
    fpr: { estimacion: 1 - esp.estimacion, icInferior: 1 - esp.icSuperior, icSuperior: 1 - esp.icInferior },
    logVerosimilitud: logVer,
    convergio,
    limites,
    advertencias,
    rangoFprObservado: [Math.min(...fprs), Math.max(...fprs)],
    auc: areaBajoCurva(p),
  };
}

function limitesAlcanzados(x: number[]): string[] {
  const out: string[] = [];
  if (x[2] <= LIM.logSd[0] + 1e-3) out.push("sdSens");
  if (x[3] <= LIM.logSd[0] + 1e-3) out.push("sdEsp");
  if (Math.abs(x[4]) >= LIM.atanhRho[1] - 1e-3) out.push("rho");
  return out;
}

/** Valida las tablas: conteos finitos, no negativos y con enfermos y sanos en cada estudio. Devuelve las utilizables. */
export function tablasUtilizables(tablas: TablaDiagnostica[]): TablaDiagnostica[] {
  return tablas.filter(
    (t) => [t.vp, t.fp, t.fn, t.vn].every((x) => Number.isFinite(x) && x >= 0 && Number.isInteger(x)) && t.vp + t.fn > 0 && t.vn + t.fp > 0,
  );
}

/**
 * Ajusta el modelo bivariado. `null` con menos de 4 estudios utilizables o si el ajuste no se puede calcular (matriz singular).
 * Con «binomial», si la optimización no converge cae en la aproximación normal e indica el motivo en `advertencias`.
 */
export function ajustarBivariado(tablasEntrada: TablaDiagnostica[], opc: OpcionesAjuste = {}): ResultadoBivariado | null {
  const tablas = tablasUtilizables(tablasEntrada);
  const k = tablas.length;
  if (k < 4) return null;
  const metodo = opc.metodo ?? "binomial";
  const advertencias: string[] = [];
  if (k < 10) advertencias.push(`Con solo ${k} estudios, las varianzas y la correlación entre estudios se estiman con muy poca precisión: interpreta la región de predicción y la curva con cautela.`);

  const d = datosLogit(tablas);
  const normal = ajusteNormal(d);
  if (!normal) return null;
  const pNormal: ParametrosBivariado = { muSens: normal.mu[0], muEsp: normal.mu[1], sdSens: normal.sdA, sdEsp: normal.sdB, rho: normal.rho };
  if (metodo === "normal") {
    const limites = limitesAlcanzados(normal.theta);
    return resumen("normal", tablas, pNormal, normal.vcov, -normal.objetivo, normal.convergio, limites, advertencias);
  }

  // ---- binomial exacta ----
  const qOpt = opc.nodosOptimizacion ?? 9;
  const qFin = opc.nodosFinales ?? 21;
  const x0 = empaquetar(pNormal);
  const nll = (q: number) => (x: number[]) => {
    const v = -logVerosimilitudBinomial(tablas, desempaquetar(x), q);
    return Number.isFinite(v) ? v : 1e12;
  };
  let sol = nelderMead(nll(qOpt), x0, { paso: 0.3, maxEval: 3500, reinicios: 3, tolX: 1e-7, tolF: 1e-9 });
  sol = nelderMead(nll(qFin), sol.x, { paso: 0.05, maxEval: 2500, reinicios: 4, tolX: 1e-10, tolF: 1e-12 });
  // Criterio de convergencia: un reinicio más desde el óptimo no mejora la verosimilitud (diferencia < 1e-8), no el tamaño del símplex.
  const reintento = nelderMead(nll(qFin), sol.x, { paso: 0.01, maxEval: 600, reinicios: 0, tolX: 1e-10, tolF: 1e-12 });
  const estable = Math.abs(reintento.fx - sol.fx) < 1e-8;
  if (reintento.fx < sol.fx) sol = { ...reintento, convergio: sol.convergio };
  const x = sol.x.map((v, i) => (i === 0 || i === 1 ? limitar(v, ...LIM.mu) : v));
  const limites = limitesAlcanzados(x);
  const p = desempaquetar(x);
  const ll = -nll(qFin)(x);
  // Covarianza de (μA, μB): inversa del Hessiano sobre los parámetros libres (los que no quedaron en el límite, que se fijan: es el
  // submodelo con ρ = ±1 o varianza nula). Sin parámetros en el límite es el Hessiano completo 5×5. Si no es invertible, la condicionada
  // a las varianzas (2×2).
  const fin = nll(qFin);
  const nombres = ["muSens", "muEsp", "sdSens", "sdEsp", "rho"];
  const libres = [0, 1, 2, 3, 4].filter((i) => i < 2 || !limites.includes(nombres[i]));
  const hessLibres = (idx: number[]) => hessianoNumerico((m: number[]) => { const xx = x.slice(); idx.forEach((j, a) => (xx[j] = m[a])); return fin(xx); }, idx.map((j) => x[j]), 2e-3);
  const bloqueMu = (H: number[][]): [[number, number], [number, number]] | null => {
    const inv = invertir(H);
    return inv && inv[0][0] > 0 && inv[1][1] > 0 && inv[0][0] * inv[1][1] - inv[0][1] ** 2 > 0 ? [[inv[0][0], inv[0][1]], [inv[1][0], inv[1][1]]] : null;
  };
  let vcov = bloqueMu(hessLibres(libres));
  if (!vcov && libres.length > 2) vcov = bloqueMu(hessLibres([0, 1]));
  if (!vcov) {
    advertencias.push("No se pudo calcular la covarianza de las medias con el modelo binomial; se muestra la de la aproximación normal.");
    return resumen("normal", tablas, pNormal, normal.vcov, -normal.objetivo, false, limitesAlcanzados(normal.theta), advertencias);
  }
  if (limites.length > 0) advertencias.push("Algún parámetro de variación entre estudios quedó en el límite: los errores estándar de las medias se calcularon fijándolo en ese límite (con él, no se estima su incertidumbre).");
  if (limites.includes("rho")) advertencias.push("La correlación entre sensibilidad y especificidad quedó en ±1 (con tan pocos estudios no se distingue de un caso extremo).");
  if (limites.includes("sdSens") || limites.includes("sdEsp")) advertencias.push("La variación entre estudios de la sensibilidad o de la especificidad es prácticamente nula.");
  return resumen("binomial", tablas, p, vcov, ll, estable, limites, advertencias);
}

const cacheAjustes = new Map<string, ResultadoBivariado | null>();

/**
 * Igual que `ajustarBivariado`, pero recuerda los últimos ajustes por (método + tablas): el Paso 7, el manuscrito y la cobertura piden el
 * mismo SROC varias veces y el ajuste binomial tarda ~1 s con 40 estudios.
 */
export function ajustarBivariadoConCache(tablas: TablaDiagnostica[], metodo: MetodoBivariado = "binomial"): ResultadoBivariado | null {
  const clave = `${metodo}|${tablas.map((t) => `${t.vp},${t.fp},${t.fn},${t.vn}`).join(";")}`;
  if (cacheAjustes.has(clave)) return cacheAjustes.get(clave) ?? null;
  const r = ajustarBivariado(tablas, { metodo });
  if (cacheAjustes.size >= 16) cacheAjustes.delete(cacheAjustes.keys().next().value as string);
  cacheAjustes.set(clave, r);
  return r;
}

/** Aproximación normal de Reitsma (REML) para comparar con el modelo binomial. */
export const ajustarReitsmaNormal = (tablas: TablaDiagnostica[]): ResultadoBivariado | null => ajustarBivariado(tablas, { metodo: "normal" });
