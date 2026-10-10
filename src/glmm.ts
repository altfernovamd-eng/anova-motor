// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * **Modelos lineales generalizados mixtos (GLMM) para datos binarios** (Stijnen, Hamza y Özdemir 2010; Jackson y cols. 2018):
 * en vez de pasar cada estudio a un efecto (con corrección de continuidad) y combinar esos efectos con pesos, se modela el recuento de
 * eventos directamente con su distribución exacta y un efecto aleatorio entre estudios. Es la alternativa que la literatura propone para
 * eventos raros, muestras pequeñas y proporciones (donde el doble arcoseno puede engañar: Schwarzer y cols. 2019).
 *
 * Dos modelos, los dos de un solo efecto aleatorio, con máxima verosimilitud:
 *  · **Proporciones, binomial-normal:**  xᵢ ~ Binomial(nᵢ, pᵢ),  logit(pᵢ) = μ + uᵢ,  uᵢ ~ N(0, τ²).
 *  · **Odds ratio, hipergeométrica-normal** (verosimilitud condicional exacta; «CM.EL» en `metafor::rma.glmm`):  dado el total de
 *    eventos sᵢ de cada tabla 2×2, los eventos del brazo experimental siguen la hipergeométrica no central con log-OR θᵢ,
 *    θᵢ ~ N(θ, τ²). No tiene parámetros de molestia (el riesgo basal de cada estudio desaparece al condicionar), no necesita corrección
 *    de continuidad y los estudios sin información (sin eventos en ningún brazo, o con todos los participantes con evento) aportan 1 a
 *    la verosimilitud: se declaran, no se «corrigen».
 *
 * La verosimilitud marginal de cada estudio, ∫ P(xᵢ | η = μ + τz) φ(z) dz, se integra por cuadratura de Gauss-Hermite adaptativa de
 * una dimensión (se centra en el máximo del integrando y se escala con su curvatura; el logaritmo de P es cóncavo en η en las dos
 * familias, así que el máximo es único). Los parámetros (μ, ln τ) se estiman con Nelder-Mead y el error estándar sale del Hessiano
 * numérico completo (inversa de la información observada), no condicionado a τ. Si τ queda en el límite inferior se declara y el
 * error estándar de μ se calcula con τ fijo.
 *
 * Sin dependencias de la app (solo `funciones-especiales` y `numerica-compartida`). Lógica pura.
 *
 * Fuentes de las fórmulas: Stijnen y cols. 2010 (hipergeométrica-normal y binomial-normal; comparación de siete modelos para el OR en Jackson y cols. 2018); la hipergeométrica no central es la familia
 * exponencial de Fisher (su derivada primera es x − E[X] y la segunda −Var[X]); Golub y Welsch 1969; Nelder y Mead 1965.
 * Detalle de procedencia y de la validación contra R en docs/PROCEDENCIA.md y VALIDACION_ESTADISTICA.md.
 */
import { normalCuantilDosColas, normalDosColas } from "./funciones-especiales.ts";
import { hessianoNumerico, invertir, limitar, lnChoose, logSigmoide, nelderMead, nodosGaussHermite, sigmoide } from "./numerica-compartida.ts";

// ---------------------------------------------------------------------------------------------
// Tipos
// ---------------------------------------------------------------------------------------------

/** Recuentos de una proporción: eventos y total. */
export interface DatosProporcion {
  etiqueta?: string;
  x: number;
  n: number;
}

/** Tabla 2×2 con el brazo experimental primero: eventos y total de cada brazo. */
export interface DatosTabla {
  etiqueta?: string;
  eventosExp: number;
  totalExp: number;
  eventosCtrl: number;
  totalCtrl: number;
}

export type ModeloGLMM = "binomial-normal" | "hipergeometrica-normal";

export interface ResultadoGLMM {
  modelo: ModeloGLMM;
  /** Estudios que entran en el análisis (recuentos válidos); incluye los que no aportan información. */
  k: number;
  /** Estudios que sí aportan información a la verosimilitud (en la hipergeométrica, los que tienen más de un valor posible). */
  kInformativos: number;
  /** Estudios descartados por datos no válidos (recuentos no enteros, negativos o imposibles), con el motivo. */
  excluidos: { etiqueta: string; motivo: string }[];
  /** Estudios válidos que no aportan información (sin eventos en ningún brazo, o todos con evento): no se corrigen, se declaran. */
  sinInformacion: string[];
  /** Efecto combinado en la escala del modelo: logit de la proporción o ln(OR). */
  mu: number;
  se: number;
  /** IC 95 % de Wald en la escala del modelo. */
  muIc: [number, number];
  /** Efecto combinado en la escala natural (proporción u OR) con su IC 95 %. */
  estimacion: number;
  icInferior: number;
  icSuperior: number;
  /** Valor p de Wald para H₀: ln(OR) = 0; `null` en proporciones (no hay un nulo natural). */
  pValor: number | null;
  tau2: number;
  tau: number;
  logVerosimilitud: number;
  convergio: boolean;
  /** Límites alcanzados: «tau» (τ prácticamente nula) o «tauAlto» (τ en el tope de 30: heterogeneidad no estimable). */
  limites: string[];
  advertencias: string[];
}

// ---------------------------------------------------------------------------------------------
// Familias: logaritmo de la probabilidad del dato del estudio i dado el predictor lineal η, y sus dos primeras derivadas en η
// ---------------------------------------------------------------------------------------------

interface Familia {
  /** log P(datos del estudio | η). */
  logP(eta: number): number;
  /** d/dη log P. */
  d1(eta: number): number;
  /** d²/dη² log P (≤ 0). */
  d2(eta: number): number;
}

function familiaBinomial(x: number, n: number): Familia {
  const lc = lnChoose(n, x);
  return {
    logP: (eta) => lc + x * logSigmoide(eta) + (n - x) * logSigmoide(-eta),
    d1: (eta) => x - n * sigmoide(eta),
    d2: (eta) => {
      const p = sigmoide(eta);
      return -n * p * (1 - p);
    },
  };
}

/** Hipergeométrica no central de la tabla (a, n1, c, n2): P(X = a | total s = a + c) con parámetro de log-OR η. */
function familiaHipergeometrica(a: number, n1: number, c: number, n2: number): Familia {
  const s = a + c;
  const lo = Math.max(0, s - n2);
  const hi = Math.min(n1, s);
  const u: number[] = [];
  const lw: number[] = []; // ln[C(n1,u) C(n2,s−u)]
  for (let v = lo; v <= hi; v++) {
    u.push(v);
    lw.push(lnChoose(n1, v) + lnChoose(n2, s - v));
  }
  const base = lnChoose(n1, a) + lnChoose(n2, s - a);
  /** ln Σᵤ w_u e^{ηu}, y los momentos de U bajo η. */
  const momentos = (eta: number) => {
    let m = -Infinity;
    for (let j = 0; j < u.length; j++) m = Math.max(m, lw[j] + eta * u[j]);
    let z = 0, e1 = 0, e2 = 0;
    for (let j = 0; j < u.length; j++) {
      const w = Math.exp(lw[j] + eta * u[j] - m);
      z += w;
      e1 += w * u[j];
      e2 += w * u[j] * u[j];
    }
    const media = e1 / z;
    return { lnZ: m + Math.log(z), media, varianza: Math.max(0, e2 / z - media * media) };
  };
  return {
    logP: (eta) => base + eta * a - momentos(eta).lnZ,
    d1: (eta) => a - momentos(eta).media,
    d2: (eta) => -momentos(eta).varianza,
  };
}

// ---------------------------------------------------------------------------------------------
// Verosimilitud marginal con cuadratura de Gauss-Hermite adaptativa (1 dimensión)
// ---------------------------------------------------------------------------------------------

/**
 * ln ∫ P(η = mu + tau·z) φ(z) dz para un estudio. El integrando exp(g(z)), g(z) = ln P(mu + tau z) − z²/2, es log-cóncavo: se halla su
 * máximo z* por Newton con retroceso, se escala con s = 1/√(−g''(z*)) y se aplica Gauss-Hermite normalizada centrada en z*:
 * ∫ exp(g) dz = √(2π)·s·Σ wₖ exp(g(z* + s·xₖ) + xₖ²/2); el 1/√(2π) de φ cancela el factor.
 */
function logIntegralEstudio(f: Familia, mu: number, tau: number, nodos: { z: number[]; logw: number[] }): number {
  const g = (z: number) => f.logP(mu + tau * z) - 0.5 * z * z;
  let z = 0;
  let gz = g(z);
  for (let it = 0; it < 60; it++) {
    const d1 = tau * f.d1(mu + tau * z) - z;
    const d2 = tau * tau * f.d2(mu + tau * z) - 1;
    const paso = -d1 / d2;
    let t = 1;
    let zn = z + paso;
    let gn = g(zn);
    while (!(gn >= gz - 1e-14) && t > 1e-10) {
      t /= 2;
      zn = z + t * paso;
      gn = g(zn);
    }
    const mover = Math.abs(zn - z);
    z = zn;
    gz = gn;
    if (mover < 1e-12) break;
  }
  const curvatura = 1 - tau * tau * f.d2(mu + tau * z);
  const s = 1 / Math.sqrt(curvatura);
  const terminos = nodos.z.map((x, k) => nodos.logw[k] + g(z + s * x) + 0.5 * x * x);
  const m = Math.max(...terminos);
  let suma = 0;
  for (const t of terminos) suma += Math.exp(t - m);
  return Math.log(s) + m + Math.log(suma);
}

function logVerosimilitud(familias: Familia[], mu: number, tau: number, nodos: { z: number[]; logw: number[] }): number {
  let total = 0;
  for (const f of familias) total += logIntegralEstudio(f, mu, tau, nodos);
  return total;
}

// ---------------------------------------------------------------------------------------------
// Ajuste
// ---------------------------------------------------------------------------------------------

const LIM_LN_TAU: [number, number] = [Math.log(1e-3), Math.log(30)];
const NODOS = 21;

interface Ajuste {
  mu: number;
  tau: number;
  ll: number;
  se: number;
  convergio: boolean;
  limites: string[];
  advertencias: string[];
}

function ajustar(familias: Familia[], mu0: number): Ajuste | null {
  const nodos = nodosGaussHermite(NODOS);
  const nll = (x: number[]) => {
    const v = -logVerosimilitud(familias, limitar(x[0], -30, 30), Math.exp(limitar(x[1], ...LIM_LN_TAU)), nodos);
    return Number.isFinite(v) ? v : 1e12;
  };
  let mejor: { x: number[]; fx: number } | null = null;
  for (const lt0 of [Math.log(0.3), Math.log(1), Math.log(0.05)]) {
    let sol = nelderMead(nll, [mu0, lt0], { paso: 0.3, maxEval: 1500, reinicios: 3, tolX: 1e-8, tolF: 1e-11 });
    sol = nelderMead(nll, sol.x, { paso: 0.02, maxEval: 1500, reinicios: 4, tolX: 1e-11, tolF: 1e-13 });
    if (!mejor || sol.fx < mejor.fx) mejor = sol;
  }
  if (!mejor) return null;
  const reintento = nelderMead(nll, mejor.x, { paso: 0.005, maxEval: 600, reinicios: 0, tolX: 1e-11, tolF: 1e-13 });
  const estable = Math.abs(reintento.fx - mejor.fx) < 1e-8;
  if (reintento.fx < mejor.fx) mejor = reintento;
  const lnTau = limitar(mejor.x[1], ...LIM_LN_TAU);
  const enLimiteBajo = lnTau <= LIM_LN_TAU[0] + 1e-3;
  const enLimiteAlto = lnTau >= LIM_LN_TAU[1] - 1e-3;
  const enLimite = enLimiteBajo || enLimiteAlto;
  const x = [limitar(mejor.x[0], -30, 30), lnTau];
  const advertencias: string[] = [];
  // Error estándar de μ: inversa del Hessiano completo (μ, ln τ); con τ en el límite, sobre μ solo (τ fijo)
  let se: number | null = null;
  if (!enLimite) {
    const inv = invertir(hessianoNumerico(nll, x, 2e-3));
    if (inv && inv[0][0] > 0 && inv[1][1] > 0) se = Math.sqrt(inv[0][0]);
  }
  if (se === null) {
    const inv = invertir(hessianoNumerico((m) => nll([m[0], x[1]]), [x[0]], 2e-3));
    if (!inv || !(inv[0][0] > 0)) return null;
    se = Math.sqrt(inv[0][0]);
    if (enLimiteBajo) advertencias.push("La variación entre estudios (τ) quedó en el límite inferior, prácticamente nula: el error estándar se calculó con τ fijo en ese límite.");
    if (enLimiteAlto) advertencias.push("La variación entre estudios (τ) creció hasta el límite superior (30): con tan pocos datos no se puede estimar la heterogeneidad y el resultado no es fiable.");
  }
  return { mu: x[0], tau: Math.exp(x[1]), ll: -mejor.fx, se, convergio: estable, limites: enLimiteBajo ? ["tau"] : enLimiteAlto ? ["tauAlto"] : [], advertencias };
}

function envolver(
  modelo: ModeloGLMM,
  a: Ajuste,
  k: number,
  kInformativos: number,
  excluidos: { etiqueta: string; motivo: string }[],
  sinInformacion: string[],
  esOR: boolean,
): ResultadoGLMM {
  const z = normalCuantilDosColas(0.05);
  const lo = a.mu - z * a.se;
  const hi = a.mu + z * a.se;
  const nat = (v: number) => (esOR ? Math.exp(v) : sigmoide(v));
  const advertencias = [...a.advertencias];
  if (kInformativos < 5) advertencias.push(`Con solo ${kInformativos} estudios que aportan información, el intervalo de Wald de un modelo mixto puede quedarse corto y τ se estima con poca precisión.`);
  if (!a.convergio) advertencias.push("La optimización no mostró convergencia estricta (un reinicio más todavía mejoraba la verosimilitud): contrasta con otro programa antes de publicar.");
  return {
    modelo,
    k,
    kInformativos,
    excluidos,
    sinInformacion,
    mu: a.mu,
    se: a.se,
    muIc: [lo, hi],
    estimacion: nat(a.mu),
    icInferior: nat(lo),
    icSuperior: nat(hi),
    pValor: esOR ? normalDosColas(a.mu / a.se) : null,
    tau2: a.tau * a.tau,
    tau: a.tau,
    logVerosimilitud: a.ll,
    convergio: a.convergio,
    limites: a.limites,
    advertencias,
  };
}

const esEntero = (v: number) => Number.isFinite(v) && Number.isInteger(v);

/** Proporción combinada con el modelo binomial-normal (logit). `null` con menos de 2 estudios válidos. */
export function ajustarGLMMProporcion(datos: DatosProporcion[]): ResultadoGLMM | null {
  const excluidos: { etiqueta: string; motivo: string }[] = [];
  const validos: DatosProporcion[] = [];
  datos.forEach((d, i) => {
    const etiqueta = d.etiqueta ?? `Estudio ${i + 1}`;
    if (!esEntero(d.x) || !esEntero(d.n)) excluidos.push({ etiqueta, motivo: "recuentos no enteros" });
    else if (d.n <= 0 || d.x < 0 || d.x > d.n) excluidos.push({ etiqueta, motivo: "recuentos imposibles (eventos fuera de 0–N)" });
    else validos.push({ ...d, etiqueta });
  });
  if (validos.length < 2) return null;
  const familias = validos.map((d) => familiaBinomial(d.x, d.n));
  const mu0 = validos.reduce((s, d) => s + Math.log((d.x + 0.5) / (d.n - d.x + 0.5)), 0) / validos.length;
  const a = ajustar(familias, mu0);
  if (!a) return null;
  return envolver("binomial-normal", a, validos.length, validos.length, excluidos, [], false);
}

/** Odds ratio combinado con el modelo hipergeométrico-normal (verosimilitud condicional exacta). `null` con menos de 2 estudios con información. */
export function ajustarGLMMOdds(datos: DatosTabla[]): ResultadoGLMM | null {
  const excluidos: { etiqueta: string; motivo: string }[] = [];
  const validos: (DatosTabla & { etiqueta: string })[] = [];
  datos.forEach((d, i) => {
    const etiqueta = d.etiqueta ?? `Estudio ${i + 1}`;
    if (![d.eventosExp, d.totalExp, d.eventosCtrl, d.totalCtrl].every(esEntero)) excluidos.push({ etiqueta, motivo: "recuentos no enteros (¿grupo control dividido?)" });
    else if (d.totalExp <= 0 || d.totalCtrl <= 0 || d.eventosExp < 0 || d.eventosCtrl < 0 || d.eventosExp > d.totalExp || d.eventosCtrl > d.totalCtrl)
      excluidos.push({ etiqueta, motivo: "recuentos imposibles (eventos fuera de 0–N)" });
    else validos.push({ ...d, etiqueta });
  });
  const sinInformacion: string[] = [];
  const informativos: (DatosTabla & { etiqueta: string })[] = [];
  for (const d of validos) {
    const s = d.eventosExp + d.eventosCtrl;
    const soporte = Math.min(d.totalExp, s) - Math.max(0, s - d.totalCtrl);
    if (soporte <= 0) sinInformacion.push(d.etiqueta);
    else informativos.push(d);
  }
  if (informativos.length < 2) return null;
  const familias = informativos.map((d) => familiaHipergeometrica(d.eventosExp, d.totalExp, d.eventosCtrl, d.totalCtrl));
  // Punto de partida: media de los ln(OR) con +0,5 en las cuatro celdas
  const mu0 = informativos.reduce((s, d) => s + Math.log(((d.eventosExp + 0.5) * (d.totalCtrl - d.eventosCtrl + 0.5)) / ((d.totalExp - d.eventosExp + 0.5) * (d.eventosCtrl + 0.5))), 0) / informativos.length;
  const a = ajustar(familias, mu0);
  if (!a) return null;
  return envolver("hipergeometrica-normal", a, validos.length, informativos.length, excluidos, sinInformacion, true);
}

/** Log-verosimilitud marginal (con los coeficientes binomiales) de las proporciones en (μ, τ); expuesta para validar la cuadratura. */
export function logVerosimilitudProporcion(datos: DatosProporcion[], mu: number, tau: number, nodosPorEstudio = NODOS): number {
  return logVerosimilitud(datos.map((d) => familiaBinomial(d.x, d.n)), mu, tau, nodosGaussHermite(nodosPorEstudio));
}

/** Log-verosimilitud marginal de las tablas 2×2 (hipergeométrica no central) en (ln OR, τ); expuesta para validar la cuadratura. */
export function logVerosimilitudOdds(datos: DatosTabla[], mu: number, tau: number, nodosPorEstudio = NODOS): number {
  return logVerosimilitud(datos.map((d) => familiaHipergeometrica(d.eventosExp, d.totalExp, d.eventosCtrl, d.totalCtrl)), mu, tau, nodosGaussHermite(nodosPorEstudio));
}
