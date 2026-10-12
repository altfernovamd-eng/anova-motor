// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
import type { CampoSubgrupo, CertezaGrade, DatosCrudos, Estudio, Medida } from "./tipos.ts";
import { chiCuadradoCola, lnGamma, tStudentDosColas } from "./funciones-especiales.ts";
import { ajustarBivariadoConCache, type ResultadoBivariado } from "./sroc-bivariado.ts";
import { curvaMosesLittenberg, intervaloWilson, puntoDiagnostico, regresionMosesLittenberg, type PuntoDiagnostico } from "./moses-littenberg.ts";

const Z95 = 1.959963985;

/** `fisher`: z de Fisher, atanh(r), para combinar coeficientes de correlación. */
export type EscalaAnalisis = "log" | "lineal" | "logit" | "fisher";

/** Measures that can be pooled with a generic inverse-variance model, and the scale to pool them on. */
export const escalaPorMedida: Partial<Record<Medida, EscalaAnalisis>> = {
  RR: "log",
  OR: "log",
  HR: "log",
  PR: "log",
  DM: "lineal",
  DME: "lineal",
  RAR: "lineal",
  Proporción: "logit",
  RTI: "log",
  DTI: "lineal",
  RM: "log",
  r: "fisher",
};

// --- Funciones de probabilidad (χ², t de Student, normal): implementadas en `funciones-especiales.ts` ---

/** p de la cola superior de χ² con `gl` grados de libertad. */
export function chiSquarePValue(q: number, gl: number): number {
  return chiCuadradoCola(q, gl);
}

/** P(|T| > t) para una t de Student con `gl` grados de libertad: la cola de AMBOS lados (ya es el p bilateral). */
const studentTDosColas = tStudentDosColas;

/** Percentil 97,5 de la t de Student con `gl` grados de libertad (cuantil bilateral del 95 %), por bisección. */
export function studentTQuantile975(gl: number): number {
  if (gl <= 0) return NaN;
  let lo = 0;
  let hi = 1000;
  for (let i = 0; i < 200; i++) {
    const mid = (lo + hi) / 2;
    if (studentTDosColas(mid, gl) > 0.05) lo = mid;
    else hi = mid;
  }
  return (lo + hi) / 2;
}

// --- Meta-analysis pooling (generic inverse-variance, fixed-effect + DerSimonian-Laird random-effects) ---

export interface EstudioPooled {
  id: number;
  etiqueta: string;
  estimacion: number;
  icInferior: number;
  icSuperior: number;
  pesoFijo: number;
  pesoAleatorio: number;
  /** Transformed effect (log/logit/linear) and its standard error — used by the funnel plot & Egger's test. */
  y: number;
  se: number;
}

export interface ResultadoMetaAnalisis {
  medida: Medida;
  escala: EscalaAnalisis;
  estudios: EstudioPooled[];
  fijo: { estimacion: number; icInferior: number; icSuperior: number };
  aleatorio: {
    estimacion: number;
    icInferior: number;
    icSuperior: number;
    prediccionInferior: number | null;
    prediccionSuperior: number | null;
    /** Knapp-Hartung (HKSJ) adjusted 95% CI: t-based, more conservative — recommended with few studies. */
    hkIcInferior: number;
    hkIcSuperior: number;
  };
  /** Etiquetas de los estudios excluidos del OR/RR/PR por no tener eventos en ningún brazo (Handbook 10.4.4.2). */
  excluidosSinEventos: string[];
  heterogeneidad: {
    q: number;
    gl: number;
    pValor: number;
    i2: number;
    tau2: number;
    tau: number;
    metodoTau2: MetodoTau2;
  };
}

/** Lleva un valor a la escala de análisis: log (OR, RR, HR), logit (proporciones), z de Fisher (correlaciones) o el mismo valor (escala lineal). */
export function transformar(valor: number, escala: EscalaAnalisis): number {
  if (escala === "log") return Math.log(valor);
  if (escala === "logit") return Math.log(valor / (1 - valor));
  if (escala === "fisher") return Math.atanh(valor);
  return valor;
}

/** Inversa de `transformar`: devuelve un valor de la escala de análisis a su escala natural. */
export function destransformar(valor: number, escala: EscalaAnalisis): number {
  if (escala === "log") return Math.exp(valor);
  if (escala === "logit") {
    const e = Math.exp(valor);
    return e / (1 + e);
  }
  if (escala === "fisher") return Math.tanh(valor);
  return valor;
}

// --- Raw-data entry: compute the effect + 95% CI + SE from what the paper actually reports ---

export interface EfectoCalculado {
  /** Effect on its natural scale (OR/RR value, or the risk/mean difference). */
  valor: number;
  icInferior: number;
  icSuperior: number;
  /** SE on the analysis scale (log for OR/RR, linear for RAR/DM) — matches what the pooling engine uses. */
  se: number;
  /** Whether the Haldane-Anscombe 0.5 continuity correction was applied (a zero cell was present). */
  correccion: boolean;
}

/**
 * Binary 2×2 table → effect size. `medida`: "OR", "RR" or "RAR" (risk difference).
 * `eventos`/`total` are per arm. Applies the Haldane-Anscombe (+0.5) correction if any cell is zero.
 */
/**
 * Texto con el que se GUARDA un efecto calculado desde datos crudos. La estimación y los límites del IC son lo único que
 * lee el motor (de ahí reconstruye el error estándar), así que redondearlos a 4 decimales —como se muestran— hacía que el
 * resultado combinado se apartara de R hasta ~1,6 % (validación contra metafor, oct 2026). 9 cifras significativas dejan
 * el error por debajo de 1e-8 y siguen leyéndose bien.
 */
export const formatoParaGuardar = (x: number): string => String(Number(x.toPrecision(9)));

/**
 * Calcula la medida de efecto (OR, RR, RAR) a partir de una tabla 2×2.
 *
 * @param eventosTrat - Eventos en el grupo tratamiento (a)
 * @param totalTrat - Total en el grupo tratamiento (a+b)
 * @param eventosCtrl - Eventos en el grupo control (c)
 * @param totalCtrl - Total en el grupo control (c+d)
 * @param medida - "OR" (odds ratio), "RR" (riesgo relativo), "RAR" (reducción de riesgo absoluto), o
 *                 "PR" (razón de prevalencia — matemáticamente idéntica a RR, para transversales)
 * @returns Efecto calculado con IC 95% y error estándar, o null si datos inválidos
 *
 * @note Aplica corrección de continuidad (+0.5) si hay celdas vacías (Cochrane Handbook, 6.4.3.1)
 *
 * @example
 * const resultado = efectoDesde2x2(10, 50, 5, 50, "OR");
 * // {valor: 2.25, icInferior: 0.709, icSuperior: 7.141, se: 0.589, correccion: false}   (coincide con metafor::escalc)
 */
export function efectoDesde2x2(
  eventosTrat: number,
  totalTrat: number,
  eventosCtrl: number,
  totalCtrl: number,
  medida: "OR" | "RR" | "RAR" | "PR",
): EfectoCalculado | null {
  const a = eventosTrat;
  const c = eventosCtrl;
  const n1 = totalTrat;
  const n2 = totalCtrl;
  if (![a, n1, c, n2].every(Number.isFinite)) return null;
  if (n1 <= 0 || n2 <= 0) return null;
  if (a < 0 || c < 0 || a > n1 || c > n2) return null;

  let A = a;
  let B = n1 - a; // no-eventos tratamiento
  let C = c;
  let D = n2 - c; // no-eventos control
  let correccion = false;
  if (A === 0 || B === 0 || C === 0 || D === 0) {
    A += 0.5;
    B += 0.5;
    C += 0.5;
    D += 0.5;
    correccion = true;
  }
  const N1 = A + B;
  const N2 = C + D;

  const esLog = medida === "OR" || medida === "RR" || medida === "PR";
  let y: number;
  let se: number;
  if (medida === "OR") {
    y = Math.log((A * D) / (B * C));
    se = Math.sqrt(1 / A + 1 / B + 1 / C + 1 / D);
  } else if (medida === "RR" || medida === "PR") {
    // PR (razón de prevalencia) se calcula idéntico a RR — misma estructura 2×2, distinto nombre semántico.
    y = Math.log((A / N1) / (C / N2));
    se = Math.sqrt(1 / A - 1 / N1 + 1 / C - 1 / N2);
  } else {
    y = A / N1 - C / N2; // diferencia de riesgos (escala lineal)
    se = Math.sqrt((A * B) / N1 ** 3 + (C * D) / N2 ** 3);
  }
  if (!Number.isFinite(y) || !Number.isFinite(se) || se <= 0) return null;

  return {
    valor: esLog ? Math.exp(y) : y,
    icInferior: esLog ? Math.exp(y - Z95 * se) : y - Z95 * se,
    icSuperior: esLog ? Math.exp(y + Z95 * se) : y + Z95 * se,
    se,
    correccion,
  };
}

/**
 * Calcula una proporción de un solo brazo (p.ej. prevalencia) a partir de eventos/total — sin grupo
 * de comparación. Transforma en escala logit (igual que `escalaPorMedida.Proporción`) para el error
 * estándar, y aplica corrección de continuidad cuando el resultado es 0% o 100% (donde logit(0) y
 * logit(1) no están definidos — sin esto, esos estudios, perfectamente legítimos en la práctica,
 * quedarían excluidos del meta-análisis o forzarían a calcular el IC a mano).
 *
 * @param eventos - Nº de casos observados
 * @param total - Tamaño de la muestra
 * @returns Efecto calculado (proporción + IC 95% en escala natural 0–1), o null si datos inválidos
 *
 * @note Misma filosofía que `efectoDesde2x2`: con 0 o 100%, se suma 0.5 a eventos y 1 al total antes
 *       de calcular tanto el punto como el IC — consistente con la corrección de Haldane-Anscombe ya
 *       usada para tablas 2×2 con casillas en cero.
 *
 * @example
 * const resultado = efectoDesdeProporcion(0, 20); // 0 de 20 — caso límite
 * // {valor: ~0.024, icInferior: ~0.001, icSuperior: ~0.15, se: ..., correccion: true}
 */
export function efectoDesdeProporcion(eventos: number, total: number): EfectoCalculado | null {
  if (![eventos, total].every(Number.isFinite)) return null;
  if (total <= 0 || eventos < 0 || eventos > total) return null;

  let e = eventos;
  let n = total;
  let correccion = false;
  if (e === 0 || e === n) {
    e += 0.5;
    n += 1;
    correccion = true;
  }
  const p = e / n;
  const y = Math.log(p / (1 - p)); // logit
  const se = Math.sqrt(1 / e + 1 / (n - e));
  if (!Number.isFinite(y) || !Number.isFinite(se) || se <= 0) return null;

  const destransformarLogit = (t: number) => {
    const ex = Math.exp(t);
    return ex / (1 + ex);
  };

  return {
    valor: p,
    icInferior: destransformarLogit(y - Z95 * se),
    icSuperior: destransformarLogit(y + Z95 * se),
    se,
    correccion,
  };
}

/** Continuous outcome (mean, SD, N per arm) → mean difference (DM) + 95% CI + SE. */
/**
 * Calcula la diferencia de medias (DM) a partir de medias, DE y n por grupo.
 *
 * @param mediaTrat - Media del grupo tratamiento
 * @param deTrat - Desviación estándar del grupo tratamiento
 * @param nTrat - Tamaño de muestra del grupo tratamiento
 * @param mediaCtrl - Media del grupo control
 * @param deCtrl - Desviación estándar del grupo control
 * @param nCtrl - Tamaño de muestra del grupo control
 * @returns DM con IC 95% y error estándar, o null si datos inválidos
 *
 * @note Escala original (no estandarizada). Usado para desenlaces continuos medidos en unidades iguales.
 */
export function efectoDesdeMedias(
  mediaTrat: number,
  deTrat: number,
  nTrat: number,
  mediaCtrl: number,
  deCtrl: number,
  nCtrl: number,
): EfectoCalculado | null {
  if (![mediaTrat, deTrat, nTrat, mediaCtrl, deCtrl, nCtrl].every(Number.isFinite)) return null;
  if (nTrat < 1 || nCtrl < 1 || deTrat < 0 || deCtrl < 0) return null;
  const y = mediaTrat - mediaCtrl;
  const se = Math.sqrt((deTrat * deTrat) / nTrat + (deCtrl * deCtrl) / nCtrl);
  if (!Number.isFinite(se) || se <= 0) return null;
  return { valor: y, icInferior: y - Z95 * se, icSuperior: y + Z95 * se, se, correccion: false };
}

/**
 * Medias (2 grupos) → **diferencia de medias estandarizada (g de Hedges)**. Réplica de
 * `metafor::escalc(measure="SMD")`: `d = (m̄₁−m̄₂)/sₚ`, corrección de sesgo de muestra pequeña
 * `J = Γ(m/2) / (√(m/2)·Γ((m−1)/2))` con `m = n₁+n₂−2`, `g = J·d`, y varianza (Hedges 1982, ec. 8):
 * `v = 1/n₁ + 1/n₂ + g²/(2(n₁+n₂))`. Útil cuando los estudios miden el mismo constructo con escalas
 * distintas. `null` si los datos no son válidos.
 */
/**
 * Calcula la g de Hedges (diferencia de medias estandarizada, DME) a partir de medias, DE y n.
 *
 * @param mediaTrat - Media del grupo tratamiento
 * @param deTrat - Desviación estándar del grupo tratamiento
 * @param nTrat - Tamaño de muestra del grupo tratamiento
 * @param mediaCtrl - Media del grupo control
 * @param deCtrl - Desviación estándar del grupo control
 * @param nCtrl - Tamaño de muestra del grupo control
 * @returns g de Hedges con IC 95%, o null si datos inválidos
 *
 * @note DME = (Media1 - Media2) / DE agrupada. Sin unidades, comparable entre escalas distintas.
 *        Interpretación: 0.2=pequeño, 0.5=medio, 0.8=grande (Cohen's d).
 *        Fórmula: Hedges' g (versión insesgada, con corrección J).
 */
export function efectoDesdeMediasEstandarizada(
  mediaTrat: number,
  deTrat: number,
  nTrat: number,
  mediaCtrl: number,
  deCtrl: number,
  nCtrl: number,
): EfectoCalculado | null {
  if (![mediaTrat, deTrat, nTrat, mediaCtrl, deCtrl, nCtrl].every(Number.isFinite)) return null;
  if (nTrat < 2 || nCtrl < 2 || deTrat < 0 || deCtrl < 0) return null;
  const m = nTrat + nCtrl - 2;
  if (m <= 1) return null;
  const sdp = Math.sqrt(((nTrat - 1) * deTrat * deTrat + (nCtrl - 1) * deCtrl * deCtrl) / m);
  if (!(sdp > 0)) return null;
  const d = (mediaTrat - mediaCtrl) / sdp;
  const J = Math.exp(lnGamma(m / 2) - Math.log(Math.sqrt(m / 2)) - lnGamma((m - 1) / 2));
  const g = J * d;
  const n = nTrat + nCtrl;
  const se = Math.sqrt(1 / nTrat + 1 / nCtrl + (g * g) / (2 * n));
  if (!Number.isFinite(se) || se <= 0) return null;
  return { valor: g, icInferior: g - Z95 * se, icSuperior: g + Z95 * se, se, correccion: false };
}

/**
 * Tasas de incidencia (eventos y tiempo-persona de cada brazo) → **razón de tasas** (`RTI`) o **diferencia de tasas** (`DTI`)
 * (Cochrane Handbook 6.7, 10.8). Réplica de `metafor::escalc(measure = "IRR" / "IRD", add = 1/2, to = "only0")`:
 *  · RTI: ln(RTI) = ln[(e₁/t₁)/(e₀/t₀)], con error estándar √(1/e₁ + 1/e₀) del logaritmo (eventos de Poisson).
 *  · DTI: e₁/t₁ − e₀/t₀, con error estándar √(e₁/t₁² + e₀/t₀²).
 * Si algún brazo tiene cero eventos se suma 0,5 a los dos conteos (no al tiempo), igual que en la tabla 2×2. Un estudio sin
 * eventos en ningún brazo no informa la razón y queda fuera del combinado (ver `excluidoPorSinEventos`).
 * `null` si los datos no son válidos (tiempos ≤ 0, eventos negativos o no finitos).
 */
export function efectoDesdeTasas(
  eventosTrat: number,
  tiempoTrat: number,
  eventosCtrl: number,
  tiempoCtrl: number,
  medida: "RTI" | "DTI",
): EfectoCalculado | null {
  if (![eventosTrat, tiempoTrat, eventosCtrl, tiempoCtrl].every(Number.isFinite)) return null;
  if (tiempoTrat <= 0 || tiempoCtrl <= 0 || eventosTrat < 0 || eventosCtrl < 0) return null;
  let e1 = eventosTrat;
  let e0 = eventosCtrl;
  const correccion = e1 === 0 || e0 === 0;
  if (correccion) {
    e1 += 0.5;
    e0 += 0.5;
  }
  let y: number;
  let se: number;
  if (medida === "RTI") {
    y = Math.log(e1 / tiempoTrat / (e0 / tiempoCtrl));
    se = Math.sqrt(1 / e1 + 1 / e0);
  } else {
    y = e1 / tiempoTrat - e0 / tiempoCtrl;
    se = Math.sqrt(e1 / tiempoTrat ** 2 + e0 / tiempoCtrl ** 2);
  }
  if (!Number.isFinite(y) || !Number.isFinite(se) || se <= 0) return null;
  const esLog = medida === "RTI";
  return {
    valor: esLog ? Math.exp(y) : y,
    icInferior: esLog ? Math.exp(y - Z95 * se) : y - Z95 * se,
    icSuperior: esLog ? Math.exp(y + Z95 * se) : y + Z95 * se,
    se,
    correccion,
  };
}

/**
 * Medias (2 grupos) → **razón de medias** (Cochrane Handbook 6.5.1.3; Friedrich, Adhikari y Beyene 2008): ln(RM) = ln(m₁/m₂), con
 * varianza (de primer orden) sd₁²/(n₁·m₁²) + sd₂²/(n₂·m₂²). Solo tiene sentido con medidas que no pueden ser negativas (las dos
 * medias deben ser positivas). No incluye la corrección de sesgo de muestra pequeña de `metafor::escalc("ROM")` (Hedges y cols.
 * 1999), que cambia el logaritmo en una fracción de 1/n: la del Handbook es la fórmula simple. `null` si no son válidos.
 */
export function efectoDesdeRazonDeMedias(
  mediaTrat: number,
  deTrat: number,
  nTrat: number,
  mediaCtrl: number,
  deCtrl: number,
  nCtrl: number,
): EfectoCalculado | null {
  if (![mediaTrat, deTrat, nTrat, mediaCtrl, deCtrl, nCtrl].every(Number.isFinite)) return null;
  if (mediaTrat <= 0 || mediaCtrl <= 0 || nTrat < 1 || nCtrl < 1 || deTrat < 0 || deCtrl < 0) return null;
  const y = Math.log(mediaTrat / mediaCtrl);
  const se = Math.sqrt((deTrat * deTrat) / (nTrat * mediaTrat * mediaTrat) + (deCtrl * deCtrl) / (nCtrl * mediaCtrl * mediaCtrl));
  if (!Number.isFinite(y) || !Number.isFinite(se) || se <= 0) return null;
  return { valor: Math.exp(y), icInferior: Math.exp(y - Z95 * se), icSuperior: Math.exp(y + Z95 * se), se, correccion: false };
}

/**
 * Coeficiente de correlación r de un estudio → su efecto en la **escala z de Fisher**: z = atanh(r), con error estándar
 * 1/√(n − 3) (`metafor::escalc("ZCOR")`). El IC se calcula en z y se devuelve a la escala de r con tanh. Necesita n > 3 y
 * −1 < r < 1. `se` está en la escala z, como el de las demás medidas en su escala de análisis.
 */
export function efectoDesdeCorrelacion(r: number, n: number): EfectoCalculado | null {
  if (![r, n].every(Number.isFinite) || n <= 3 || r <= -1 || r >= 1) return null;
  const z = Math.atanh(r);
  const se = 1 / Math.sqrt(n - 3);
  return { valor: r, icInferior: Math.tanh(z - Z95 * se), icSuperior: Math.tanh(z + Z95 * se), se, correccion: false };
}

/**
 * Tiempo hasta el evento: **O−E y V** (observados menos esperados y su varianza, como los da una prueba de rangos logarítmicos)
 * → razón de riesgos instantáneos HR (Cochrane Handbook 6.8.2 y 10.9; Peto): ln(HR) = (O−E)/V, con error estándar 1/√V.
 * Combinar estos ln(HR) con la varianza inversa en un modelo fijo equivale exactamente al método de Peto, Σ(O−E)/ΣV. `null` si V ≤ 0.
 */
export function efectoDesdeOEV(oMenosE: number, varianza: number): EfectoCalculado | null {
  if (![oMenosE, varianza].every(Number.isFinite) || varianza <= 0) return null;
  const y = oMenosE / varianza;
  const se = 1 / Math.sqrt(varianza);
  return { valor: Math.exp(y), icInferior: Math.exp(y - Z95 * se), icSuperior: Math.exp(y + Z95 * se), se, correccion: false };
}

/**
 * Ensayo **cruzado (crossover)**, desenlace continuo → diferencia de medias PAREADA (Cochrane
 * Handbook cap. 23.2.7.1). A diferencia de `efectoDesdeMedias` (dos grupos INDEPENDIENTES), un
 * cruzado mide al MISMO grupo de `n` participantes bajo ambas condiciones — usar la fórmula de grupos
 * independientes aquí ignora la correlación entre periodos y calcula un error estándar incorrecto
 * (el "error de unidad de análisis" del cap. 23), sea de más o de menos según el signo de la
 * correlación real.
 *
 * `SE(DM) = SDdiff / √n`. Si el estudio no reporta `SDdiff` directamente, se imputa desde las DE de
 * cada periodo y una correlación (reportada por el estudio, o tomada de otro estudio pareado de la
 * revisión con análisis de sensibilidad — el Handbook exige declarar esto, no aNova inventarlo):
 * `SDdiff = √(SDE² + SDC² − 2·Corr·SDE·SDC)` (verbatim del Handbook).
 *
 * @param sdDiff - DE de las diferencias intra-paciente, si el estudio la reporta directamente (NaN si no)
 * @param sdTrat / sdCtrl / correlacion - alternativa cuando `sdDiff` no se reporta
 * @returns `null` si no hay suficiente información por NINGUNA de las dos vías — nunca inventa un
 *          valor de correlación ni de SDdiff.
 */
export function efectoDesdeCruzado(
  mediaTrat: number,
  mediaCtrl: number,
  n: number,
  sdDiff: number,
  sdTrat: number,
  sdCtrl: number,
  correlacion: number,
): EfectoCalculado | null {
  if (![mediaTrat, mediaCtrl, n].every(Number.isFinite)) return null;
  if (n < 2) return null;

  let sd = sdDiff;
  if (!Number.isFinite(sd) || sd <= 0) {
    if (![sdTrat, sdCtrl, correlacion].every(Number.isFinite)) return null;
    if (sdTrat < 0 || sdCtrl < 0 || correlacion < -1 || correlacion > 1) return null;
    const varianza = sdTrat * sdTrat + sdCtrl * sdCtrl - 2 * correlacion * sdTrat * sdCtrl;
    if (!(varianza > 0)) return null;
    sd = Math.sqrt(varianza);
  }

  const y = mediaTrat - mediaCtrl;
  const se = sd / Math.sqrt(n);
  if (!Number.isFinite(se) || se <= 0) return null;
  return { valor: y, icInferior: y - Z95 * se, icSuperior: y + Z95 * se, se, correccion: false };
}

/**
 * **Efecto de diseño** de un ensayo clúster-aleatorizado (Cochrane Handbook cap. 23.1): infla el
 * error estándar por `√(efecto de diseño)` en vez de dividir N/eventos por él ("effective sample
 * size" method) — el propio Handbook los llama equivalentes (cap. 23.1.5, "a slightly more flexible
 * approach"), y este evita redondear N a un entero, funcionando con la MISMA fórmula tanto para
 * desenlaces binarios como continuos (a diferencia del método de N efectiva, que divide también los
 * eventos SOLO para binarios).
 *
 * `Efecto de diseño = 1 + (M−1)·ICC`, con `M` = tamaño medio del clúster = `n / número de clústeres`.
 *
 * Ignorar este ajuste hace que el IC quede "artificialmente estrecho" y el estudio reciba más peso
 * del que le corresponde en el combinado — el Handbook lo llama la causa más común de falsos
 * positivos en revisiones que incluyen ensayos clúster sin ajustar (error de unidad de análisis).
 *
 * @returns El efecto de diseño (≥1), o `null` si `n`/`numClusters`/`icc` no son válidos (ICC fuera de 0–1, M < 1…) — nunca se
 *          aplica un ajuste con datos incompletos; el estudio debe excluirse del combinado hasta que
 *          se complete (igual criterio que el resto de datos inválidos en `extraerEstudiosValidos`).
 */
export function efectoDeDisenoCluster(n: number, numClusters: number, icc: number): number | null {
  if (![n, numClusters, icc].every(Number.isFinite)) return null;
  // El ICC es una correlación: va de 0 a 1 (un valor mayor no es un dato, es un error de captura).
  if (n <= 0 || numClusters <= 0 || icc < 0 || icc > 1) return null;
  const m = n / numClusters;
  if (m < 1) return null;
  const de = 1 + (m - 1) * icc;
  return de >= 1 ? de : null;
}

/** Medidas de razón que se calculan desde una tabla 2×2 y no pueden estimarse sin algún evento que difiera. */
const RAZONES_DE_TABLA_2X2: ReadonlySet<Medida> = new Set<Medida>(["OR", "RR", "PR"]);

/**
 * Tabla 2×2 que no informa NINGÚN efecto relativo: sin eventos en los dos brazos, o con el evento en todos los
 * participantes de los dos brazos (sin «no eventos»). Es el criterio de `drop00` de metafor y el que ya usa
 * Mantel-Haenszel en aNova. `false` si no hay tabla 2×2 completa y numérica.
 */
export function esTablaSinInformacion(crudos: DatosCrudos | undefined): boolean {
  if (!crudos || crudos.tipo !== "2x2") return false;
  const a = parseFloat(crudos.eventosTrat);
  const n1 = parseFloat(crudos.totalTrat);
  const c = parseFloat(crudos.eventosCtrl);
  const n2 = parseFloat(crudos.totalCtrl);
  if (![a, n1, c, n2].every(Number.isFinite)) return false;
  return (a === 0 && c === 0) || (a === n1 && c === n2);
}

/**
 * ¿Queda este estudio fuera del OR/RR/PR por no tener eventos en ningún brazo? Cochrane Handbook 10.4.4.2: «Studies with
 * no events in either arm … do not provide any information about the direction of the effect and are excluded from
 * odds ratio or risk ratio analyses». No aplica a la diferencia de riesgos, donde sí aportan (RD = 0). Solo se detecta
 * con datos crudos 2×2; un estudio capturado solo con su efecto e IC no se puede reconocer.
 */
export function excluidoPorSinEventos(e: Estudio): boolean {
  if (RAZONES_DE_TABLA_2X2.has(e.medida)) return esTablaSinInformacion(e.crudos);
  // Razón de tasas: sin eventos en ningún brazo tampoco informa la dirección del efecto (mismo criterio que el OR/RR).
  if (e.medida === "RTI" && e.crudos?.tipo === "tasas") {
    return parseFloat(e.crudos.eventosTrat) === 0 && parseFloat(e.crudos.eventosCtrl) === 0;
  }
  return false;
}

/**
 * Divisor válido del grupo control compartido (Handbook 23.3.4): un entero ≥ 2, o 1 (sin dividir) si no hay división
 * indicada o el valor no es válido. Lo usan el cálculo del efecto y todas las pruebas que releen los datos crudos
 * (Mantel-Haenszel, Peto, Peters, riesgo basal) para que el mismo control no cuente dos veces en ninguna.
 */
export function divisorControl(dividirControlEn: number | undefined): number {
  return Number.isInteger(dividirControlEn) && (dividirControlEn as number) >= 2 ? (dividirControlEn as number) : 1;
}

/** Whether a raw (untransformed) value is inside the domain the scale can transform. */
function dentroDeRango(valor: number, escala: EscalaAnalisis): boolean {
  if (escala === "log") return valor > 0;
  if (escala === "logit") return valor > 0 && valor < 1;
  if (escala === "fisher") return valor > -1 && valor < 1;
  return true;
}

interface EstudioValido {
  id: number;
  etiqueta: string;
  estimacion: number;
  icInferior: number;
  icSuperior: number;
  y: number;
  se: number;
}

/**
 * Parses, validates and transforms the studies of one measure into `{y (transformed effect), se}`
 * pairs — the single source of truth used by pooling, sensitivity, Egger's test and the funnel plot.
 *
 * Ensayos clúster-aleatorizados (`varianteEnsayo.tipo === "cluster"`, Handbook cap. 23.1): el SE se
 * infla por `√(efecto de diseño)` ANTES de entrar al pool, así todo lo que consume `EstudioValido`
 * (peso, heterogeneidad, Egger, GRADE) ya ve la precisión correcta sin tener que saber de clústeres.
 * El IC natural también se recalcula desde el `y`/`se` ajustados, para que el forest plot muestre el
 * mismo intervalo que realmente entra al cálculo — nunca uno más estrecho que el usado internamente.
 * Los estudios con tabla 2×2 sin eventos en ningún brazo se excluyen del OR/RR/PR (`excluidoPorSinEventos`).
 * Si faltan los datos del clúster (`numClusters`/`icc`) o son inválidos, el estudio se EXCLUYE del
 * combinado (mismo criterio que cualquier otro dato inválido aquí) en vez de incluirlo sin ajustar,
 * que sería aplicar silenciosamente el error de unidad de análisis que este ajuste existe para evitar.
 */
function extraerEstudiosValidos(estudios: Estudio[], medida: Medida): EstudioValido[] {
  const escala = escalaPorMedida[medida];
  if (!escala) return [];

  const validos: EstudioValido[] = [];
  estudios
    .filter((e) => e.medida === medida)
    .forEach((e) => {
      const est = parseFloat(e.valor);
      const inf = parseFloat(e.icInferior);
      const sup = parseFloat(e.icSuperior);
      if (!Number.isFinite(est) || !Number.isFinite(inf) || !Number.isFinite(sup)) return;
      if (inf >= sup) return;
      if (!dentroDeRango(est, escala) || !dentroDeRango(inf, escala) || !dentroDeRango(sup, escala)) return;
      if (excluidoPorSinEventos(e)) return; // sin eventos en ningún brazo: no informa el OR/RR (Handbook 10.4.4.2)

      const y = transformar(est, escala);
      const yInf = transformar(inf, escala);
      const ySup = transformar(sup, escala);
      let se = (ySup - yInf) / (2 * Z95);
      if (!Number.isFinite(se) || se <= 0) return;

      let icInferiorFinal = inf;
      let icSuperiorFinal = sup;
      if (e.varianteEnsayo?.tipo === "cluster") {
        const de = efectoDeDisenoCluster(
          parseFloat(e.n),
          parseFloat(e.varianteEnsayo.numClusters),
          parseFloat(e.varianteEnsayo.icc),
        );
        if (de === null) return; // datos de clúster incompletos/inválidos: excluir, no adivinar.
        se = se * Math.sqrt(de);
        icInferiorFinal = destransformar(y - Z95 * se, escala);
        icSuperiorFinal = destransformar(y + Z95 * se, escala);
      }

      validos.push({
        id: e.id,
        etiqueta: `${e.autor || "Estudio"}${e.año ? ` (${e.año})` : ""}`,
        estimacion: est,
        icInferior: icInferiorFinal,
        icSuperior: icSuperiorFinal,
        y,
        se,
      });
    });
  return validos;
}

/** Estimador de la varianza entre estudios τ². DL = DerSimonian-Laird (cerrado), PM = Paule-Mandel, REML. */
export type MetodoTau2 = "DL" | "PM" | "REML";

export const ETIQUETA_TAU2: Record<MetodoTau2, string> = {
  DL: "DerSimonian-Laird",
  PM: "Paule-Mandel",
  REML: "Máxima verosimilitud restringida (REML)",
};

/** Un efecto en escala transformada y su error estándar — lo mínimo que necesitan el pooling y los estimadores de τ². */
type PuntoEfecto = { y: number; se: number };

/** Media ponderada y sumas de pesos aleatorios para un τ² dado (usado por PM y REML). */
function ponderarParaTau2(validos: readonly PuntoEfecto[], t2: number): { mu: number; sw: number; sw2: number } {
  let sw = 0;
  let swy = 0;
  let sw2 = 0;
  for (const v of validos) {
    const w = 1 / (v.se * v.se + t2);
    sw += w;
    swy += w * v.y;
    sw2 += w * w;
  }
  return { mu: swy / sw, sw, sw2 };
}

/** Q generalizado Σ (yᵢ−ȳ*)²/(vᵢ+τ²); Paule-Mandel busca el τ² donde vale k−1. Monótono decreciente en τ². */
function qGeneralizado(validos: readonly PuntoEfecto[], t2: number): number {
  const { mu } = ponderarParaTau2(validos, t2);
  let s = 0;
  for (const v of validos) s += (v.y - mu) ** 2 / (v.se * v.se + t2);
  return s;
}

/**
 * Estima τ² por el método elegido. DL es cerrado (lo recibe ya calculado). PM resuelve Q generalizado = k−1
 * por bisección (Paule & Mandel 1982; ≡ empirical Bayes). REML itera el punto fijo
 * τ²=Σwᵢ²[(yᵢ−ȳ)²−vᵢ]/Σwᵢ² + 1/Σwᵢ (DerSimonian & Kacker 2007; Viechtbauer 2005). Ambos se truncan en 0.
 */
function estimarTau2(validos: readonly PuntoEfecto[], metodo: MetodoTau2, tau2DL: number): number {
  if (metodo === "DL") return tau2DL;
  const gl = validos.length - 1;

  if (metodo === "PM") {
    if (qGeneralizado(validos, 0) <= gl) return 0; // Q ≤ k−1: sin heterogeneidad detectable
    let lo = 0;
    let hi = 1;
    let guard = 0;
    while (qGeneralizado(validos, hi) > gl && guard++ < 200) hi *= 2;
    for (let it = 0; it < 200; it++) {
      const mid = (lo + hi) / 2;
      if (qGeneralizado(validos, mid) > gl) lo = mid;
      else hi = mid;
    }
    return (lo + hi) / 2;
  }

  // REML: punto fijo, arrancando desde DL. La tolerancia es RELATIVA a la escala de las varianzas de los estudios: con efectos de
  // magnitud pequeña (p. ej. una diferencia de tasas, del orden de 0,001) τ² es del orden de 1e-8 y una tolerancia absoluta de
  // 1e-12 detenía la iteración antes de converger.
  const escalaVarianza = validos.reduce((s, v) => s + v.se * v.se, 0) / validos.length;
  let t2 = tau2DL;
  for (let it = 0; it < 20000; it++) {
    const { mu, sw, sw2 } = ponderarParaTau2(validos, t2);
    let numer = 0;
    for (const v of validos) {
      const w = 1 / (v.se * v.se + t2);
      numer += w * w * ((v.y - mu) ** 2 - v.se * v.se);
    }
    let t2new = numer / sw2 + 1 / sw;
    if (t2new < 0) t2new = 0;
    if (Math.abs(t2new - t2) <= 1e-14 * Math.max(escalaVarianza, t2new)) {
      t2 = t2new;
      break;
    }
    t2 = t2new;
  }
  return t2;
}

// --- IC de τ²/I²/τ por el método Q-profile (Viechtbauer 2007) ---

/** Cuantil de la χ² con `gl` grados de libertad para probabilidad acumulada `pAcum`, por bisección sobre la cola. */
function chiCuadradoCuantil(pAcum: number, gl: number): number {
  const colaObjetivo = 1 - pAcum; // chiSquarePValue devuelve la cola superior P(χ² > q)
  let lo = 0;
  let hi = Math.max(10, gl * 10);
  let guard = 0;
  while (chiSquarePValue(hi, gl) > colaObjetivo && guard++ < 200) hi *= 2;
  for (let it = 0; it < 200; it++) {
    const mid = (lo + hi) / 2;
    if (chiSquarePValue(mid, gl) > colaObjetivo) lo = mid;
    else hi = mid;
  }
  return (lo + hi) / 2;
}

/** τ² ≥ 0 tal que Q_gen(τ²) = objetivo (Q_gen es decreciente en τ²); 0 si ni siquiera τ²=0 lo alcanza. */
function raizQPerfil(validos: readonly PuntoEfecto[], objetivo: number): number {
  if (qGeneralizado(validos, 0) <= objetivo) return 0;
  let lo = 0;
  let hi = 1;
  let guard = 0;
  while (qGeneralizado(validos, hi) > objetivo && guard++ < 200) hi *= 2;
  for (let it = 0; it < 200; it++) {
    const mid = (lo + hi) / 2;
    if (qGeneralizado(validos, mid) > objetivo) lo = mid;
    else hi = mid;
  }
  return (lo + hi) / 2;
}

export interface ICHeterogeneidad {
  tau2Inferior: number;
  tau2Superior: number;
  tauInferior: number;
  tauSuperior: number;
  i2Inferior: number;
  i2Superior: number;
}

/**
 * **Intervalos de confianza del 95% para τ², τ e I² por el método Q-profile** (Viechtbauer 2007). Invierte el
 * estadístico Q generalizado `Q_gen(τ²)=Σ(yᵢ−ȳ)²/(vᵢ+τ²)` (que sigue χ²_{k−1}): los límites de τ² son los τ²
 * que igualan Q_gen a los cuantiles χ² 0.975 y 0.025. I² se deriva con la varianza intra-estudio típica
 * `s²=(k−1)Σw/((Σw)²−Σw²)`. Réplica de `metafor::confint` / `meta` (método por defecto). `null` si <2 estudios.
 */
export function confianzaHeterogeneidad(resultado: ResultadoMetaAnalisis): ICHeterogeneidad | null {
  const pts: PuntoEfecto[] = resultado.estudios.map((e) => ({ y: e.y, se: e.se }));
  const k = pts.length;
  if (k < 2) return null;
  const gl = k - 1;
  const { sw, sw2 } = ponderarParaTau2(pts, 0); // pesos fijos w=1/v
  const s2 = (gl * sw) / (sw * sw - sw2); // varianza intra-estudio "típica" (Higgins & Thompson)

  const tau2Inferior = raizQPerfil(pts, chiCuadradoCuantil(0.975, gl)); // objetivo alto → τ² bajo
  const tau2Superior = raizQPerfil(pts, chiCuadradoCuantil(0.025, gl)); // objetivo bajo → τ² alto
  const i2 = (t2: number) => (100 * t2) / (t2 + s2);

  return {
    tau2Inferior,
    tau2Superior,
    tauInferior: Math.sqrt(tau2Inferior),
    tauSuperior: Math.sqrt(tau2Superior),
    i2Inferior: i2(tau2Inferior),
    i2Superior: i2(tau2Superior),
  };
}

// --- Meta-regresión general (moderadores múltiples y categóricos), efectos mixtos con τ² por REML ---

/** Inversa de una matriz p×p por Gauss-Jordan; `null` si es singular. */
function invertirMatriz(m: number[][]): number[][] | null {
  const n = m.length;
  const a = m.map((fila, i) => [...fila, ...Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))]);
  for (let col = 0; col < n; col++) {
    let piv = col;
    for (let r = col + 1; r < n; r++) if (Math.abs(a[r][col]) > Math.abs(a[piv][col])) piv = r;
    if (Math.abs(a[piv][col]) < 1e-12) return null;
    [a[col], a[piv]] = [a[piv], a[col]];
    const d = a[col][col];
    for (let j = 0; j < 2 * n; j++) a[col][j] /= d;
    for (let r = 0; r < n; r++) {
      if (r === col) continue;
      const f = a[r][col];
      for (let j = 0; j < 2 * n; j++) a[r][j] -= f * a[col][j];
    }
  }
  return a.map((fila) => fila.slice(n));
}

/** WLS con pesos wᵢ=1/(vᵢ+τ²): devuelve la inversa de X'WX, los coeficientes β y los residuos. `null` si singular. */
function ajusteWLS(
  y: number[],
  v: number[],
  X: number[][],
  t2: number,
): { A: number[][]; beta: number[]; e: number[]; w: number[] } | null {
  const n = y.length;
  const p = X[0].length;
  const w = v.map((vi) => 1 / (vi + t2));
  const C = Array.from({ length: p }, () => new Array(p).fill(0));
  const XtWy = new Array(p).fill(0);
  for (let i = 0; i < n; i++) {
    for (let a = 0; a < p; a++) {
      XtWy[a] += w[i] * X[i][a] * y[i];
      for (let b = a; b < p; b++) C[a][b] += w[i] * X[i][a] * X[i][b];
    }
  }
  for (let a = 0; a < p; a++) for (let b = 0; b < a; b++) C[a][b] = C[b][a];
  const A = invertirMatriz(C);
  if (!A) return null;
  const beta = A.map((fila) => fila.reduce((s, val, j) => s + val * XtWy[j], 0));
  const e = y.map((yi, i) => yi - X[i].reduce((s, xij, j) => s + xij * beta[j], 0));
  return { A, beta, e, w };
}

/** τ² residual por REML (punto fijo con la matriz proyección P), para la meta-regresión de efectos mixtos. */
function remlTau2Regresion(y: number[], v: number[], X: number[][]): number {
  const n = y.length;
  const p = X[0].length;
  let t2 = 0.05;
  for (let it = 0; it < 500; it++) {
    const aj = ajusteWLS(y, v, X, t2);
    if (!aj) break;
    const { A, e, w } = aj;
    // Py = W e ; y'PPy = Σ (wᵢeᵢ)²
    let yPPy = 0;
    for (let i = 0; i < n; i++) yPPy += (w[i] * e[i]) ** 2;
    // P = diag(w) − M A Mᵀ  con M = W X ; trazas de P y de P²
    const M = X.map((fila, i) => fila.map((xij) => w[i] * xij));
    const P = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        let s = 0;
        for (let a = 0; a < p; a++) for (let b = 0; b < p; b++) s += M[i][a] * A[a][b] * M[j][b];
        P[i][j] = (i === j ? w[i] : 0) - s;
      }
    }
    let trP = 0;
    let trPP = 0;
    for (let i = 0; i < n; i++) {
      trP += P[i][i];
      for (let j = 0; j < n; j++) trPP += P[i][j] * P[i][j];
    }
    if (trPP === 0) break;
    let t2n = t2 + (yPPy - trP) / trPP;
    if (t2n < 0) t2n = 0;
    if (Math.abs(t2n - t2) < 1e-12) {
      t2 = t2n;
      break;
    }
    t2 = t2n;
  }
  return t2;
}

export interface CoeficienteRegresion {
  nombre: string;
  estimacion: number;
  se: number;
  z: number;
  pValor: number;
  icInferior: number;
  icSuperior: number;
}

export interface ResultadoMetaRegresionMulti {
  /** Coeficientes (el primero es el intercepto). */
  coeficientes: CoeficienteRegresion[];
  tau2Residual: number;
  /** Test ómnibus de los moderadores (todos menos el intercepto): estadístico Wald QM ~ χ²_gl. */
  qm: number;
  qmGl: number;
  qmPValor: number;
  /** % de la heterogeneidad explicada por los moderadores. */
  r2: number;
  k: number;
}

/**
 * **Meta-regresión de efectos mixtos** con una matriz de diseño arbitraria `X` (primera columna = intercepto),
 * τ² residual por **REML** (réplica de `metafor::rma(mods=…, method="REML")`). Soporta moderadores múltiples y
 * categóricos (codificados como columnas dummy por quien construye `X`). Devuelve coeficientes con IC, el test
 * ómnibus QM y R². `null` si hay menos estudios que parámetros o el diseño es singular.
 */
export function ajustarMetaRegresion(
  puntos: readonly PuntoEfecto[],
  X: number[][],
  nombres: string[],
): ResultadoMetaRegresionMulti | null {
  const y = puntos.map((p) => p.y);
  const v = puntos.map((p) => p.se * p.se);
  const n = y.length;
  const p = X[0].length;
  if (n <= p) return null;

  const t2 = remlTau2Regresion(y, v, X);
  const aj = ajusteWLS(y, v, X, t2);
  if (!aj) return null;
  const { A, beta } = aj;

  const coeficientes: CoeficienteRegresion[] = beta.map((b, j) => {
    const se = Math.sqrt(Math.max(0, A[j][j]));
    const z = se > 0 ? b / se : 0;
    return {
      nombre: nombres[j] ?? `x${j}`,
      estimacion: b,
      se,
      z,
      pValor: se > 0 ? Math.min(1, colaNormalDosLados(z)) : 1,
      icInferior: b - Z95 * se,
      icSuperior: b + Z95 * se,
    };
  });

  // QM: test ómnibus de los moderadores (columnas 1..p−1), Wald conjunto = β_m' (A_mm)⁻¹ β_m.
  const gl = p - 1;
  let qm = 0;
  if (gl >= 1) {
    const Amm = Array.from({ length: gl }, (_, a) => Array.from({ length: gl }, (_, b) => A[a + 1][b + 1]));
    const AmmInv = invertirMatriz(Amm);
    if (AmmInv) {
      const bm = beta.slice(1);
      for (let a = 0; a < gl; a++) for (let b = 0; b < gl; b++) qm += bm[a] * AmmInv[a][b] * bm[b];
    }
  }

  const t2Total = remlTau2Regresion(y, v, X.map(() => [1]));
  const r2 = t2Total > 0 ? Math.max(0, ((t2Total - t2) / t2Total) * 100) : 0;

  return {
    coeficientes,
    tau2Residual: t2,
    qm,
    qmGl: gl,
    qmPValor: chiSquarePValue(qm, gl),
    r2,
    k: n,
  };
}

export interface EspecModerador {
  id: string;
  etiqueta: string;
  tipo: "continuo" | "categorico";
}

/** Moderadores disponibles para la meta-regresión, tomados de los campos de cada estudio. */
export const MODERADORES_DISPONIBLES: EspecModerador[] = [
  { id: "año", etiqueta: "Año de publicación", tipo: "continuo" },
  { id: "n", etiqueta: "Tamaño muestral (N)", tipo: "continuo" },
  { id: "diseño", etiqueta: "Diseño del estudio", tipo: "categorico" },
  { id: "calidad", etiqueta: "Riesgo de sesgo", tipo: "categorico" },
  { id: "evidencia", etiqueta: "Certeza (GRADE)", tipo: "categorico" },
];

function valorModerador(e: Estudio, id: string): string | undefined {
  const map: Record<string, string | undefined> = {
    año: e.año,
    n: e.n,
    diseño: e.diseño,
    calidad: e.calidad,
    evidencia: e.evidencia,
  };
  return map[id];
}

/**
 * Meta-regresión desde los campos de los estudios: construye la matriz de diseño (intercepto + moderadores
 * continuos centrados + categóricos codificados como dummies, referencia = primer nivel alfabético) y llama a
 * `ajustarMetaRegresion` (REML). Excluye estudios sin todos los moderadores seleccionados. `null` si no hay
 * datos suficientes (menos estudios que parámetros).
 */
export function metaRegresionMulti(
  estudios: Estudio[],
  medida: Medida,
  moderadores: EspecModerador[],
): ResultadoMetaRegresionMulti | null {
  if (moderadores.length === 0 || !escalaPorMedida[medida]) return null;
  const validos = extraerEstudiosValidos(estudios, medida);
  const porId = new Map(estudios.map((e) => [e.id, e]));

  const filas: { y: number; se: number; vals: string[] }[] = [];
  for (const v of validos) {
    const est = porId.get(v.id);
    if (!est) continue;
    const vals: string[] = [];
    let ok = true;
    for (const m of moderadores) {
      const raw = valorModerador(est, m.id);
      if (raw == null || raw === "" || (m.tipo === "continuo" && !Number.isFinite(parseFloat(raw)))) {
        ok = false;
        break;
      }
      vals.push(raw);
    }
    if (ok) filas.push({ y: v.y, se: v.se, vals });
  }
  if (filas.length < 2) return null;

  const nombres = ["Intercepto"];
  const X: number[][] = filas.map(() => [1]);
  moderadores.forEach((m, mi) => {
    if (m.tipo === "continuo") {
      const xs = filas.map((f) => parseFloat(f.vals[mi]));
      const media = xs.reduce((a, b) => a + b, 0) / xs.length;
      nombres.push(m.etiqueta);
      filas.forEach((_, i) => X[i].push(xs[i] - media));
    } else {
      const niveles = [...new Set(filas.map((f) => f.vals[mi]))].sort();
      niveles.slice(1).forEach((nivel) => {
        nombres.push(`${m.etiqueta}: ${nivel}`);
        filas.forEach((f, i) => X[i].push(f.vals[mi] === nivel ? 1 : 0));
      });
    }
  });
  if (filas.length <= nombres.length) return null;

  return ajustarMetaRegresion(
    filas.map((f) => ({ y: f.y, se: f.se })),
    X,
    nombres,
  );
}

/**
 * **Motor del meta-análisis.** Calcula efecto combinado, heterogeneidad, IC, etc. para una medida.
 *
 * Modelos:
 * - **Fijo:** asume un único efecto verdadero; usado si heterogeneidad baja (I² < 25).
 * - **Aleatorio:** asume efectos varían entre estudios; recomendado si I² > 25. Estima τ² (varianza entre-estudios).
 *
 * Método τ²:
 * - **DL** (DerSimonian-Laird): rápido, default. Usa Q de Cochran.
 * - **PM** (Paule-Mandel): alternativa.
 * - **REML** (máxima verosimilitud restringida): más preciso, más lento. Recomendado si k < 5.
 *
 * IC:
 * - Fijo: usa z-normal (1.96 para α=0.05).
 * - Aleatorio: usa t-distribution (Knapp-Hartung), más conservador.
 *
 * @param estudios - Array de Estudio con medida, y (efecto log), se (error estándar)
 * @param medida - La medida a analizar (OR, RR, DM, DME, etc.)
 * @param metodoTau2 - "DL", "PM", o "REML". Default: "DL"
 * @returns ResultadoMetaAnalisis (fijo, aleatorio, estudios pooled, I², Q, τ², etc.) o null si <2 estudios válidos
 *
 * @note Validado contra R `meta::metagen(method.tau="REML", hakn=TRUE)` con coincidencia a todos los decimales.
 *       Los números del forest plot y GRADE dependen de esta función.
 *
 * @example
 * const resultado = calcularMetaAnalisis(estudios, "OR", "REML");
 * console.log(resultado?.aleatorio.estimacion, resultado?.heterogeneidad.i2, resultado?.heterogeneidad.tau2);
 */
export function calcularMetaAnalisis(
  estudios: Estudio[],
  medida: Medida,
  metodoTau2: MetodoTau2 = "DL",
): ResultadoMetaAnalisis | null {
  const escala = escalaPorMedida[medida];
  if (!escala) return null;

  const validos = extraerEstudiosValidos(estudios, medida);
  if (validos.length < 2) return null;

  const pesosFijos = validos.map((v) => 1 / (v.se * v.se));
  const sumaPesosFijos = pesosFijos.reduce((a, b) => a + b, 0);
  const yFijo = validos.reduce((acc, v, i) => acc + pesosFijos[i] * v.y, 0) / sumaPesosFijos;
  const seFijo = Math.sqrt(1 / sumaPesosFijos);

  const q = validos.reduce((acc, v, i) => acc + pesosFijos[i] * (v.y - yFijo) ** 2, 0);
  const gl = validos.length - 1;
  const sumaPesos2 = pesosFijos.reduce((a, b) => a + b * b, 0);
  const c = sumaPesosFijos - sumaPesos2 / sumaPesosFijos;
  const tau2DL = c > 0 ? Math.max(0, (q - gl) / c) : 0;
  const tau2 = estimarTau2(validos, metodoTau2, tau2DL);
  const i2 = q > 0 ? Math.max(0, ((q - gl) / q) * 100) : 0;
  const pValor = chiSquarePValue(q, gl);

  const pesosAleatorios = validos.map((v) => 1 / (v.se * v.se + tau2));
  const sumaPesosAleatorios = pesosAleatorios.reduce((a, b) => a + b, 0);
  const yAleatorio = validos.reduce((acc, v, i) => acc + pesosAleatorios[i] * v.y, 0) / sumaPesosAleatorios;
  const seAleatorio = Math.sqrt(1 / sumaPesosAleatorios);

  // Knapp-Hartung (HKSJ) adjustment: a t-based CI for the random-effects estimate using a modified
  // SE (Hartung & Knapp 2001; IntHout 2014). qHK scales the variance; the t_{k-1} quantile widens
  // the CI, guarding against over-confidence — especially with few studies. gl = k-1 ≥ 1 here.
  const qHK = validos.reduce((acc, v, i) => acc + pesosAleatorios[i] * (v.y - yAleatorio) ** 2, 0) / gl;
  const seHK = Math.sqrt(qHK / sumaPesosAleatorios);
  const tCritHK = studentTQuantile975(gl);
  const hkIcInferior = destransformar(yAleatorio - tCritHK * seHK, escala);
  const hkIcSuperior = destransformar(yAleatorio + tCritHK * seHK, escala);

  // 95% prediction interval (Higgins/IntHout): needs at least 3 studies (gl = k-2).
  let prediccionInferior: number | null = null;
  let prediccionSuperior: number | null = null;
  if (validos.length >= 3) {
    const glPrediccion = validos.length - 2;
    const tCrit = studentTQuantile975(glPrediccion);
    const margen = tCrit * Math.sqrt(tau2 + seAleatorio * seAleatorio);
    prediccionInferior = destransformar(yAleatorio - margen, escala);
    prediccionSuperior = destransformar(yAleatorio + margen, escala);
  }

  const estudiosPooled: EstudioPooled[] = validos.map((v, i) => ({
    id: v.id,
    etiqueta: v.etiqueta,
    estimacion: v.estimacion,
    icInferior: v.icInferior,
    icSuperior: v.icSuperior,
    pesoFijo: (pesosFijos[i] / sumaPesosFijos) * 100,
    pesoAleatorio: (pesosAleatorios[i] / sumaPesosAleatorios) * 100,
    y: v.y,
    se: v.se,
  }));

  return {
    medida,
    escala,
    estudios: estudiosPooled,
    fijo: {
      estimacion: destransformar(yFijo, escala),
      icInferior: destransformar(yFijo - Z95 * seFijo, escala),
      icSuperior: destransformar(yFijo + Z95 * seFijo, escala),
    },
    aleatorio: {
      estimacion: destransformar(yAleatorio, escala),
      icInferior: destransformar(yAleatorio - Z95 * seAleatorio, escala),
      icSuperior: destransformar(yAleatorio + Z95 * seAleatorio, escala),
      prediccionInferior,
      prediccionSuperior,
      hkIcInferior,
      hkIcSuperior,
    },
    excluidosSinEventos: estudios
      .filter((e) => e.medida === medida && excluidoPorSinEventos(e))
      .map((e) => `${e.autor || "Estudio"}${e.año ? ` (${e.año})` : ""}`),
    heterogeneidad: { q, gl, pValor, i2, tau2, tau: Math.sqrt(tau2), metodoTau2 },
  };
}

/** Runs calcularMetaAnalisis for every poolable measure present, returning only groups with enough valid data. */
/**
 * **Orquestador del meta-análisis.** Calcula meta-análisis para TODAS las medidas presentes en los estudios.
 *
 * Flujo:
 * 1. Itera sobre cada medida poolable (OR, RR, DM, DME, etc.)
 * 2. Filtra estudios que tengan esa medida
 * 3. Llama a `calcularMetaAnalisis()` para esa medida
 * 4. Retorna array de resultados (uno por medida)
 *
 * Usado en:
 * - AnalysisStep: renderiza forest/funnel para cada resultado
 * - RetoGuiado: calcula demo en vivo
 * - Reportes: resume hallazgos por medida
 *
 * @param estudios - Array de Estudio con medida y efecto
 * @param metodoTau2 - Método para estimar τ². Default: "DL"
 * @returns Array de ResultadoMetaAnalisis (filtrados, solo medidas con ≥2 estudios)
 *
 * @example
 * const resultados = agruparPorMedida(estudios, "REML");
 * // [{medida: "OR", fijo: {...}, aleatorio: {...}, ...}, {medida: "DME", ...}]
 * const orResult = resultados.find(r => r.medida === "OR");
 */
export function agruparPorMedida(estudios: Estudio[], metodoTau2: MetodoTau2 = "DL"): ResultadoMetaAnalisis[] {
  const medidasPoolables = Object.keys(escalaPorMedida) as Medida[];
  return medidasPoolables
    .map((m) => calcularMetaAnalisis(estudios, m, metodoTau2))
    .filter((r): r is ResultadoMetaAnalisis => r !== null);
}

export interface FilaAcumulativa {
  /** Estudio que se acaba de incorporar en este paso. */
  etiqueta: string;
  año: string;
  /** Nº de estudios acumulados hasta aquí. */
  k: number;
  estimacion: number;
  icInferior: number;
  icSuperior: number;
  i2: number;
}

/**
 * **Meta-análisis acumulativo**: ordena los estudios (por defecto por año, ascendente) y recalcula el efecto
 * combinado añadiendo un estudio a la vez, mostrando cómo se ha ido acumulando la evidencia. Reusa el motor ya
 * validado (`calcularMetaAnalisis`) sobre cada prefijo; el primer paso (k=1) es el estudio solo. Replica
 * `metafor::cumul`. `null` si hay <2 estudios válidos con esta medida.
 */
export function metaAnalisisAcumulativo(
  estudios: Estudio[],
  medida: Medida,
  modelo: "fijo" | "aleatorio",
  metodoTau2: MetodoTau2 = "DL",
): FilaAcumulativa[] | null {
  const escala = escalaPorMedida[medida];
  if (!escala) return null;
  const validos = extraerEstudiosValidos(estudios, medida);
  if (validos.length < 2) return null;

  const idsValidos = new Set(validos.map((v) => v.id));
  const ordenados = estudios
    .filter((e) => e.medida === medida && idsValidos.has(e.id))
    .map((e, i) => ({ e, i }))
    .sort((a, b) => {
      const ya = parseFloat(a.e.año);
      const yb = parseFloat(b.e.año);
      const na = Number.isFinite(ya) ? ya : Infinity;
      const nb = Number.isFinite(yb) ? yb : Infinity;
      return na - nb || a.i - b.i; // orden estable: empates conservan el orden original
    })
    .map((x) => x.e);

  const filas: FilaAcumulativa[] = [];
  for (let k = 1; k <= ordenados.length; k++) {
    const añadido = ordenados[k - 1];
    const base = { etiqueta: añadido.autor || "Estudio", año: añadido.año, k };
    if (k === 1) {
      const uno = extraerEstudiosValidos([añadido], medida)[0];
      filas.push({
        ...base,
        estimacion: destransformar(uno.y, escala),
        icInferior: destransformar(uno.y - Z95 * uno.se, escala),
        icSuperior: destransformar(uno.y + Z95 * uno.se, escala),
        i2: 0,
      });
    } else {
      const r = calcularMetaAnalisis(ordenados.slice(0, k), medida, metodoTau2)!;
      const pooled = modelo === "aleatorio" ? r.aleatorio : r.fijo;
      filas.push({
        ...base,
        estimacion: pooled.estimacion,
        icInferior: pooled.icInferior,
        icSuperior: pooled.icSuperior,
        i2: r.heterogeneidad.i2,
      });
    }
  }
  return filas;
}

export interface FilaSensibilidad {
  idOmitido: number;
  etiqueta: string;
  estimacion: number;
  icInferior: number;
  icSuperior: number;
  i2: number;
}

/**
 * Análisis de sensibilidad «dejar uno fuera»: repite el meta-análisis omitiendo cada estudio por turno, para ver si un solo
 * estudio determina el resultado global. Requiere ≥ 3 estudios válidos (así quedan ≥ 2 al quitar uno). Usa el mismo motor
 * y el modelo fijo o aleatorio elegido.
 */
export function analisisSensibilidad(
  estudios: Estudio[],
  medida: Medida,
  modelo: "fijo" | "aleatorio",
  metodoTau2: MetodoTau2 = "DL",
): FilaSensibilidad[] | null {
  const base = calcularMetaAnalisis(estudios, medida, metodoTau2);
  if (!base || base.estudios.length < 3) return null;

  const filas: FilaSensibilidad[] = [];
  for (const s of base.estudios) {
    const parcial = calcularMetaAnalisis(
      estudios.filter((e) => e.id !== s.id),
      medida,
      metodoTau2,
    );
    if (!parcial) continue;
    const pooled = modelo === "aleatorio" ? parcial.aleatorio : parcial.fijo;
    filas.push({
      idOmitido: s.id,
      etiqueta: s.etiqueta,
      estimacion: pooled.estimacion,
      icInferior: pooled.icInferior,
      icSuperior: pooled.icSuperior,
      i2: parcial.heterogeneidad.i2,
    });
  }
  return filas;
}

export interface ResultadoEgger {
  intercepto: number;
  se: number;
  t: number;
  gl: number;
  pValor: number;
  k: number;
}

/**
 * Egger's regression test for funnel-plot asymmetry / small-study effects (Egger 1997). Regresses the
 * standard normal deviate (yᵢ/SEᵢ) on precision (1/SEᵢ) by OLS; a non-zero intercept signals asymmetry.
 * Returns null with fewer than 3 studies (needs df = k−2 ≥ 1). Note: for binary effect measures Egger
 * has an inflated false-positive rate — the funnel and this test should be read with that caveat.
 */
/**
 * **Test de Egger:** regresión lineal para detectar sesgo de publicación en funnel plots.
 *
 * Hipótesis nula: asimetría del funnel = 0 (no hay sesgo).
 * Regresiona: (efecto / SE) ~ 1/SE. Si p < 0.1, evidencia de asimetría.
 *
 * Interpretación:
 * - Intercepto (bias): magnitud de la asimetría
 * - p-valor < 0.1: evidencia estadística de asimetría
 * - Solo válido con k ≥ 10 estudios (Cochrane Handbook)
 *
 * Limitación: puede ser falso positivo si hay heterogeneidad o estudios pequeños de baja calidad.
 *
 * @param resultado - Meta-análisis ya calculado
 * @returns {bias (intercepto), seBias, tStat, pValor} o null si k < 3
 *
 * @example
 * const egger = testEgger(resultado);
 * if (egger?.pValor < 0.1) console.log("Evidencia de asimetría (posible sesgo)");
 */
export function testEgger(resultado: ResultadoMetaAnalisis): ResultadoEgger | null {
  const est = resultado.estudios;
  const k = est.length;
  if (k < 3) return null;

  const xs = est.map((e) => 1 / e.se); // precision
  const ys = est.map((e) => e.y / e.se); // standard normal deviate
  const mediaX = xs.reduce((a, b) => a + b, 0) / k;
  const mediaY = ys.reduce((a, b) => a + b, 0) / k;

  let sxx = 0;
  let sxy = 0;
  for (let i = 0; i < k; i++) {
    sxx += (xs[i] - mediaX) ** 2;
    sxy += (xs[i] - mediaX) * (ys[i] - mediaY);
  }
  // Precisión (casi) constante → la regresión no tiene pendiente estimable. `sxx === 0` no lo detecta
  // por el ruido de coma flotante, así que se compara contra la escala de la propia precisión.
  if (!(sxx / k > 1e-12 * mediaX * mediaX)) return null;

  const pendiente = sxy / sxx;
  const intercepto = mediaY - pendiente * mediaX;

  let sse = 0;
  for (let i = 0; i < k; i++) {
    const pred = intercepto + pendiente * xs[i];
    sse += (ys[i] - pred) ** 2;
  }
  const gl = k - 2;
  const sigma2 = sse / gl;
  const se = Math.sqrt(sigma2 * (1 / k + (mediaX * mediaX) / sxx));
  const t = se > 0 ? intercepto / se : 0;
  const pValor = studentTDosColas(Math.abs(t), gl);

  return { intercepto, se, t, gl, pValor, k };
}

// --- Test de Begg & Mazumdar (correlación de rangos de Kendall) para asimetría del funnel ---

/**
 * p exacta de dos lados de la τ de Kendall vía la distribución de inversiones (números de Mahon),
 * como `cor.test(method="kendall")` de R para k pequeño sin empates. Segura en doble precisión hasta k≈18.
 */
function kendallPExacta(S: number, k: number): number {
  const total = (k * (k - 1)) / 2; // nº de pares
  const D = Math.round((total - S) / 2); // discordancias observadas
  let dist = [1]; // I(1, ·)
  for (let m = 2; m <= k; m++) {
    const maxk = (m * (m - 1)) / 2;
    const prefix: number[] = [];
    let acc = 0;
    for (let i = 0; i < dist.length; i++) {
      acc += dist[i];
      prefix[i] = acc;
    }
    const sumaHasta = (idx: number) => (idx < 0 ? 0 : prefix[Math.min(idx, prefix.length - 1)]);
    const next = new Array<number>(maxk + 1);
    for (let key = 0; key <= maxk; key++) next[key] = sumaHasta(key) - sumaHasta(key - m);
    dist = next;
  }
  const nFact = dist.reduce((a, b) => a + b, 0); // = k!
  const menor = Math.min(D, total - D);
  let cum = 0;
  for (let key = 0; key <= menor; key++) cum += dist[key];
  return Math.min(1, (2 * cum) / nFact);
}

export interface ResultadoBegg {
  /** τ de Kendall entre las varianzas y los efectos estandarizados (τ-b si hay empates). */
  tau: number;
  pValor: number;
  /** true si la p es exacta (k≤18 y sin empates); false si es la aproximación normal. */
  exacta: boolean;
  /** Hay valores repetidos (p. ej. estudios con el mismo tamaño de muestra y por tanto la misma varianza): se usa τ-b y la varianza con empates. */
  empates: boolean;
  k: number;
}

/** Dos valores que difieren solo por el ruido de redondeo de doble precisión cuentan como iguales (empate). */
const IGUALES = (a: number, b: number) => Math.abs(a - b) <= 1e-9 * Math.max(Math.abs(a), Math.abs(b), 1e-300);

/**
 * Tamaños de los grupos de valores iguales (los de tamaño 1 se omiten): lo que necesitan la τ-b y la varianza de Kendall con empates.
 * Las varianzas de estudios con el mismo tamaño se reconstruyen desde sus IC y pueden diferir en el último dígito; se tratan como iguales.
 */
function tamanosDeEmpates(valores: number[]): number[] {
  const orden = [...valores].sort((a, b) => a - b);
  const tamanos: number[] = [];
  let actual = 1;
  for (let i = 1; i < orden.length; i++) {
    if (IGUALES(orden[i], orden[i - 1])) actual++;
    else {
      if (actual > 1) tamanos.push(actual);
      actual = 1;
    }
  }
  if (actual > 1) tamanos.push(actual);
  return tamanos;
}

/**
 * **Test de rangos de Begg & Mazumdar (1994)** para asimetría del funnel plot. Estandariza cada efecto
 * `(yᵢ−θ̂_fijo)/√(vᵢ−v̄)` (con `v̄=1/Σwᵢ`) y calcula la **τ de Kendall** entre esos valores y las varianzas
 * `vᵢ`. Réplica de `metafor::ranktest`. p exacta (distribución de inversiones) para k≤18, aproximación normal
 * para k mayor. `null` si <3 estudios. Menos potente que Egger con pocos estudios; complementa, no reemplaza.
 */
/**
 * **Test de Begg:** correlación de Kendall entre efecto y varianza.
 *
 * Hipótesis nula: no hay correlación (no hay sesgo de publicación).
 * Si p < 0.1, posible sesgo: estudios pequeños subestimados / no publicados.
 *
 * Ventaja vs. Egger: menos sensible a heterogeneidad.
 * Limitación: requiere k ≥ 10 para potencia adecuada.
 *
 * @param resultado - Meta-análisis ya calculado
 * @returns {tau (correlación), zStat, pValor} o null si k < 3
 */
export function testBegg(resultado: ResultadoMetaAnalisis): ResultadoBegg | null {
  const pts = resultado.estudios.map((e) => ({ y: e.y, v: e.se * e.se }));
  const k = pts.length;
  if (k < 3) return null;

  const sw = pts.reduce((s, p) => s + 1 / p.v, 0);
  const vBar = 1 / sw;
  const theta = pts.reduce((s, p) => s + p.y / p.v, 0) / sw;
  const star = pts.map((p) => ({ v: p.v, e: (p.y - theta) / Math.sqrt(p.v - vBar) }));

  let S = 0;
  for (let i = 0; i < k; i++) {
    for (let j = i + 1; j < k; j++) {
      const dv = IGUALES(star[j].v, star[i].v) ? 0 : star[j].v - star[i].v;
      const de = IGUALES(star[j].e, star[i].e) ? 0 : star[j].e - star[i].e;
      S += Math.sign(dv * de);
    }
  }
  const pares = (k * (k - 1)) / 2;
  const tx = tamanosDeEmpates(star.map((s) => s.v));
  const ty = tamanosDeEmpates(star.map((s) => s.e));
  const empates = tx.length > 0 || ty.length > 0;
  // Sin empates, τ = S / nº de pares. Con empates, τ-b = S / √((pares − Σtₓ(tₓ−1)/2)(pares − Σt_y(t_y−1)/2)), como `cor(method = "kendall")`.
  const n1 = tx.reduce((s, c) => s + (c * (c - 1)) / 2, 0);
  const n2 = ty.reduce((s, c) => s + (c * (c - 1)) / 2, 0);
  const denominador = Math.sqrt((pares - n1) * (pares - n2));
  const tau = denominador > 0 ? S / denominador : 0;

  let pValor: number;
  let exacta: boolean;
  if (k <= 18 && !empates) {
    pValor = kendallPExacta(S, k);
    exacta = true;
  } else {
    // Varianza de S bajo H₀ con empates (Kendall 1970, la que usa `cor.test(method = "kendall")`):
    //   [k(k−1)(2k+5) − Σt(t−1)(2t+5) − Σu(u−1)(2u+5)]/18 + Σt(t−1)(t−2)·Σu(u−1)(u−2)/(9k(k−1)(k−2)) + Σt(t−1)·Σu(u−1)/(2k(k−1)).
    const f = (xs: number[], g: (c: number) => number) => xs.reduce((s, c) => s + g(c), 0);
    const varianza =
      (k * (k - 1) * (2 * k + 5) - f(tx, (c) => c * (c - 1) * (2 * c + 5)) - f(ty, (c) => c * (c - 1) * (2 * c + 5))) / 18 +
      (k > 2 ? (f(tx, (c) => c * (c - 1) * (c - 2)) * f(ty, (c) => c * (c - 1) * (c - 2))) / (9 * k * (k - 1) * (k - 2)) : 0) +
      (f(tx, (c) => c * (c - 1)) * f(ty, (c) => c * (c - 1))) / (2 * k * (k - 1));
    const z = varianza > 0 ? S / Math.sqrt(varianza) : 0;
    pValor = varianza > 0 ? Math.min(1, colaNormalDosLados(z)) : 1;
    exacta = false;
  }
  return { tau, pValor, exacta, empates, k };
}

export interface ResultadoPeters {
  /** Slope β₁ of ln(OR) ~ 1/n; asymmetry when it differs from 0. */
  pendiente: number;
  se: number;
  t: number;
  gl: number;
  pValor: number;
  /** Nº de estudios con tabla 2×2 utilizable. */
  k: number;
}

/**
 * Peters' regression test for publication bias in binary (OR) meta-analyses.
 * Weighted regression of ln(OR) on 1/n (total sample size), weights w = 1/(1/(a+c) + 1/(b+d));
 * tests the SLOPE. Requires the studies' raw 2×2 counts (`Estudio.crudos`). Only for the OR measure.
 * Ref: Harrer et al., Doing Meta-Analysis in R (ec. 9.3–9.4); Peters et al. 2006.
 */
/**
 * **Test de Peters:** regresión robusta para sesgo de publicación (alternativa a Egger).
 *
 * Regresiona efecto vs. 1/N (tamaño inverso). Menos sensible a outliers que Egger.
 * Recomendado por Cochrane Handbook como complemento a Egger/Begg.
 *
 * @param estudios - Array de estudios con n (tamaño de muestra total)
 * @param medida - Medida a analizar
 * @returns {bias, seBias, tStat, pValor} o null si k < 3 o estudios sin n
 */
export function testPeters(estudios: Estudio[], medida: Medida): ResultadoPeters | null {
  if (medida !== "OR") return null;

  const pts: { y: number; x: number; w: number }[] = [];
  estudios
    .filter((e) => e.medida === medida)
    .forEach((e) => {
      const cr = e.crudos;
      if (!cr || cr.tipo !== "2x2") return;
      const kCtrl = divisorControl(e.dividirControlEn);
      const a = parseFloat(cr.eventosTrat);
      const n1 = parseFloat(cr.totalTrat);
      const c = parseFloat(cr.eventosCtrl) / kCtrl;
      const n2 = parseFloat(cr.totalCtrl) / kCtrl;
      if (![a, n1, c, n2].every(Number.isFinite)) return;
      if (n1 <= 0 || n2 <= 0 || a < 0 || c < 0 || a > n1 || c > n2) return;

      let A = a;
      let B = n1 - a;
      let C = c;
      let D = n2 - c;
      // Peters et al. (2006) y `meta::metabias`: el +0,5 corrige SOLO el ln(OR); el peso 1/(1/(a+c)+1/(b+d)) usa los conteos
      // originales (con una celda cero sigue bien definido; un doble-cero queda con peso 0 y no aporta).
      const w = 1 / (1 / (A + C) + 1 / (B + D));
      if (A === 0 || B === 0 || C === 0 || D === 0) {
        A += 0.5;
        B += 0.5;
        C += 0.5;
        D += 0.5;
      }
      const y = Math.log((A * D) / (B * C));
      if (!Number.isFinite(y) || !Number.isFinite(w) || w <= 0) return;
      pts.push({ y, x: 1 / (n1 + n2), w });
    });

  const k = pts.length;
  if (k < 3) return null;

  const sumaW = pts.reduce((s, p) => s + p.w, 0);
  const mediaX = pts.reduce((s, p) => s + p.w * p.x, 0) / sumaW;
  const mediaY = pts.reduce((s, p) => s + p.w * p.y, 0) / sumaW;
  let sxx = 0;
  let sxy = 0;
  for (const p of pts) {
    sxx += p.w * (p.x - mediaX) ** 2;
    sxy += p.w * (p.x - mediaX) * (p.y - mediaY);
  }
  // Si todos los estudios tienen (casi) el mismo tamaño, 1/n es constante y la pendiente no se puede
  // estimar. Comparar `sxx === 0` no basta: el ruido de coma flotante deja sxx ≈ 1e-19 y la regresión
  // devolvería una pendiente enorme con p = 1 ("sin asimetría") que en realidad no significa nada.
  if (!(sxx / sumaW > 1e-12 * mediaX * mediaX)) return null;

  const pendiente = sxy / sxx;
  const intercepto = mediaY - pendiente * mediaX;
  let sse = 0;
  for (const p of pts) {
    const r = p.y - (intercepto + pendiente * p.x);
    sse += p.w * r * r;
  }
  const gl = k - 2;
  const se = Math.sqrt(sse / gl / sxx);
  const t = se > 0 ? pendiente / se : 0;
  const pValor = studentTDosColas(Math.abs(t), gl);

  return { pendiente, se, t, gl, pValor, k };
}

// --- Mantel-Haenszel: pooling para datos binarios 2×2, recomendado con eventos raros / muestras pequeñas ---

/** Tabla 2×2 de un estudio, ya parseada y validada. */
interface Tabla2x2 {
  a: number; // eventos en tratamiento
  b: number; // no-eventos en tratamiento
  c: number; // eventos en control
  d: number; // no-eventos en control
  n1: number; // total tratamiento (a+b)
  n2: number; // total control (c+d)
  N: number; // total (n1+n2)
}

export interface ResultadoMH {
  medida: Medida;
  escala: EscalaAnalisis;
  /** Efecto combinado en su escala natural (OR/RR o diferencia de riesgos). */
  estimacion: number;
  icInferior: number;
  icSuperior: number;
  /** SE en la escala de análisis (log para OR/RR, lineal para RD). */
  se: number;
  pValor: number;
  k: number;
  /** Estudios descartados (sin crudos 2×2 válidos, o doble-cero no informativos). */
  excluidos: number;
  /** Heterogeneidad por inverso-varianza sobre los efectos individuales (igual que metafor junto a MH). */
  heterogeneidad: { q: number; gl: number; pValor: number; i2: number };
}

/** Extrae y valida las tablas 2×2 de los estudios de una medida binaria. Descarta dobles-ceros no informativos. */
function extraerTablas(estudios: Estudio[], medida: Medida): Tabla2x2[] {
  const tablas: Tabla2x2[] = [];
  estudios
    .filter((e) => e.medida === medida)
    .forEach((e) => {
      const cr = e.crudos;
      if (!cr || cr.tipo !== "2x2") return;
      const kCtrl = divisorControl(e.dividirControlEn);
      const a = parseFloat(cr.eventosTrat);
      const n1 = parseFloat(cr.totalTrat);
      const c = parseFloat(cr.eventosCtrl) / kCtrl;
      const n2 = parseFloat(cr.totalCtrl) / kCtrl;
      if (![a, n1, c, n2].every(Number.isFinite)) return;
      if (n1 <= 0 || n2 <= 0 || a < 0 || c < 0 || a > n1 || c > n2) return;
      const b = n1 - a;
      const d = n2 - c;
      // Doble-cero (sin eventos en ambos grupos, o sin no-eventos en ambos). Para OR y RR no aporta nada al estimador de
      // Mantel-Haenszel (todos sus términos valen 0), así que se descarta. Para la diferencia de riesgos SÍ cuenta: es
      // un estudio con RD = 0 y peso n1·n2/N, y excluirlo desplaza el resultado lejos de 0 (Cochrane Handbook cap. 10:
      // el estimador de MH está definido con doble-cero y no es necesario ni aconsejable excluirlos; metafor los incluye).
      if (medida !== "RAR" && ((a === 0 && c === 0) || (b === 0 && d === 0))) return;
      tablas.push({ a, b, c, d, n1, n2, N: n1 + n2 });
    });
  return tablas;
}

/**
 * Heterogeneidad (Q, I²) sobre los efectos individuales, con corrección de continuidad +0,5 aplicada SOLO a
 * los estudios con alguna celda cero (como metafor: `to="only0"`). La Q se **centra en el estimador de
 * Mantel-Haenszel** (`centro`, en escala de análisis), no en la media inverso-varianza: así se reproduce
 * exactamente el Q que `metafor::rma.mh` reporta junto a la estimación MH.
 */
function heterogeneidadBinaria(
  tablas: Tabla2x2[],
  medida: Medida,
  centro: number,
): { q: number; gl: number; pValor: number; i2: number } {
  const pts: PuntoEfecto[] = tablas.map((t) => {
    let { a, b, c, d } = t;
    if (a === 0 || b === 0 || c === 0 || d === 0) {
      a += 0.5;
      b += 0.5;
      c += 0.5;
      d += 0.5;
    }
    const n1 = a + b;
    const n2 = c + d;
    if (medida === "OR") return { y: Math.log((a * d) / (b * c)), se: Math.sqrt(1 / a + 1 / b + 1 / c + 1 / d) };
    if (medida === "RR")
      return { y: Math.log((a / n1) / (c / n2)), se: Math.sqrt(1 / a - 1 / n1 + 1 / c - 1 / n2) };
    // RAR = diferencia de riesgos (RD)
    const rd = a / n1 - c / n2;
    return { y: rd, se: Math.sqrt((a * b) / n1 ** 3 + (c * d) / n2 ** 3) };
  });
  const q = pts.reduce((s, p) => s + (1 / (p.se * p.se)) * (p.y - centro) ** 2, 0);
  const gl = pts.length - 1;
  const i2 = q > 0 ? Math.max(0, ((q - gl) / q) * 100) : 0;
  return { q, gl, pValor: chiSquarePValue(q, gl), i2 };
}

/**
 * Meta-análisis por el método de **Mantel-Haenszel** (efecto fijo) para datos binarios 2×2.
 * Preferible al inverso-varianza cuando hay **eventos raros o muestras pequeñas**, porque no depende de
 * varianzas intra-estudio de muestra grande y maneja bien las celdas cero (sin corrección de continuidad).
 * - OR: `OR_MH = Σ(aᵢdᵢ/Nᵢ) / Σ(bᵢcᵢ/Nᵢ)`; varianza de ln(OR) por Robins-Breslow-Greenland.
 * - RR: `RR_MH = Σ(aᵢn2ᵢ/Nᵢ) / Σ(cᵢn1ᵢ/Nᵢ)`; varianza de ln(RR) por Greenland-Robins.
 * - RD (RAR): `RD_MH = Σwᵢ(RDᵢ)/Σwᵢ`, `wᵢ=n1ᵢn2ᵢ/Nᵢ`; varianza por Greenland-Robins.
 * Fuentes: Cochrane Handbook §10.4 (Deeks & Higgins); Robins, Breslow & Greenland (1986); verificado contra
 * `metafor::rma.mh` (Rothman et al. 2008). `null` si <2 estudios binarios con tabla 2×2 válida.
 */
export function metaAnalisisMH(estudios: Estudio[], medida: Medida): ResultadoMH | null {
  if (medida !== "OR" && medida !== "RR" && medida !== "RAR") return null;
  const escala = escalaPorMedida[medida];
  if (!escala) return null;

  const totales = estudios.filter((e) => e.medida === medida && e.crudos?.tipo === "2x2").length;
  const tablas = extraerTablas(estudios, medida);
  if (tablas.length < 2) return null;

  let estimacionT: number; // en escala de análisis
  let se: number;

  if (medida === "OR") {
    let R = 0;
    let S = 0;
    let sPR = 0; // Σ Pᵢ Rᵢ
    let sPSQR = 0; // Σ (Pᵢ Sᵢ + Qᵢ Rᵢ)
    let sQS = 0; // Σ Qᵢ Sᵢ
    for (const t of tablas) {
      const Ri = (t.a * t.d) / t.N;
      const Si = (t.b * t.c) / t.N;
      const Pi = (t.a + t.d) / t.N;
      const Qi = (t.b + t.c) / t.N;
      R += Ri;
      S += Si;
      sPR += Pi * Ri;
      sPSQR += Pi * Si + Qi * Ri;
      sQS += Qi * Si;
    }
    estimacionT = Math.log(R / S);
    const varLog = sPR / (2 * R * R) + sPSQR / (2 * R * S) + sQS / (2 * S * S);
    se = Math.sqrt(varLog);
  } else if (medida === "RR") {
    let R = 0;
    let S = 0;
    let num = 0; // Σ (n1ᵢ n2ᵢ (aᵢ+cᵢ) − aᵢ cᵢ Nᵢ)/Nᵢ²
    for (const t of tablas) {
      R += (t.a * t.n2) / t.N;
      S += (t.c * t.n1) / t.N;
      num += (t.n1 * t.n2 * (t.a + t.c) - t.a * t.c * t.N) / (t.N * t.N);
    }
    estimacionT = Math.log(R / S);
    se = Math.sqrt(num / (R * S));
  } else {
    // RD (RAR): media ponderada de las diferencias de riesgo con pesos n1·n2/N.
    let sw = 0;
    let swRD = 0;
    let varNum = 0; // Σ (aᵢbᵢn2ᵢ³ + cᵢdᵢn1ᵢ³)/(n1ᵢ n2ᵢ Nᵢ²)
    for (const t of tablas) {
      const w = (t.n1 * t.n2) / t.N;
      sw += w;
      swRD += w * (t.a / t.n1 - t.c / t.n2);
      varNum += (t.a * t.b * t.n2 ** 3 + t.c * t.d * t.n1 ** 3) / (t.n1 * t.n2 * t.N * t.N);
    }
    estimacionT = swRD / sw;
    se = Math.sqrt(varNum / (sw * sw));
  }

  const estimacion = destransformar(estimacionT, escala);
  const icInferior = destransformar(estimacionT - Z95 * se, escala);
  const icSuperior = destransformar(estimacionT + Z95 * se, escala);
  const pValor = se > 0 ? Math.min(1, colaNormalDosLados(estimacionT / se)) : 1;

  return {
    medida,
    escala,
    estimacion,
    icInferior,
    icSuperior,
    se,
    pValor,
    k: tablas.length,
    excluidos: totales - tablas.length,
    heterogeneidad: heterogeneidadBinaria(tablas, medida, estimacionT),
  };
}

// --- Peto: OR de "un paso" para datos binarios 2×2 (eventos raros, grupos balanceados) ---

export interface ResultadoPeto {
  medida: "OR";
  escala: EscalaAnalisis;
  estimacion: number;
  icInferior: number;
  icSuperior: number;
  /** SE del ln(OR) de Peto. */
  se: number;
  pValor: number;
  k: number;
  excluidos: number;
  heterogeneidad: { q: number; gl: number; pValor: number; i2: number };
}

/**
 * Meta-análisis por el método de **Peto** (efecto fijo, solo OR) para datos binarios 2×2.
 * Recomendado con **eventos raros y grupos de tamaño similar** (Cochrane Handbook §10.4.2). Por estudio usa
 * la diferencia observado−esperado y la varianza hipergeométrica:
 *   Eᵢ = n1ᵢ(aᵢ+cᵢ)/Nᵢ ; Vᵢ = n1ᵢ·n2ᵢ·(aᵢ+cᵢ)·(bᵢ+dᵢ) / (Nᵢ²·(Nᵢ−1)) ; (O−E)ᵢ = aᵢ − Eᵢ
 *   ln(OR_Peto) = Σ(O−E)/ΣV ; var = 1/ΣV ; Q = Σ(O−E)²/V − (Σ(O−E))²/ΣV
 * No necesita corrección de continuidad (maneja celdas cero). Verificado contra `metafor::rma.peto`.
 * `null` si <2 estudios OR con tabla 2×2 utilizable. Advertencia: Peto se sesga con OR lejos de 1 o grupos
 * desbalanceados; en ese caso preferir Mantel-Haenszel o inverso-varianza.
 */
export function metaAnalisisPeto(estudios: Estudio[]): ResultadoPeto | null {
  const escala = escalaPorMedida.OR!;
  const totales = estudios.filter((e) => e.medida === "OR" && e.crudos?.tipo === "2x2").length;
  const tablas = extraerTablas(estudios, "OR");

  const puntos: { oe: number; v: number }[] = [];
  for (const t of tablas) {
    const eventos = t.a + t.c;
    const noEventos = t.b + t.d;
    const E = (t.n1 * eventos) / t.N;
    const V = (t.n1 * t.n2 * eventos * noEventos) / (t.N * t.N * (t.N - 1));
    if (!(V > 0)) continue; // sin varianza informativa (p. ej. todos eventos o todos no-eventos)
    puntos.push({ oe: t.a - E, v: V });
  }
  if (puntos.length < 2) return null;

  const sumOE = puntos.reduce((s, p) => s + p.oe, 0);
  const sumV = puntos.reduce((s, p) => s + p.v, 0);
  const lnOR = sumOE / sumV;
  const se = Math.sqrt(1 / sumV);

  const q = puntos.reduce((s, p) => s + (p.oe * p.oe) / p.v, 0) - (sumOE * sumOE) / sumV;
  const gl = puntos.length - 1;
  const i2 = q > 0 ? Math.max(0, ((q - gl) / q) * 100) : 0;

  return {
    medida: "OR",
    escala,
    estimacion: destransformar(lnOR, escala),
    icInferior: destransformar(lnOR - Z95 * se, escala),
    icSuperior: destransformar(lnOR + Z95 * se, escala),
    se,
    pValor: se > 0 ? Math.min(1, colaNormalDosLados(lnOR / se)) : 1,
    k: puntos.length,
    excluidos: totales - puntos.length,
    heterogeneidad: { q, gl, pValor: chiSquarePValue(q, gl), i2 },
  };
}

// --- Trim-and-fill (Duval & Tweedie 2000): estimate & adjust for studies missing due to publication bias ---

/** Rangos (1..n) de |vals|, con rangos promedio para empates. */
function rangosAbs(vals: number[]): number[] {
  const idx = vals.map((v, i) => ({ a: Math.abs(v), i })).sort((x, y) => x.a - y.a);
  const ranks = new Array<number>(vals.length);
  let i = 0;
  while (i < idx.length) {
    let j = i;
    while (j + 1 < idx.length && idx[j + 1].a === idx[i].a) j++;
    const promedio = (i + j) / 2 + 1;
    for (let t = i; t <= j; t++) ranks[idx[t].i] = promedio;
    i = j + 1;
  }
  return ranks;
}

/** Media ponderada (fija) de un conjunto de puntos {y, se}. */
function mediaFija(pts: readonly PuntoEfecto[]): number {
  const w = pts.map((p) => 1 / (p.se * p.se));
  const sw = w.reduce((a, b) => a + b, 0);
  return pts.reduce((acc, p, i) => acc + w[i] * p.y, 0) / sw;
}

/** Combina puntos {y, se} en escala transformada por el modelo elegido; devuelve {y, se} del efecto combinado. */
function combinarYSe(
  pts: readonly PuntoEfecto[],
  modelo: "fijo" | "aleatorio",
  metodoTau2: MetodoTau2,
): { y: number; se: number } {
  const wf = pts.map((p) => 1 / (p.se * p.se));
  const swf = wf.reduce((a, b) => a + b, 0);
  const yF = pts.reduce((acc, p, i) => acc + wf[i] * p.y, 0) / swf;
  if (modelo === "fijo") return { y: yF, se: Math.sqrt(1 / swf) };

  const gl = pts.length - 1;
  const q = pts.reduce((acc, p, i) => acc + wf[i] * (p.y - yF) ** 2, 0);
  const sw2 = wf.reduce((a, b) => a + b * b, 0);
  const c = swf - sw2 / swf;
  const tau2DL = c > 0 ? Math.max(0, (q - gl) / c) : 0;
  const tau2 = estimarTau2(pts, metodoTau2, tau2DL);
  const wr = pts.map((p) => 1 / (p.se * p.se + tau2));
  const swr = wr.reduce((a, b) => a + b, 0);
  const yR = pts.reduce((acc, p, i) => acc + wr[i] * p.y, 0) / swr;
  return { y: yR, se: Math.sqrt(1 / swr) };
}

export interface ResultadoTrimFill {
  /** Nº estimado de estudios "faltantes" por sesgo de publicación (estimador L0). */
  k0: number;
  /** Lado del funnel donde se imputan los estudios faltantes. */
  lado: "izquierda" | "derecha";
  estimacionAjustada: number;
  icInferiorAjustado: number;
  icSuperiorAjustado: number;
  /** Estudios imputados (efecto transformado y su SE), para dibujarlos en el funnel plot. */
  imputados: PuntoEfecto[];
}

/**
 * Trim-and-fill (Duval & Tweedie 2000, estimador L0). Estima cuántos estudios faltan por asimetría,
 * los imputa reflejando los más extremos del lado sobrerrepresentado (y* = 2μ − y) y recalcula el
 * efecto combinado ajustado. L0 = (4·Tₙ − k(k+1)) / (2k−1). Requiere ≥3 estudios.
 */
/**
 * **Trim-and-Fill:** método de Duval & Tweedie para ajustar sesgo de publicación en funnel plots.
 *
 * Detección de asimetría:
 * 1. Calcula suma de rangos signados de los estudios
 * 2. Si hay asimetría (lado sobrerrepresentado), asume estudios "faltantes"
 * 3. Iterativamente recorta (trim) estudios del lado + fuerte, estima media simetría
 * 4. Imputa (fill) reflejos de los cortados en el lado opuesto
 *
 * Uso en AnalysisStep:
 * - Dibuja círculos violeta huecos en el funnel (puntos imputados)
 * - Permite evaluar "¿cuánto cambiaría el resultado si se publicaran los estudios faltantes?"
 *
 * @param resultado - Meta-análisis ya calculado (con estudios, escala, etc.)
 * @param modelo - "fijo" o "aleatorio" (usado para recalcular efecto imputed)
 * @param metodoTau2 - Método para τ² si modelo="aleatorio". Default: "DL"
 * @returns {imputados, estimacionAjustada} o null si k < 3 o sin asimetría detectada
 *
 * @note Conservador: solo detecta asimetría con ≥3 estudios. No implica sesgo de publicación definitivo.
 *       (También puede deberse a heterogeneidad o estudios pequeños de peor calidad.)
 *
 * @example
 * const tfill = trimAndFill(resultado, "aleatorio", "REML");
 * if (tfill) {
 *   // Dibujar círculos huecos en funnel: tfill.imputados.map(p => ({y: p.y, se: p.se}))
 *   console.log("Efecto ajustado:", tfill.estimacionAjustada);
 * }
 */
export function trimAndFill(
  resultado: ResultadoMetaAnalisis,
  modelo: "fijo" | "aleatorio",
  metodoTau2: MetodoTau2 = "DL",
): ResultadoTrimFill | null {
  const escala = resultado.escala;
  const base: PuntoEfecto[] = resultado.estudios.map((e) => ({ y: e.y, se: e.se }));
  const k = base.length;
  if (k < 3) return null;

  // Centro de simetría: la media del MISMO modelo con el que se combina al final (efectos fijos: media ponderada; efectos
  // aleatorios: con su τ² estimado en cada paso), como hacen `metafor::trimfill` y `metatrim` de Stata.
  const centro = (pts: readonly PuntoEfecto[]): number =>
    modelo === "fijo" ? mediaFija(pts) : combinarYSe(pts, "aleatorio", metodoTau2).y;

  // Lado sobrerrepresentado por la suma de rangos signados respecto al centro.
  const mu0 = centro(base);
  const eBase = base.map((p) => p.y - mu0);
  const rBase = rangosAbs(eBase);
  const sr = eBase.reduce((acc, ei, i) => acc + Math.sign(ei) * rBase[i], 0);
  const overRight = sr >= 0;
  const s = overRight ? 1 : -1;
  // Orientamos para que el lado sobrerrepresentado sea siempre el positivo/derecho.
  const oriented: PuntoEfecto[] = base.map((p) => ({ y: s * p.y, se: p.se }));

  let k0 = 0;
  for (let iter = 0; iter < 100; iter++) {
    const sorted = [...oriented].sort((a, b) => a.y - b.y);
    const mu = centro(sorted.slice(0, k - k0)); // recorta los k0 mayores
    const e = oriented.map((p) => p.y - mu);
    const ranks = rangosAbs(e);
    let tn = 0;
    for (let i = 0; i < k; i++) if (e[i] > 0) tn += ranks[i];
    const l0 = (4 * tn - k * (k + 1)) / (2 * k - 1);
    const k0new = Math.max(0, Math.round(l0));
    if (k0new === k0) break;
    k0 = Math.min(k0new, k - 1);
  }

  const sortedFinal = [...oriented].sort((a, b) => a.y - b.y);
  const muFinal = centro(sortedFinal.slice(0, k - k0));
  const imputados: PuntoEfecto[] = sortedFinal
    .slice(k - k0)
    .map((p) => ({ y: s * (2 * muFinal - p.y), se: p.se }));

  const comb = combinarYSe([...base, ...imputados], modelo, metodoTau2);
  return {
    k0,
    lado: overRight ? "izquierda" : "derecha",
    estimacionAjustada: destransformar(comb.y, escala),
    icInferiorAjustado: destransformar(comb.y - Z95 * comb.se, escala),
    icSuperiorAjustado: destransformar(comb.y + Z95 * comb.se, escala),
    imputados,
  };
}

// --- Subgroup analysis: split by a categorical moderator and test whether the subgroups differ ---

export interface SubgrupoResultado {
  etiqueta: string;
  k: number;
  estimacion: number;
  icInferior: number;
  icSuperior: number;
  /** I² dentro del subgrupo; null si tiene <2 estudios. */
  i2: number | null;
}

export interface ResultadoSubgrupos {
  campo: CampoSubgrupo;
  subgrupos: SubgrupoResultado[];
  /** Q de diferencias entre subgrupos y su prueba (χ² con gl = nº subgrupos − 1). */
  qBet: number;
  gl: number;
  pValor: number;
}

/**
 * Análisis de subgrupos: combina los estudios dentro de cada nivel de un moderador categórico y prueba
 * si los subgrupos difieren con la Q entre subgrupos `Q_bet = Σ Wg(θg − θ̄)²` (χ², gl = G−1), como RevMan.
 * Requiere ≥2 subgrupos con datos válidos.
 */
export function analisisSubgrupos(
  estudios: Estudio[],
  medida: Medida,
  campo: CampoSubgrupo,
  modelo: "fijo" | "aleatorio",
  metodoTau2: MetodoTau2 = "DL",
): ResultadoSubgrupos | null {
  const escala = escalaPorMedida[medida];
  if (!escala) return null;
  const validos = extraerEstudiosValidos(estudios, medida);
  if (validos.length < 2) return null;

  const grupoDe = new Map<number, string>();
  for (const e of estudios) {
    const val = e[campo];
    grupoDe.set(e.id, val ? String(val) : "Sin especificar");
  }

  const grupos = new Map<string, EstudioValido[]>();
  for (const v of validos) {
    const g = grupoDe.get(v.id) ?? "Sin especificar";
    const arr = grupos.get(g);
    if (arr) arr.push(v);
    else grupos.set(g, [v]);
  }
  if (grupos.size < 2) return null;

  const puntos: PuntoEfecto[] = [];
  const subgrupos: SubgrupoResultado[] = [];
  for (const [etiqueta, vs] of grupos) {
    const comb = combinarYSe(vs, modelo, metodoTau2);
    puntos.push(comb);

    let i2: number | null = null;
    if (vs.length >= 2) {
      const wf = vs.map((v) => 1 / (v.se * v.se));
      const swf = wf.reduce((a, b) => a + b, 0);
      const yf = vs.reduce((a, v, i) => a + wf[i] * v.y, 0) / swf;
      const qw = vs.reduce((a, v, i) => a + wf[i] * (v.y - yf) ** 2, 0);
      const glw = vs.length - 1;
      i2 = qw > 0 ? Math.max(0, ((qw - glw) / qw) * 100) : 0;
    }

    subgrupos.push({
      etiqueta,
      k: vs.length,
      estimacion: destransformar(comb.y, escala),
      icInferior: destransformar(comb.y - Z95 * comb.se, escala),
      icSuperior: destransformar(comb.y + Z95 * comb.se, escala),
      i2,
    });
  }

  const wg = puntos.map((p) => 1 / (p.se * p.se));
  const sumWg = wg.reduce((a, b) => a + b, 0);
  const ybar = puntos.reduce((a, p, i) => a + wg[i] * p.y, 0) / sumWg;
  const qBet = puntos.reduce((a, p, i) => a + wg[i] * (p.y - ybar) ** 2, 0);
  const gl = puntos.length - 1;
  const pValor = chiSquarePValue(qBet, gl);

  return { campo, subgrupos, qBet, gl, pValor };
}

// --- Meta-regression: relate the effect size to a continuous moderator (weighted least squares) ---

/**
 * Cola normal de dos lados, P(|Z| > z). Exacta (a ~1e-15) vía la χ² con 1 grado de libertad, Z² ~ χ²₁. Antes usaba la
 * aproximación de Abramowitz-Stegun 7.1.26 (error absoluto hasta 1,5e-7), que la validación contra R dejó al descubierto.
 */
export const colaNormalDosLados = (z: number): number => chiSquarePValue(z * z, 1);

// --- Baujat plot: which studies drive the heterogeneity and which drive the pooled result ---

export interface PuntoBaujat {
  id: number;
  etiqueta: string;
  /** Contribución del estudio a la heterogeneidad (Q): w(y−θF)². La suma de todas es Q. */
  x: number;
  /** Influencia del estudio en el resultado combinado (cambio estandarizado al quitarlo). */
  y: number;
}

/**
 * Diagrama de Baujat (Baujat et al. 2002): por cada estudio, X = su contribución a la Q de heterogeneidad
 * (`w(y−θF)²`) e Y = su influencia en el efecto combinado (`(θF − θF₍₋ᵢ₎)² · Σw₍₋ᵢ₎`), en efectos fijos.
 * Los estudios arriba-derecha son los que más conviene revisar. Requiere ≥3 estudios.
 */
export function baujat(resultado: ResultadoMetaAnalisis): PuntoBaujat[] | null {
  const est = resultado.estudios;
  const k = est.length;
  if (k < 3) return null;

  const w = est.map((e) => 1 / (e.se * e.se));
  const sw = w.reduce((a, b) => a + b, 0);
  const swy = est.reduce((acc, e, i) => acc + w[i] * e.y, 0);
  const thetaF = swy / sw;

  return est.map((e, i) => {
    const x = w[i] * (e.y - thetaF) ** 2;
    const swSin = sw - w[i];
    const thetaSin = (swy - w[i] * e.y) / swSin;
    const y = (thetaF - thetaSin) ** 2 * swSin; // (Δθ)² / Var(−i), con Var(−i) = 1/swSin
    return { id: e.id, etiqueta: e.etiqueta, x, y };
  });
}

// --- Diagnóstico de influencia y gráfico de Galbraith ---

/** Un estudio en el diagnóstico de influencia (todo en la escala de análisis: log para razones, etc.). */
export interface DiagnosticoEstudio {
  id: number;
  etiqueta: string;
  /** Residuo estandarizado externamente: (yᵢ − μ̂₍₋ᵢ₎) / √(vᵢ + τ̂²₍₋ᵢ₎ + Var μ̂₍₋ᵢ₎). */
  rstudent: number;
  /** DFFITS: cuántas desviaciones estándar cambia el valor ajustado del estudio al quitarlo. */
  dffits: number;
  /** Distancia de Cook: (μ̂ − μ̂₍₋ᵢ₎)² / Var μ̂. */
  cookD: number;
  /** Razón de covarianzas: Var μ̂₍₋ᵢ₎ / Var μ̂ (< 1: quitarlo daría una estimación más precisa). */
  covR: number;
  /** τ² estimado sin este estudio. */
  tau2Sin: number;
  /** Q de heterogeneidad (pesos de efecto fijo) sin este estudio. */
  qSin: number;
  /** Valor de palanca (hat): peso relativo del estudio en el modelo. */
  hat: number;
  /** Peso del estudio en el modelo, en %. */
  pesoPct: number;
  /** DFBETAS: cuántas desviaciones estándar cambia el efecto combinado al quitarlo. */
  dfbetas: number;
  /** Supera alguno de los umbrales convencionales (marcado con asterisco en metafor). */
  influyente: boolean;
  /** Cuáles: «DFFITS», «Cook», «hat», «DFBETAS». */
  motivos: string[];
}

export interface DiagnosticosInfluencia {
  k: number;
  modelo: "fijo" | "aleatorio";
  metodoTau2: MetodoTau2;
  /** τ² del modelo con todos los estudios (0 en el modelo fijo). */
  tau2: number;
  estudios: DiagnosticoEstudio[];
  /** Umbrales convencionales usados para marcar estudios influyentes (el propio criterio de `metafor::influence`). */
  umbrales: { dffits: number; cook: number; hat: number; dfbetas: number };
}

/** Ajuste de efecto común de un conjunto de puntos: τ² (0 en el modelo fijo), efecto, su varianza y la Q de efecto fijo. */
function ajusteComun(puntos: readonly PuntoEfecto[], modelo: "fijo" | "aleatorio", metodo: MetodoTau2) {
  const w = puntos.map((p) => 1 / (p.se * p.se));
  const sw = w.reduce((a, b) => a + b, 0);
  const muFijo = puntos.reduce((s, p, i) => s + w[i] * p.y, 0) / sw;
  const q = puntos.reduce((s, p, i) => s + w[i] * (p.y - muFijo) ** 2, 0);
  const gl = puntos.length - 1;
  const c = sw - w.reduce((a, b) => a + b * b, 0) / sw;
  const tau2DL = c > 0 ? Math.max(0, (q - gl) / c) : 0;
  const tau2 = modelo === "fijo" ? 0 : estimarTau2(puntos, metodo, tau2DL);
  const { mu, sw: swRE } = ponderarParaTau2(puntos, tau2);
  return { tau2, mu, varMu: 1 / swRE, q };
}

/**
 * **Diagnóstico de influencia** (Viechtbauer y Cheung 2010, como `metafor::influence`): para cada estudio, cuánto cambian el efecto
 * combinado y la heterogeneidad al quitarlo, y cuán atípico es respecto al resto. Todo se calcula en la escala de análisis
 * (log para razones). Con modelo aleatorio el τ² se vuelve a estimar sin cada estudio con el mismo método. Los umbrales para
 * marcar un estudio como influyente son los convencionales del propio metafor (son orientativos, no una prueba): |DFFITS| >
 * 3·√(1/(k−1)), Cook > mediana de χ²₁ (0,455), hat > 3/k, |DFBETAS| > 1. `null` con menos de 3 estudios.
 */
export function diagnosticosInfluencia(
  resultado: ResultadoMetaAnalisis,
  modelo: "fijo" | "aleatorio",
  metodoTau2: MetodoTau2 = "DL",
): DiagnosticosInfluencia | null {
  const puntos: PuntoEfecto[] = resultado.estudios.map((e) => ({ y: e.y, se: e.se }));
  const k = puntos.length;
  if (k < 3) return null;

  const todo = ajusteComun(puntos, modelo, metodoTau2);
  const w = puntos.map((p) => 1 / (p.se * p.se + todo.tau2));
  const sw = w.reduce((a, b) => a + b, 0);
  const umbrales = { dffits: 3 * Math.sqrt(1 / (k - 1)), cook: chiCuadradoCuantil(0.5, 1), hat: 3 / k, dfbetas: 1 };

  const estudios = puntos.map((p, i) => {
    const sin = puntos.filter((_, j) => j !== i);
    const ajuste = ajusteComun(sin, modelo, metodoTau2);
    const v = p.se * p.se;
    const hat = w[i] / sw;
    const delta = todo.mu - ajuste.mu; // cambio del efecto combinado al quitar el estudio
    const rstudent = (p.y - ajuste.mu) / Math.sqrt(v + ajuste.tau2 + ajuste.varMu);
    const dffits = delta / Math.sqrt(hat * (v + ajuste.tau2));
    const cookD = (delta * delta) / todo.varMu;
    const covR = ajuste.varMu / todo.varMu;
    // DFBETAS: el error estándar sale de TODOS los estudios con el τ² estimado sin este (así lo define metafor).
    const seDfb = Math.sqrt(1 / puntos.reduce((s, q) => s + 1 / (q.se * q.se + ajuste.tau2), 0));
    const dfbetas = delta / seDfb;
    const motivos: string[] = [];
    if (Math.abs(dffits) > umbrales.dffits) motivos.push("DFFITS");
    if (cookD > umbrales.cook) motivos.push("Cook");
    if (hat > umbrales.hat) motivos.push("hat");
    if (Math.abs(dfbetas) > umbrales.dfbetas) motivos.push("DFBETAS");
    return {
      id: resultado.estudios[i].id,
      etiqueta: resultado.estudios[i].etiqueta,
      rstudent,
      dffits,
      cookD,
      covR,
      tau2Sin: ajuste.tau2,
      qSin: ajuste.q,
      hat,
      pesoPct: hat * 100,
      dfbetas,
      influyente: motivos.length > 0,
      motivos,
    };
  });
  return { k, modelo, metodoTau2, tau2: todo.tau2, estudios, umbrales };
}

export interface PuntoGalbraith {
  id: number;
  etiqueta: string;
  /** Precisión: 1/√(vᵢ + τ²). */
  x: number;
  /** Efecto estandarizado: yᵢ/√(vᵢ + τ²). */
  z: number;
}

export interface GraficoGalbraith {
  puntos: PuntoGalbraith[];
  /** Pendiente de la recta por el origen = efecto combinado en la escala de análisis. */
  pendiente: number;
  tau2: number;
  /** Semiancho de la banda (1,96): los puntos fuera de ella se apartan del efecto combinado más de lo esperable por azar. */
  banda: number;
}

/**
 * **Gráfico de Galbraith (radial)**: cada estudio se dibuja como (precisión, efecto estandarizado) = (1/√(vᵢ+τ²), yᵢ/√(vᵢ+τ²)).
 * Con el efecto común μ̂ como pendiente, un estudio compatible con él cae a menos de 1,96 unidades verticales de la recta por el
 * origen. Igual que `metafor::radial`: con el modelo aleatorio se usa el τ² del ajuste. `null` con menos de 2 estudios.
 */
export function galbraith(resultado: ResultadoMetaAnalisis, modelo: "fijo" | "aleatorio"): GraficoGalbraith | null {
  const est = resultado.estudios;
  if (est.length < 2) return null;
  const tau2 = modelo === "aleatorio" ? resultado.heterogeneidad.tau2 : 0;
  const pesos = est.map((e) => 1 / (e.se * e.se + tau2));
  const sw = pesos.reduce((a, b) => a + b, 0);
  const pendiente = est.reduce((s, e, i) => s + pesos[i] * e.y, 0) / sw;
  return {
    puntos: est.map((e, i) => ({ id: e.id, etiqueta: e.etiqueta, x: Math.sqrt(pesos[i]), z: e.y * Math.sqrt(pesos[i]) })),
    pendiente,
    tau2,
    banda: Z95,
  };
}

// --- GRADE certainty of evidence (Summary of Findings) ---

/** Señales derivadas del meta-análisis que alimentan los dominios GRADE de datos. */
export interface SeñalesGrade {
  i2: number;
  /** El IC combinado cruza el valor nulo (imprecisión). */
  cruzaNulo: boolean;
  /** false para medidas sin nulo (proporciones): la imprecisión no se juzga automáticamente. */
  imprecisionEvaluable: boolean;
  /** true = sesgo detectado; false = no; null = no evaluable (<10 estudios) o ambiguo (ROB-ME "Algunas dudas"). */
  sesgoPublicacion: boolean | null;
  /** Motivo textual del veredicto de `sesgoPublicacion` (típicamente de `sesgoPublicacionParaGrade`,
   *  que prioriza ROB-ME sobre Egger/Peters). Si se omite, se usa un texto genérico. */
  sesgoPublicacionMotivo?: string;
  /** true si `sesgoPublicacion` es null por AMBIGÜEDAD real (ROB-ME "Algunas dudas"), no por falta de
   *  datos — la UI debe seguir marcando el dominio para revisión manual aunque no reste niveles. */
  sesgoPublicacionAdvertencia?: boolean;
  /**
   * % del peso total (0-100) del modelo elegido aportado por estudios de riesgo Alto o Poco Claro.
   * Reemplaza el conteo simple de estudios: un estudio pequeño de alto riesgo debe pesar poco en el
   * juicio, uno grande debe pesar mucho — igual que pesa en el meta-análisis mismo.
   */
  pesoAltoRiesgo: number;
  /**
   * % del peso total (0-100) aportado por estudios SIN riesgo de sesgo evaluado con una herramienta (Paso 6).
   * 100 = ningún estudio evaluado: el dominio no se puede juzgar y GRADE debe decirlo, no darlo por «sin descenso».
   * Si se omite se asume 0 (todos evaluados).
   */
  pesoSinEvaluar?: number;
  /** Estimación combinada en escala natural (destransformada), para el factor de aumento "efecto grande". */
  estimacionNatural: number;
  /** Valor nulo en escala natural (1 para razones, 0 para diferencias). null si no aplica (p.ej. Proporción). */
  nuloNatural: number | null;
  /** Límites del IC combinado en escala natural, para el chequeo de amplitud (imprecisión). */
  icInferiorNatural: number;
  icSuperiorNatural: number;
}

export interface DominioGrade {
  nombre: string;
  /** Niveles que baja este dominio (0, 1 o 2). */
  baja: number;
  motivo: string;
  /** true = este dominio no bajó automáticamente (baja=0) pero sigue pendiente de revisión manual
   *  por ambigüedad real (no por falta de datos) — la UI debe mostrarlo igual que uno que sí bajó. */
  advertencia?: boolean;
}

export interface ResultadoGrade {
  nivelInicial: CertezaGrade;
  nivelFinal: CertezaGrade;
  dominios: DominioGrade[];
  nEstudios: number;
  participantes: number | null;
  /** true si NINGÚN estudio del resultado tiene el riesgo de sesgo evaluado: la certeza mostrada no descuenta ese dominio. */
  riesgoSesgoNoEvaluado: boolean;
}

const NIVELES_GRADE: Record<number, CertezaGrade> = { 4: "Alta", 3: "Moderada", 2: "Baja", 1: "Muy Baja" };
const nivelDesde = (n: number): CertezaGrade => NIVELES_GRADE[Math.max(1, Math.min(4, n))];

/**
 * Sugiere la certeza GRADE de un desenlace: arranca en Alta (mayoría ECA) o Baja (observacional) y
 * ajusta por los 5 dominios de descenso (riesgo de sesgo, inconsistencia, imprecisión, evidencia
 * indirecta, sesgo de publicación) y, si arrancó en Baja, los factores de aumento (efecto grande
 * automático; dosis-respuesta y confusión residual quedan siempre a tu juicio). Es una sugerencia
 * editable, con el motivo de cada dominio — no un veredicto final.
 *
 * @note Fuente: GRADE Handbook (Guyatt et al., BMJ 2008-2011 series). Validación conceptual:
 *       - Riesgo de sesgo pondera por el % del peso del meta-análisis (no cuenta estudios: un
 *         estudio pequeño de alto riesgo no debería pesar igual que uno grande de bajo riesgo).
 *       - Imprecisión baja tanto si el IC cruza el nulo como si es muy ancho sin cruzarlo (un IC
 *         de [1.05, 9.8] no cruza el nulo pero es demasiado amplio para ser "preciso").
 *       - El aumento por "efecto grande" (RR/OR ≥2 o ≤0.5 → +1; ≥5 o ≤0.2 → +2) solo aplica si el
 *         nivel de partida es Baja (evidencia mayormente observacional) y requiere que NO haya
 *         confusores residuales plausibles — se muestra el candidato pero el aumento queda anotado
 *         para que lo confirmes, igual que Cochrane pide para este factor (es el más disputado de los tres).
 *       - Sesgo de publicación: NO se calcula aquí — el llamador arma `sesgoPublicacion` con
 *         `sesgoPublicacionParaGrade` (rob-me.ts), que prioriza el veredicto ROB-ME de ese resultado
 *         sobre el proxy estadístico (Egger/Peters) cuando ya está completo.
 */
export function evaluarGrade(estudios: Estudio[], medida: Medida, señales: SeñalesGrade): ResultadoGrade | null {
  const delMedida = estudios.filter((e) => e.medida === medida);
  const nEstudios = delMedida.length;
  if (nEstudios === 0) return null;

  let participantes: number | null = null;
  for (const e of delMedida) {
    const n = parseInt(e.n, 10);
    if (Number.isFinite(n)) participantes = (participantes ?? 0) + n;
  }

  const nEca = delMedida.filter((e) => e.diseño === "ECA").length;
  const mayoriaEca = nEca * 2 >= nEstudios;
  let nivel = mayoriaEca ? 4 : 2;
  const nivelInicial = nivelDesde(nivel);
  const dominios: DominioGrade[] = [];

  // 1. Riesgo de sesgo — ponderado por el peso real en el meta-análisis, no por conteo de estudios.
  //    Si nadie lo evaluó, el dominio NO se da por «sin descenso»: se marca como pendiente y la certeza lo dice.
  const pesoSinEvaluar = señales.pesoSinEvaluar ?? 0;
  const riesgoSesgoNoEvaluado = pesoSinEvaluar >= 99.95;
  if (riesgoSesgoNoEvaluado) {
    dominios.push({
      nombre: "Riesgo de sesgo",
      baja: 0,
      advertencia: true,
      motivo:
        "Riesgo de sesgo no evaluado: ningún estudio de este resultado tiene la herramienta del Paso 6 completada, así que este dominio no se pudo juzgar. La certeza mostrada no lo descuenta y podría ser menor.",
    });
  } else if (señales.pesoAltoRiesgo > 50) {
    nivel -= 1;
    dominios.push({
      nombre: "Riesgo de sesgo",
      baja: 1,
      motivo: `El ${señales.pesoAltoRiesgo.toFixed(0)}% del peso del meta-análisis proviene de estudios con riesgo de sesgo alto o poco claro.`,
    });
  } else if (señales.pesoAltoRiesgo > 25) {
    dominios.push({
      nombre: "Riesgo de sesgo",
      baja: 0,
      motivo: `${señales.pesoAltoRiesgo.toFixed(0)}% del peso viene de estudios de riesgo alto/poco claro — no domina, pero revísalo.`,
    });
  } else if (pesoSinEvaluar > 25) {
    dominios.push({
      nombre: "Riesgo de sesgo",
      baja: 0,
      advertencia: true,
      motivo: `${pesoSinEvaluar.toFixed(0)}% del peso viene de estudios sin el riesgo de sesgo evaluado (Paso 6): el dominio solo se pudo juzgar con el resto.`,
    });
  } else {
    dominios.push({
      nombre: "Riesgo de sesgo",
      baja: 0,
      motivo: "La mayor parte del peso del meta-análisis proviene de estudios con bajo riesgo de sesgo.",
    });
  }

  // 2. Inconsistencia.
  if (señales.i2 > 50) {
    nivel -= 1;
    dominios.push({ nombre: "Inconsistencia", baja: 1, motivo: `I² = ${señales.i2.toFixed(0)}%: heterogeneidad importante entre estudios.` });
  } else {
    dominios.push({ nombre: "Inconsistencia", baja: 0, motivo: `I² = ${señales.i2.toFixed(0)}%: heterogeneidad baja o moderada.` });
  }

  // 3. Imprecisión — cruza el nulo, O el IC es muy ancho aunque no lo cruce (proxy de amplitud,
  //    no sustituye el juicio con un umbral clínico real, pero evita marcar como "preciso" un
  //    intervalo como OR 3.2 [1.05–9.8] solo porque técnicamente no toca 1).
  if (!señales.imprecisionEvaluable) {
    dominios.push({ nombre: "Imprecisión", baja: 0, motivo: "No se evalúa automáticamente para esta medida; revísala a mano." });
  } else if (señales.cruzaNulo) {
    nivel -= 1;
    dominios.push({ nombre: "Imprecisión", baja: 1, motivo: "El intervalo de confianza cruza el valor nulo: el resultado no es concluyente." });
  } else {
    const anchoAmplio = esIntervaloAmplio(señales);
    if (anchoAmplio) {
      nivel -= 1;
      dominios.push({
        nombre: "Imprecisión",
        baja: 1,
        motivo: `El intervalo (${señales.icInferiorNatural.toFixed(2)}–${señales.icSuperiorNatural.toFixed(2)}) no cruza el nulo, pero es muy amplio: cubre un rango demasiado grande de efectos posibles para considerarlo preciso.`,
      });
    } else {
      dominios.push({ nombre: "Imprecisión", baja: 0, motivo: "El intervalo de confianza no cruza el nulo y su amplitud es razonable." });
    }
  }

  // 4. Evidencia indirecta — no automatizable (depende de si P/I/C/O responden tu pregunta real).
  dominios.push({
    nombre: "Evidencia indirecta",
    baja: 0,
    motivo: "Evalúa a mano si población, intervención y desenlace responden tu pregunta PICOT (no automatizable).",
  });

  // 5. Sesgo de publicación — prioriza ROB-ME sobre Egger/Peters cuando está disponible, ver
  //    `sesgoPublicacionParaGrade` en rob-me.ts (quien arma `sesgoPublicacion`/`sesgoPublicacionMotivo`).
  if (señales.sesgoPublicacion === true) {
    nivel -= 1;
    dominios.push({
      nombre: "Sesgo de publicación",
      baja: 1,
      motivo: señales.sesgoPublicacionMotivo ?? "Los tests (Egger/Peters) sugieren asimetría del funnel.",
    });
  } else if (señales.sesgoPublicacion === false) {
    dominios.push({
      nombre: "Sesgo de publicación",
      baja: 0,
      motivo: señales.sesgoPublicacionMotivo ?? "Sin evidencia estadística de sesgo de publicación.",
    });
  } else {
    dominios.push({
      nombre: "Sesgo de publicación",
      baja: 0,
      advertencia: señales.sesgoPublicacionAdvertencia,
      motivo: señales.sesgoPublicacionMotivo ?? "No evaluable de forma fiable con menos de 10 estudios.",
    });
  }

  // Factores de AUMENTO — solo tienen sentido cuando se arrancó en Baja (evidencia mayormente
  // observacional). GRADE define 3: efecto grande (automatizable con cautela), gradiente
  // dosis-respuesta y confusión residual que reduciría el efecto (ambos requieren datos que esta
  // herramienta no tiene — se muestran como recordatorio, nunca se aplican solos).
  if (!mayoriaEca && señales.nuloNatural !== null) {
    const razon = señales.nuloNatural !== 0 ? señales.estimacionNatural / señales.nuloNatural : null;
    if (razon !== null && Number.isFinite(razon)) {
      const magnitud = razon >= 1 ? razon : 1 / razon;
      if (magnitud >= 5) {
        nivel += 2;
        dominios.push({
          nombre: "Efecto muy grande (aumenta)",
          baja: -2,
          motivo: `El efecto combinado es muy grande (≈${magnitud.toFixed(1)}× el nulo). GRADE permite subir 2 niveles SI no hay confusores residuales plausibles que expliquen esta magnitud — confírmalo con tu juicio antes de aceptar el aumento.`,
        });
      } else if (magnitud >= 2) {
        nivel += 1;
        dominios.push({
          nombre: "Efecto grande (aumenta)",
          baja: -1,
          motivo: `El efecto combinado es grande (≈${magnitud.toFixed(1)}× el nulo). GRADE permite subir 1 nivel SI no hay confusores residuales plausibles que expliquen esta magnitud — confírmalo con tu juicio antes de aceptar el aumento.`,
        });
      }
    }
    dominios.push({
      nombre: "Gradiente dosis-respuesta (revisa a mano)",
      baja: 0,
      motivo: "Si a mayor dosis/exposición se observa mayor efecto de forma consistente, GRADE permite subir 1 nivel. No se detecta automáticamente.",
    });
    dominios.push({
      nombre: "Confusión residual (revisa a mano)",
      baja: 0,
      motivo: "Si TODOS los confusores plausibles no controlados reducirían (no exagerarían) el efecto observado, GRADE permite subir 1 nivel. Requiere tu juicio clínico, no es automatizable.",
    });
  }

  return { nivelInicial, nivelFinal: nivelDesde(nivel), dominios, nEstudios, participantes, riesgoSesgoNoEvaluado };
}

/**
 * Heurística de "intervalo demasiado amplio" para el dominio de imprecisión, cuando el IC no cruza
 * el nulo pero aun así cubre un rango muy grande de efectos plausibles.
 * - Escala log (OR/RR/HR): razón superior/inferior ≥ 2.5 (p.ej. algo entre "poco efecto" y "efecto enorme").
 * - Escala lineal (DM/DME/RAR): semi-ancho del IC mayor que la propia estimación (la incertidumbre
 *   es tan grande como el efecto mismo).
 * No sustituye un umbral clínico real (GRADE lo pide basado en decisiones clínicas concretas); es un
 * proxy conservador para no reportar como "preciso" lo que es evidentemente inestable.
 */
function esIntervaloAmplio(señales: SeñalesGrade): boolean {
  const { icInferiorNatural, icSuperiorNatural, estimacionNatural, nuloNatural } = señales;
  if (nuloNatural === 1) {
    // Escala de razón: ambos límites deben ser positivos para que la razón tenga sentido.
    if (icInferiorNatural <= 0) return true; // límite inválido/degenerado: trátalo como impreciso.
    return icSuperiorNatural / icInferiorNatural >= 2.5;
  }
  // Escala lineal (nulo = 0): compara el semi-ancho contra la magnitud del efecto.
  const semiAncho = (icSuperiorNatural - icInferiorNatural) / 2;
  return semiAncho > Math.abs(estimacionNatural);
}

/** Efecto absoluto de un desenlace binario: riesgo basal, riesgo con la intervención, y la diferencia. */
export interface EfectoAbsoluto {
  /** Riesgo del grupo control combinado (proporción 0–1), estimado desde los datos crudos disponibles. */
  riesgoBasal: number;
  /** Riesgo estimado con la intervención, aplicando el efecto relativo combinado al riesgo basal. */
  riesgoIntervencion: number;
  /** (riesgoIntervencion − riesgoBasal) × 1000 — convención GRADE de reportar "por 1,000 pacientes". */
  diferenciaPor1000: number;
  /** Nº de estudios cuyos datos crudos 2×2 se usaron para estimar el riesgo basal. */
  nEstudiosConDatos: number;
}

/**
 * Convierte un efecto relativo (OR o RR) en efecto absoluto — lo que GRADE llama "anticipated absolute
 * effects" y que Cochrane Handbook exige en toda tabla Summary of Findings. Un OR de 0.66 no dice si el
 * beneficio es grande o trivial; "80 eventos menos por cada 1,000 pacientes" sí.
 *
 * @param estudios - Estudios de la medida (se usan los que tengan datos crudos 2×2 cargados en el Paso 5)
 * @param medida - Solo "OR" o "RR" tienen una conversión estándar simple; otras medidas retornan null
 * @param estimacionPooled - El efecto combinado (OR o RR) en escala natural, del modelo elegido
 * @returns EfectoAbsoluto, o null si la medida no aplica o no hay datos crudos para estimar el riesgo basal
 *
 * @note El riesgo basal se estima como Σeventos_control / Σtotal_control entre los estudios CON datos
 *       crudos — es una aproximación (no pondera por el peso del meta-análisis), consistente con cómo
 *       GRADEpro y las tablas SoF publicadas suelen presentar un "riesgo basal ilustrativo".
 */
export function calcularEfectoAbsoluto(
  estudios: Estudio[],
  medida: Medida,
  estimacionPooled: number,
): EfectoAbsoluto | null {
  if (medida !== "OR" && medida !== "RR") return null;

  const conCrudos = estudios.filter(
    (e): e is Estudio & { crudos: Extract<DatosCrudos, { tipo: "2x2" }> } =>
      e.medida === medida &&
      e.crudos?.tipo === "2x2" &&
      Number.isFinite(parseFloat(e.crudos.eventosCtrl)) &&
      Number.isFinite(parseFloat(e.crudos.totalCtrl)) &&
      parseFloat(e.crudos.totalCtrl) > 0,
  );
  if (conCrudos.length === 0) return null;

  let eventos = 0;
  let total = 0;
  for (const e of conCrudos) {
    const kCtrl = divisorControl(e.dividirControlEn); // un control compartido no se suma dos veces al riesgo basal
    eventos += parseFloat(e.crudos.eventosCtrl) / kCtrl;
    total += parseFloat(e.crudos.totalCtrl) / kCtrl;
  }
  if (total <= 0) return null;
  const riesgoBasal = eventos / total;

  let riesgoIntervencion: number;
  if (medida === "RR") {
    riesgoIntervencion = riesgoBasal * estimacionPooled;
  } else {
    // OR: convertir a momios, aplicar el OR, convertir de vuelta a riesgo.
    const oddsBasal = riesgoBasal / (1 - riesgoBasal);
    const oddsIntervencion = oddsBasal * estimacionPooled;
    riesgoIntervencion = oddsIntervencion / (1 + oddsIntervencion);
  }
  riesgoIntervencion = Math.max(0, Math.min(1, riesgoIntervencion));

  return {
    riesgoBasal,
    riesgoIntervencion,
    diferenciaPor1000: (riesgoIntervencion - riesgoBasal) * 1000,
    nEstudiosConDatos: conCrudos.length,
  };
}

/** Valor nulo (sin efecto) de una escala: 1 para la escala log (OR, RR, HR), 0 para la lineal y la z de Fisher; `null` si la escala no tiene un nulo con sentido (p. ej. proporciones). */
export function valorNulo(escala: EscalaAnalisis): number | null {
  if (escala === "log") return 1;
  if (escala === "lineal" || escala === "fisher") return 0;
  return null; // logit / proportions have no null-effect line
}

/** ¿El intervalo de confianza excluye el valor nulo? (es decir, ¿el efecto combinado es estadísticamente significativo?). */
export function esSignificativo(icInferior: number, icSuperior: number, nulo: number | null): boolean {
  if (nulo === null) return false;
  return nulo < icInferior || nulo > icSuperior;
}

/**
 * Validación en vivo de una fila de la tabla de extracción: devuelve `null` mientras la fila esté sin tocar (para no
 * molestar en una fila nueva en blanco) y, si no, un mensaje de error específico.
 */
export function validarEfectoIC(estudio: Pick<Estudio, "medida" | "valor" | "icInferior" | "icSuperior">): string | null {
  const medida = estudio.medida;
  const valor = estudio.valor ?? "";
  const icInferior = estudio.icInferior ?? "";
  const icSuperior = estudio.icSuperior ?? "";
  if (!valor.trim() && !icInferior.trim() && !icSuperior.trim()) return null;

  const escala = escalaPorMedida[medida];
  const est = parseFloat(valor);
  const inf = parseFloat(icInferior);
  const sup = parseFloat(icSuperior);

  if (valor.trim() && !Number.isFinite(est)) return "La estimación debe ser un número.";
  if (icInferior.trim() && !Number.isFinite(inf)) return "El IC inferior debe ser un número.";
  if (icSuperior.trim() && !Number.isFinite(sup)) return "El IC superior debe ser un número.";
  if (Number.isFinite(inf) && Number.isFinite(sup) && inf >= sup) {
    return "El IC inferior debe ser menor que el superior.";
  }
  if (
    escala &&
    ((Number.isFinite(est) && !dentroDeRango(est, escala)) ||
      (Number.isFinite(inf) && !dentroDeRango(inf, escala)) ||
      (Number.isFinite(sup) && !dentroDeRango(sup, escala)))
  ) {
    return escala === "logit"
      ? `Para ${medida}, los valores deben estar entre 0 y 1 (ej. 0.35, no 35).`
      : escala === "fisher"
        ? `Para ${medida}, los valores deben estar entre −1 y 1 (un coeficiente de correlación).`
        : `Para ${medida} (escala logarítmica), los valores deben ser mayores que cero.`;
  }
  return null;
}

// --- SROC: precisión diagnóstica (sensibilidad/especificidad combinadas), Fase 4 ---

export interface PuntoSROC {
  id: number;
  etiqueta: string;
  vp: number;
  fp: number;
  fn: number;
  vn: number;
  n: number;
  sensibilidad: number;
  sensIcInferior: number;
  sensIcSuperior: number;
  especificidad: number;
  espIcInferior: number;
  espIcSuperior: number;
  /** Tasa de falsos positivos (1−especificidad); eje X del espacio SROC. */
  fpr: number;
  correccion: boolean;
}

export interface RegresionSROC {
  /** Intercepto y pendiente de D = a + b·S (Moses-Littenberg, WLS con peso 1/varianza(D)). */
  a: number;
  b: number;
  seA: number;
  seB: number;
  /** Test de la pendiente (Wald, normal): b≠0 sugiere un efecto de umbral (SROC asimétrica). */
  zB: number;
  pValorB: number;
}

export interface ResultadoSROC {
  puntos: PuntoSROC[];
  /** null si hay <4 estudios (la regresión de 2 parámetros necesita más grados de libertad que eso para ser estable). */
  regresion: RegresionSROC | null;
  curva: { fpr: number; tpr: number }[] | null;
  /** Odds ratio diagnóstico combinado (misma tabla 2×2 que la regresión, agrupado por inverse-variance). */
  dor: {
    estimacion: number;
    icInferior: number;
    icSuperior: number;
    q: number;
    gl: number;
    pValor: number;
    i2: number;
    tau2: number;
  };
  /** Sensibilidad/especificidad combinadas por separado (DerSimonian-Laird en escala logit) — una
   *  simplificación reconocida (Cochrane Handbook, cap. diagnóstico) frente al modelo bivariado
   *  completo (Reitsma/HSROC), que no está implementado: ignora la correlación entre ambas. */
  sensibilidadPromedio: { estimacion: number; icInferior: number; icSuperior: number };
  especificidadPromedio: { estimacion: number; icInferior: number; icSuperior: number };
  modelo: "fijo" | "aleatorio";
  metodoTau2: MetodoTau2;
  /**
   * **Modelo bivariado** (Reitsma/HSROC; ver `sroc-bivariado.ts`): el resumen recomendado por el Cochrane Handbook de precisión
   * diagnóstica. Verosimilitud binomial exacta (o, si no se pudo ajustar, la aproximación normal de Reitsma: ver su `metodo`). `null`
   * con menos de 4 estudios o si se pidió no calcularlo. La sensibilidad y la especificidad de arriba, agrupadas POR SEPARADO, quedan
   * como análisis exploratorio: el manual advierte de que ignoran la correlación entre ambas.
   */
  bivariado: ResultadoBivariado | null;
  /** La aproximación normal de Reitsma (REML, como `mada::reitsma`), para comparar con `bivariado`; `null` si `bivariado` ya es ésa. */
  bivariadoNormal: ResultadoBivariado | null;
}

/** Agrupa puntos {y,se} por inverse-variance (fijo o DerSimonian-Laird/PM/REML aleatorio) con heterogeneidad completa. */
function agruparConHeterogeneidad(
  pts: readonly PuntoEfecto[],
  modelo: "fijo" | "aleatorio",
  metodoTau2: MetodoTau2,
): { y: number; se: number; q: number; gl: number; pValor: number; i2: number; tau2: number } {
  const wf = pts.map((p) => 1 / (p.se * p.se));
  const swf = wf.reduce((a, b) => a + b, 0);
  const yF = pts.reduce((acc, p, i) => acc + wf[i] * p.y, 0) / swf;
  const gl = pts.length - 1;
  const q = pts.reduce((acc, p, i) => acc + wf[i] * (p.y - yF) ** 2, 0);
  const sw2 = wf.reduce((a, b) => a + b * b, 0);
  const c = swf - sw2 / swf;
  const tau2DL = c > 0 ? Math.max(0, (q - gl) / c) : 0;
  const tau2 = modelo === "aleatorio" ? estimarTau2(pts, metodoTau2, tau2DL) : 0;
  const i2 = q > 0 ? Math.max(0, ((q - gl) / q) * 100) : 0;
  const pValor = chiSquarePValue(q, gl);
  if (modelo === "fijo") return { y: yF, se: Math.sqrt(1 / swf), q, gl, pValor, i2, tau2: 0 };
  const wr = pts.map((p) => 1 / (p.se * p.se + tau2));
  const swr = wr.reduce((a, b) => a + b, 0);
  const yR = pts.reduce((acc, p, i) => acc + wr[i] * p.y, 0) / swr;
  return { y: yR, se: Math.sqrt(1 / swr), q, gl, pValor, i2, tau2 };
}

/**
 * **SROC (Summary Receiver Operating Characteristic)** — precisión diagnóstica combinada. Lee
 * `Estudio.datosDiagnosticos` (VP/FP/FN/VN), independiente de `medida`/`valor`.
 *
 * - Por estudio: sensibilidad/especificidad + IC 95% (Wilson).
 * - Regresión de Moses-Littenberg (Moses, Shapiro & Littenberg 1993): `D=ln(DOR)` sobre `S` (WLS,
 *   peso = 1/varianza(D)); prueba si la pendiente ≠0 (asimetría / efecto de umbral) y traza la curva
 *   SROC cuando |b|<1.
 * - DOR combinado y sensibilidad/especificidad combinadas por separado (inverse-variance, mismo
 *   `modelo`/`metodoTau2` elegidos en el Paso 7 para el resto del análisis).
 *
 * @note Pooling univariado de sensibilidad/especificidad (no el modelo bivariado de Reitsma/HSROC) —
 *       ver el campo `sensibilidadPromedio`/`especificidadPromedio`. Alcance declarado en la UI.
 * @returns null si hay menos de 3 estudios con datos diagnósticos completos y válidos.
 */
export function calcularSROC(
  estudios: Estudio[],
  modelo: "fijo" | "aleatorio" = "aleatorio",
  metodoTau2: MetodoTau2 = "DL",
  opciones: { bivariado?: boolean } = {},
): ResultadoSROC | null {
  const puntos: PuntoSROC[] = [];
  const diag: PuntoDiagnostico[] = [];
  estudios.forEach((e) => {
    const dd = e.datosDiagnosticos;
    if (!dd) return;
    const vp = parseFloat(dd.vp);
    const fp = parseFloat(dd.fp);
    const fn = parseFloat(dd.fn);
    const vn = parseFloat(dd.vn);
    const p = puntoDiagnostico(vp, fp, fn, vn);
    if (!p) return;
    const sens = intervaloWilson(vp, vp + fn);
    const esp = intervaloWilson(vn, fp + vn);
    puntos.push({
      id: e.id,
      etiqueta: `${e.autor || "Estudio"}${e.año ? ` (${e.año})` : ""}`,
      vp,
      fp,
      fn,
      vn,
      n: vp + fp + fn + vn,
      sensibilidad: sens.estimacion,
      sensIcInferior: sens.icInferior,
      sensIcSuperior: sens.icSuperior,
      especificidad: esp.estimacion,
      espIcInferior: esp.icInferior,
      espIcSuperior: esp.icSuperior,
      fpr: 1 - esp.estimacion,
      correccion: p.correccion,
    });
    diag.push(p);
  });
  if (puntos.length < 3) return null;

  const dorPts: PuntoEfecto[] = diag.map((p) => ({ y: p.dLn, se: p.seD }));
  const dorAgrupado = agruparConHeterogeneidad(dorPts, modelo, metodoTau2);
  const dor = {
    estimacion: Math.exp(dorAgrupado.y),
    icInferior: Math.exp(dorAgrupado.y - Z95 * dorAgrupado.se),
    icSuperior: Math.exp(dorAgrupado.y + Z95 * dorAgrupado.se),
    q: dorAgrupado.q,
    gl: dorAgrupado.gl,
    pValor: dorAgrupado.pValor,
    i2: dorAgrupado.i2,
    tau2: dorAgrupado.tau2,
  };

  const sensPts: PuntoEfecto[] = diag.map((p) => ({ y: Math.log(p.a / p.b), se: Math.sqrt(1 / p.a + 1 / p.b) }));
  const espPts: PuntoEfecto[] = diag.map((p) => ({ y: Math.log(p.d / p.c), se: Math.sqrt(1 / p.c + 1 / p.d) }));
  const sensAgrupado = agruparConHeterogeneidad(sensPts, modelo, metodoTau2);
  const espAgrupado = agruparConHeterogeneidad(espPts, modelo, metodoTau2);
  const sensibilidadPromedio = {
    estimacion: destransformar(sensAgrupado.y, "logit"),
    icInferior: destransformar(sensAgrupado.y - Z95 * sensAgrupado.se, "logit"),
    icSuperior: destransformar(sensAgrupado.y + Z95 * sensAgrupado.se, "logit"),
  };
  const especificidadPromedio = {
    estimacion: destransformar(espAgrupado.y, "logit"),
    icInferior: destransformar(espAgrupado.y - Z95 * espAgrupado.se, "logit"),
    icSuperior: destransformar(espAgrupado.y + Z95 * espAgrupado.se, "logit"),
  };

  let regresion: RegresionSROC | null = null;
  let curva: { fpr: number; tpr: number }[] | null = null;
  if (diag.length >= 4) {
    const r = regresionMosesLittenberg(diag);
    if (r) {
      regresion = { a: r.a, b: r.b, seA: r.seA, seB: r.seB, zB: r.zB, pValorB: r.seB > 0 ? Math.min(1, colaNormalDosLados(r.zB)) : 1 };
      curva = curvaMosesLittenberg(regresion);
    }
  }

  // Modelo bivariado (no depende de `modelo` ni de `metodoTau2`: es siempre de efectos aleatorios, ML o REML).
  let bivariado: ResultadoBivariado | null = null;
  let bivariadoNormal: ResultadoBivariado | null = null;
  if (opciones.bivariado !== false && puntos.length >= 4) {
    const tablas = puntos.map((p) => ({ vp: p.vp, fp: p.fp, fn: p.fn, vn: p.vn }));
    bivariado = ajustarBivariadoConCache(tablas, "binomial");
    if (bivariado && bivariado.metodo === "binomial") bivariadoNormal = ajustarBivariadoConCache(tablas, "normal");
  }

  return { puntos, regresion, curva, dor, sensibilidadPromedio, especificidadPromedio, modelo, metodoTau2, bivariado, bivariadoNormal };
}
