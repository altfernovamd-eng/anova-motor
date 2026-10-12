// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * **Regresión de Moses-Littenberg y piezas por estudio de la curva SROC univariada** (precisión diagnóstica), escritas desde los artículos
 * (reescrito el 7 de octubre de 2026; ver docs/PROCEDENCIA.md): la versión anterior se había escrito mientras se leía `mada::madauni`
 * (GPL); esta parte de las ecuaciones publicadas y usa `mada` solo como caja negra de comparación.
 *
 * Por estudio, con la tabla 2×2 (a = VP, b = FN, c = FP, d = VN), y +0,5 en las cuatro celdas si alguna es cero (Handbook 10.4.4.1):
 *  · D = ln(DOR) = logit(sensibilidad) − logit(FPR) = ln(a·d / b·c), con varianza 1/a + 1/b + 1/c + 1/d (Glas y cols. 2003).
 *  · S = logit(sensibilidad) + logit(FPR) = ln(a/b) + ln(c/d): un indicador del umbral de positividad (Moses, Shapiro y Littenberg 1993).
 *  · Regresión D = α + β·S por mínimos cuadrados ponderados con peso 1/Var(D) (Littenberg y Moses 1993). Si β ≠ 0, el DOR cambia con el
 *    umbral y la curva es asimétrica. Los errores estándar salen de (XᵀWX)⁻¹ con los pesos conocidos, y la prueba de β es de Wald (normal).
 *  · Curva SROC: de D = α + β·S, con D = L_T − L_F y S = L_T + L_F (L = logit): L_T = α/(1 − β) + ((1 + β)/(1 − β))·L_F, para |β| < 1.
 *
 * Sin dependencias de la app. Lógica pura.
 */

const Z95 = 1.959963985;

/** Intervalo de Wilson de una proporción x/n (más estable que el de Wald cuando x = 0 o x = n). Wilson 1927; Newcombe 1998. */
export function intervaloWilson(x: number, n: number): { estimacion: number; icInferior: number; icSuperior: number } {
  const p = x / n;
  const z2 = Z95 * Z95;
  const denominador = 1 + z2 / n;
  const centro = (p + z2 / (2 * n)) / denominador;
  const radio = (Z95 / denominador) * Math.sqrt((p * (1 - p)) / n + z2 / (4 * n * n));
  return { estimacion: p, icInferior: Math.max(0, centro - radio), icSuperior: Math.min(1, centro + radio) };
}

export interface PuntoDiagnostico {
  /** Celdas con la corrección de continuidad aplicada: a = VP, b = FN, c = FP, d = VN. */
  a: number;
  b: number;
  c: number;
  d: number;
  /** true si se sumó 0,5 a las cuatro celdas por tener alguna en cero. */
  correccion: boolean;
  /** D = ln(DOR) y su error estándar. */
  dLn: number;
  seD: number;
  /** S = logit(sensibilidad) + logit(FPR). */
  s: number;
}

/** Tabla de un estudio diagnóstico (VP, FP, FN, VN) → ln DOR con su error estándar y S = logit(sensibilidad) + logit(FPR), las dos variables de la regresión de Moses-Littenberg; suma 0,5 a las cuatro celdas si alguna es cero. `null` si los conteos no son válidos (negativos, no finitos, o sin enfermos o sin sanos). */
export function puntoDiagnostico(vp: number, fp: number, fn: number, vn: number): PuntoDiagnostico | null {
  if (![vp, fp, fn, vn].every((x) => Number.isFinite(x) && x >= 0)) return null;
  if (vp + fn <= 0 || fp + vn <= 0) return null;
  const hayCero = vp === 0 || fn === 0 || fp === 0 || vn === 0;
  const k = hayCero ? 0.5 : 0;
  const [a, b, c, d] = [vp + k, fn + k, fp + k, vn + k];
  const dLn = Math.log(a) + Math.log(d) - Math.log(b) - Math.log(c);
  const seD = Math.sqrt(1 / a + 1 / b + 1 / c + 1 / d);
  const s = Math.log(a / b) + Math.log(c / d);
  if (![dLn, seD, s].every(Number.isFinite) || seD <= 0) return null;
  return { a, b, c, d, correccion: hayCero, dLn, seD, s };
}

export interface RegresionMosesLittenberg {
  /** Intercepto y pendiente de D = α + β·S. */
  a: number;
  b: number;
  seA: number;
  seB: number;
  /** Prueba de Wald de la pendiente (normal) y su valor p de dos lados. */
  zB: number;
}

/**
 * Mínimos cuadrados ponderados de D sobre S con peso 1/Var(D). Con la matriz de información [[Σw, Σws], [Σws, Σws²]], su inversa se
 * escribe con el determinante Δ = Σw·Σws² − (Σws)²: Var(α) = Σws²/Δ, Var(β) = Σw/Δ. `null` si Δ ≈ 0 (todos los S iguales).
 */
export function regresionMosesLittenberg(puntos: readonly PuntoDiagnostico[]): RegresionMosesLittenberg | null {
  if (puntos.length < 3) return null;
  let sw = 0, sws = 0, sws2 = 0, swd = 0, swsd = 0;
  for (const p of puntos) {
    const w = 1 / (p.seD * p.seD);
    sw += w; sws += w * p.s; sws2 += w * p.s * p.s; swd += w * p.dLn; swsd += w * p.s * p.dLn;
  }
  const delta = sw * sws2 - sws * sws;
  if (!(Math.abs(delta) > 1e-12 * Math.max(1, sw * sws2))) return null;
  const a = (sws2 * swd - sws * swsd) / delta;
  const b = (sw * swsd - sws * swd) / delta;
  const seA = Math.sqrt(Math.max(0, sws2 / delta));
  const seB = Math.sqrt(Math.max(0, sw / delta));
  return { a, b, seA, seB, zB: seB > 0 ? b / seB : 0 };
}

/** Curva SROC en `n` puntos de FPR entre 0,001 y 0,999 (escala logit lineal en L_F); `null` si |β| ≥ 1 (la curva no está definida). */
export function curvaMosesLittenberg(reg: { a: number; b: number }, n = 61): { fpr: number; tpr: number }[] | null {
  if (Math.abs(reg.b) >= 1) return null;
  const origen = reg.a / (1 - reg.b);
  const pendiente = (1 + reg.b) / (1 - reg.b);
  const curva: { fpr: number; tpr: number }[] = [];
  for (let i = 0; i < n; i++) {
    const fpr = 0.001 + (i / (n - 1)) * 0.998;
    const tpr = 1 / (1 + Math.exp(-(origen + pendiente * Math.log(fpr / (1 - fpr)))));
    if (Number.isFinite(tpr)) curva.push({ fpr, tpr });
  }
  return curva.length > 0 ? curva : null;
}
