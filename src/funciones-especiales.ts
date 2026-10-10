// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Funciones especiales de probabilidad del motor: ln Γ, gamma incompleta regularizada y beta incompleta regularizada.
 * De ellas salen el p de χ² (heterogeneidad, Q-profile, meta-regresión), el p y los cuantiles de la t de Student
 * (Hartung-Knapp, intervalos de predicción, Egger, Peters) y la cola de la normal.
 *
 * Sin dependencias (ni de otros módulos ni de librerías): se importa también desde scripts que se ejecutan con Node
 * directamente. Implementación propia a partir de las fórmulas matemáticas publicadas; la procedencia de cada una
 * está en `PROCEDENCIA.md` y la comprobación contra R en `scripts/test-funciones-especiales.mts`.
 *
 * Fuentes de las fórmulas (todas son resultados matemáticos, no código de ningún programa):
 *  - ln Γ: serie de Stirling con los términos de Bernoulli (Abramowitz y Stegun 6.1.40) y la recurrencia Γ(x+1) = xΓ(x);
 *    reflexión de Euler para x < ½ (A&S 6.1.17).
 *  - Gamma incompleta: serie de potencias para P (A&S 6.5.29) y fracción continua para Q (DLMF 8.9.2).
 *  - Beta incompleta: fracción continua de A&S 26.5.8 (DLMF 8.17.22), con la simetría I_x(a,b) = 1 − I_{1−x}(b,a).
 *  - Evaluación de fracciones continuas: método de Lentz modificado (Lentz 1976; Thompson y Barnett 1986).
 */

const MINIMO = 1e-300;
const LN_RAIZ_2PI = 0.9189385332046727; // ln √(2π)

/**
 * Valor de la fracción continua b0 + a1/(b1 + a2/(b2 + a3/(b3 + …))), con `a(n)` y `b(n)` para n ≥ 1,
 * por el método de Lentz modificado. Devuelve el valor al que converge (o el último, si no converge).
 */
function fraccionContinua(
  b0: number,
  a: (n: number) => number,
  b: (n: number) => number,
  tolerancia: number,
  maxIteraciones: number,
): number {
  let f = b0 === 0 ? MINIMO : b0;
  let C = f;
  let D = 0;
  for (let n = 1; n <= maxIteraciones; n++) {
    const an = a(n);
    const bn = b(n);
    D = bn + an * D;
    if (D === 0) D = MINIMO;
    C = bn + an / C;
    if (C === 0) C = MINIMO;
    D = 1 / D;
    const delta = C * D;
    f *= delta;
    if (Math.abs(delta - 1) < tolerancia) break;
  }
  return f;
}

/** ln Γ(x) para x > 0 (y para x no entero negativo, por reflexión). */
export function lnGamma(x: number): number {
  if (Number.isNaN(x)) return NaN;
  if (x <= 0 && Number.isInteger(x)) return Infinity;
  if (x < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * x))) - lnGamma(1 - x);
  // Se lleva x por encima de 15 con Γ(x) = Γ(x + n) / [x(x+1)…(x+n−1)], donde la serie de Stirling es muy exacta.
  let y = x;
  let lnProducto = 0;
  while (y < 15) {
    lnProducto += Math.log(y);
    y += 1;
  }
  const y2 = 1 / (y * y);
  // Σ B₂ₖ / (2k(2k−1) y^(2k−1)) = 1/(12y) − 1/(360y³) + 1/(1260y⁵) − 1/(1680y⁷) + 1/(1188y⁹) − 691/(360360y¹¹) + 1/(156y¹³)
  const serie =
    (1 / 12 - y2 * (1 / 360 - y2 * (1 / 1260 - y2 * (1 / 1680 - y2 * (1 / 1188 - y2 * (691 / 360360 - y2 / 156)))))) / y;
  return (y - 0.5) * Math.log(y) - y + LN_RAIZ_2PI + serie - lnProducto;
}

/** P(a, x): gamma incompleta inferior regularizada, por la serie de potencias (adecuada cuando x < a + 1). */
function serieGammaP(a: number, x: number): number {
  let termino = 1;
  let suma = 1;
  for (let n = 1; n < 10000; n++) {
    termino *= x / (a + n);
    suma += termino;
    if (termino < suma * 1e-17) break;
  }
  return suma * Math.exp(-x + a * Math.log(x) - lnGamma(a + 1));
}

/** Q(a, x): gamma incompleta superior regularizada, por la fracción continua de DLMF 8.9.2 (adecuada cuando x ≥ a + 1). */
function fraccionGammaQ(a: number, x: number): number {
  // 1/(x + (1−a)/(1 + 1/(x + (2−a)/(1 + 2/(x + …)))))
  const f = fraccionContinua(
    0,
    (n) => (n === 1 ? 1 : n % 2 === 0 ? n / 2 - a : (n - 1) / 2),
    (n) => (n % 2 === 0 ? 1 : x),
    1e-16,
    20000,
  );
  return f * Math.exp(-x + a * Math.log(x) - lnGamma(a));
}

/** Q(a, x) = Γ(a, x) / Γ(a): cola superior de la gamma. Para a > 0 y x ≥ 0. Exacta también en colas extremas. */
export function gammaIncompletaQ(a: number, x: number): number {
  if (Number.isNaN(a) || Number.isNaN(x) || a <= 0 || x < 0) return NaN;
  if (x === 0) return 1;
  if (x === Infinity) return 0;
  const q = x < a + 1 ? 1 - serieGammaP(a, x) : fraccionGammaQ(a, x);
  return Math.max(0, Math.min(1, q));
}

/** P(a, x) = γ(a, x) / Γ(a): cola inferior de la gamma. */
export function gammaIncompletaP(a: number, x: number): number {
  if (Number.isNaN(a) || Number.isNaN(x) || a <= 0 || x < 0) return NaN;
  if (x === 0) return 0;
  if (x === Infinity) return 1;
  const p = x < a + 1 ? serieGammaP(a, x) : 1 - fraccionGammaQ(a, x);
  return Math.max(0, Math.min(1, p));
}

/** Fracción continua de A&S 26.5.8: 1/(1 + d₁/(1 + d₂/(1 + …))), con d₂ₘ₊₁ y d₂ₘ de la fórmula de la beta incompleta. */
function fraccionBeta(x: number, a: number, b: number): number {
  return fraccionContinua(
    0,
    (n) => {
      if (n === 1) return 1;
      const k = n - 1; // índice de d_k
      const m = Math.floor(k / 2);
      return k % 2 === 1
        ? (-(a + m) * (a + b + m) * x) / ((a + 2 * m) * (a + 2 * m + 1))
        : (m * (b - m) * x) / ((a + 2 * m - 1) * (a + 2 * m));
    },
    () => 1,
    1e-16,
    20000,
  );
}

/**
 * I_x(a, b) con el complemento `xc = 1 − x` ya calculado por quien llama. Cuando x está muy cerca de 1, restar `1 − x`
 * pierde cifras; pasar `xc` directamente las conserva (lo necesita la t de Student para |t| muy pequeños).
 */
function betaIncompletaConComplemento(x: number, xc: number, a: number, b: number): number {
  if (x <= 0) return 0;
  if (xc <= 0) return 1;
  const prefijo = (xx: number, xxc: number, aa: number, bb: number) =>
    Math.exp(lnGamma(aa + bb) - lnGamma(aa) - lnGamma(bb) + aa * Math.log(xx) + bb * Math.log(xxc));
  // La fracción continua converge deprisa cuando x < (a + 1)/(a + b + 2) (A&S 26.5.8); si no, se usa la simetría.
  if (x < (a + 1) / (a + b + 2)) return Math.min(1, (prefijo(x, xc, a, b) * fraccionBeta(x, a, b)) / a);
  return Math.max(0, 1 - Math.min(1, (prefijo(xc, x, b, a) * fraccionBeta(xc, b, a)) / b));
}

/** I_x(a, b): beta incompleta regularizada, para 0 ≤ x ≤ 1, a > 0, b > 0. */
export function betaIncompleta(x: number, a: number, b: number): number {
  if (Number.isNaN(x) || Number.isNaN(a) || Number.isNaN(b) || a <= 0 || b <= 0 || x < 0 || x > 1) return NaN;
  return betaIncompletaConComplemento(x, 1 - x, a, b);
}

/** p de la cola superior de χ² con `gl` grados de libertad: P(χ² ≥ q). */
export function chiCuadradoCola(q: number, gl: number): number {
  if (gl <= 0 || Number.isNaN(gl) || Number.isNaN(q)) return NaN;
  if (q <= 0) return 1;
  return gammaIncompletaQ(gl / 2, q / 2);
}

/** P(|T| > t) para una t de Student con `gl` grados de libertad: la cola de AMBOS lados (ya es el p bilateral). */
export function tStudentDosColas(t: number, gl: number): number {
  if (gl <= 0 || Number.isNaN(gl) || Number.isNaN(t)) return NaN;
  if (t === 0) return 1;
  if (!Number.isFinite(t)) return 0;
  const t2 = t * t;
  return betaIncompletaConComplemento(gl / (gl + t2), t2 / (gl + t2), gl / 2, 0.5);
}

/** Cola de dos lados de la normal estándar: P(|Z| > z) = Q(½, z²/2). */
export function normalDosColas(z: number): number {
  if (Number.isNaN(z)) return NaN;
  return gammaIncompletaQ(0.5, (z * z) / 2);
}

/**
 * Valor t > 0 tal que P(|T| > t) = `alfa` para una t de Student con `gl` grados de libertad (cuantil bilateral).
 * Por bisección sobre `tStudentDosColas`, que baja de forma estricta al subir t. `alfa` en (0, 1).
 */
export function tStudentCuantilDosColas(alfa: number, gl: number): number {
  if (!(alfa > 0 && alfa < 1) || !(gl > 0)) return NaN;
  let lo = 0;
  let hi = 1;
  while (tStudentDosColas(hi, gl) > alfa && hi < 1e300) hi *= 2;
  for (let i = 0; i < 400; i++) {
    const mid = (lo + hi) / 2;
    if (tStudentDosColas(mid, gl) > alfa) lo = mid;
    else hi = mid;
    if (hi - lo <= 4e-16 * hi) break;
  }
  return (lo + hi) / 2;
}

/** Valor z > 0 tal que P(|Z| > z) = `alfa` (cuantil bilateral de la normal estándar). `alfa` en (0, 1). */
export function normalCuantilDosColas(alfa: number): number {
  if (!(alfa > 0 && alfa < 1)) return NaN;
  let lo = 0;
  let hi = 1;
  while (normalDosColas(hi) > alfa && hi < 1e3) hi *= 2;
  for (let i = 0; i < 400; i++) {
    const mid = (lo + hi) / 2;
    if (normalDosColas(mid) > alfa) lo = mid;
    else hi = mid;
    if (hi - lo <= 4e-16 * hi) break;
  }
  return (lo + hi) / 2;
}

/** Cuantil de la normal estándar Φ⁻¹(p), para p en (0, 1). */
export function normalCuantil(p: number): number {
  if (!(p > 0 && p < 1)) return NaN;
  if (p === 0.5) return 0;
  const z = normalCuantilDosColas(2 * Math.min(p, 1 - p));
  return p > 0.5 ? z : -z;
}
