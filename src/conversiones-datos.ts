// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Conversiones de datos del Cochrane Handbook, capítulo 6 y sección 10.6: obtener la DE (o el error estándar) cuando el
 * estudio reporta otra cosa, combinar dos grupos, imputar la DE del cambio desde la línea base, estimar media y DE desde la
 * mediana y los cuartiles o el rango, y convertir entre OR y diferencia de medias estandarizada.
 *
 * Son aproximaciones o identidades algebraicas, nunca datos medidos: cada función dice de qué sección sale y cuándo no
 * conviene usarla. Devuelven `null` ante datos imposibles en vez de inventar un número. Lógica pura (sin React).
 *
 * Fuentes: Handbook 6.3.1–6.3.2 (EE desde IC y P), 6.5.2.2 (DE desde EE o IC de una media), 6.5.2.3 (DE desde la diferencia
 * de medias), 6.5.2.5–6.5.2.6 (IQR y rango), 6.5.2.8 (DE del cambio), 6.5.2.9–6.5.2.10 (medias faltantes, combinar grupos),
 * 10.6 (OR ↔ DME); Wan, Wang, Liu y Tong (2014), BMC Med Res Methodol 14:135.
 */
import { normalCuantil, normalCuantilDosColas, tStudentCuantilDosColas } from "./funciones-especiales.ts";

const esNumero = (...x: number[]) => x.every((v) => Number.isFinite(v));

/** Distribución con la que el estudio calculó su IC o su p: la t (lo correcto con muestras pequeñas) o la normal. */
export type DistribucionIC = "t" | "normal";
export interface OpcionesIC {
  /** Nivel del IC reportado (por defecto 0,95). */
  nivel?: number;
  /** Por defecto la t, que es exacta si el estudio la usó y casi igual a la normal con muestras grandes. */
  distribucion?: DistribucionIC;
}

/** Multiplicador del IC: z o t bilateral para el nivel dado. */
function multiplicador(nivel: number, gl: number, distribucion: DistribucionIC): number {
  const alfa = 1 - nivel;
  return distribucion === "normal" ? normalCuantilDosColas(alfa) : tStudentCuantilDosColas(alfa, gl);
}
const nivelValido = (n: number) => n > 0 && n < 1;

/** DE de un grupo desde el error estándar de su media: SD = SE · √N (6.5.2.2). */
export function sdDesdeEE(ee: number, n: number): number | null {
  if (!esNumero(ee, n) || ee < 0 || n < 1) return null;
  return ee * Math.sqrt(n);
}

/**
 * DE de un grupo desde el IC de su media (6.5.2.2): SE = (límite superior − inferior) / (2·q), SD = SE·√N, con q el cuantil de
 * la t con N − 1 grados de libertad (o de la normal; 1,96 al 95 %: «3,92 de ancho»). El Handbook avisa de que con menos de
 * unos 60 participantes por grupo el IC se calcula con la t.
 */
export function sdDesdeICMedia(inf: number, sup: number, n: number, opc: OpcionesIC = {}): { se: number; sd: number } | null {
  const nivel = opc.nivel ?? 0.95;
  if (!esNumero(inf, sup, n) || sup <= inf || n < 2 || !nivelValido(nivel)) return null;
  const q = multiplicador(nivel, n - 1, opc.distribucion ?? "t");
  const se = (sup - inf) / (2 * q);
  return { se, sd: se * Math.sqrt(n) };
}

/** DE intragrupo desde el error estándar de la diferencia de medias: SD = SE / √(1/N₁ + 1/N₂) (6.5.2.3). */
export function sdDesdeEEDiferencia(se: number, nE: number, nC: number): number | null {
  if (!esNumero(se, nE, nC) || se < 0 || nE < 1 || nC < 1) return null;
  return se / Math.sqrt(1 / nE + 1 / nC);
}

/** Error estándar de la diferencia de medias desde su IC (6.5.2.3), con t de N₁ + N₂ − 2 grados de libertad (o normal). */
export function eeDiferenciaDesdeIC(inf: number, sup: number, nE: number, nC: number, opc: OpcionesIC = {}): number | null {
  const nivel = opc.nivel ?? 0.95;
  if (!esNumero(inf, sup, nE, nC) || sup <= inf || nE < 1 || nC < 1 || nE + nC < 3 || !nivelValido(nivel)) return null;
  return (sup - inf) / (2 * multiplicador(nivel, nE + nC - 2, opc.distribucion ?? "t"));
}

/** Error estándar de la diferencia de medias desde el estadístico t: SE = |DM| / |t| (6.5.2.3). */
export function eeDiferenciaDesdeT(diferencia: number, t: number): number | null {
  if (!esNumero(diferencia, t) || t === 0) return null;
  return Math.abs(diferencia) / Math.abs(t);
}

/**
 * Error estándar de la diferencia de medias desde el valor de P de dos colas (6.5.2.3): primero la t que corresponde a ese P
 * con N₁ + N₂ − 2 grados de libertad, luego SE = |DM| / t. Un P aproximado («p < 0,05», «p = 0,01») da un SE aproximado.
 * `null` si P no está en (0, 1) o la diferencia es 0 (con DM = 0 el SE no se puede recuperar).
 */
export function eeDiferenciaDesdeP(diferencia: number, p: number, nE: number, nC: number): number | null {
  if (!esNumero(diferencia, p, nE, nC) || !(p > 0 && p < 1) || diferencia === 0 || nE < 1 || nC < 1 || nE + nC < 3) return null;
  return Math.abs(diferencia) / tStudentCuantilDosColas(p, nE + nC - 2);
}

/**
 * Error estándar desde el IC de una medida ABSOLUTA con distribución normal (DM, DME, diferencia de riesgos), 6.3.1:
 * SE = (sup − inf) / 3,92 al 95 %, 3,29 al 90 %, 5,15 al 99 % (es decir, 2·z).
 */
export function eeDesdeICAbsoluto(inf: number, sup: number, nivel = 0.95): number | null {
  if (!esNumero(inf, sup) || sup <= inf || !nivelValido(nivel)) return null;
  return (sup - inf) / (2 * normalCuantilDosColas(1 - nivel));
}

/** Error estándar del logaritmo de una razón (RR, OR, HR) desde su IC, 6.3.2: (ln sup − ln inf) / 3,92 al 95 %. */
export function eeLogDesdeICRazon(inf: number, sup: number, nivel = 0.95): number | null {
  if (!esNumero(inf, sup) || inf <= 0 || sup <= inf || !nivelValido(nivel)) return null;
  return (Math.log(sup) - Math.log(inf)) / (2 * normalCuantilDosColas(1 - nivel));
}

/** Error estándar desde el valor de P de dos colas de una medida ABSOLUTA (6.3.1): z desde P y SE = |estimación| / z. */
export function eeDesdePAbsoluto(estimacion: number, p: number): number | null {
  if (!esNumero(estimacion, p) || !(p > 0 && p < 1) || estimacion === 0) return null;
  return Math.abs(estimacion) / normalCuantilDosColas(p);
}

/** IC (inferior, superior) de una medida ABSOLUTA desde su estimación y su error estándar, con la normal al nivel dado. */
export function icDesdeEEAbsoluto(estimacion: number, ee: number, nivel = 0.95): { inf: number; sup: number } | null {
  if (!esNumero(estimacion, ee) || ee <= 0 || !nivelValido(nivel)) return null;
  const z = normalCuantilDosColas(1 - nivel);
  return { inf: estimacion - z * ee, sup: estimacion + z * ee };
}

/** IC (inferior, superior) de una RAZÓN desde la razón y el error estándar de su logaritmo: exp(ln(razón) ± z · EE). */
export function icDesdeEELogRazon(razon: number, eeLog: number, nivel = 0.95): { inf: number; sup: number } | null {
  if (!esNumero(razon, eeLog) || razon <= 0 || eeLog <= 0 || !nivelValido(nivel)) return null;
  const z = normalCuantilDosColas(1 - nivel);
  return { inf: Math.exp(Math.log(razon) - z * eeLog), sup: Math.exp(Math.log(razon) + z * eeLog) };
}

/** Error estándar del log de una razón desde su P de dos colas (6.3.2): SE = |ln(razón)| / z. */
export function eeLogDesdePRazon(razon: number, p: number): number | null {
  if (!esNumero(razon, p) || razon <= 0 || razon === 1 || !(p > 0 && p < 1)) return null;
  return Math.abs(Math.log(razon)) / normalCuantilDosColas(p);
}

export interface GrupoResumen {
  n: number;
  media: number;
  sd: number;
}

/**
 * Combina dos grupos en uno solo (Handbook Tabla 6.5.a), como si nunca se hubieran dividido: N = N₁ + N₂, media ponderada y
 *   SD² = [ (N₁−1)·SD₁² + (N₂−1)·SD₂² + N₁N₂/(N₁+N₂) · (M₁² + M₂² − 2·M₁·M₂) ] / (N₁ + N₂ − 1).
 * Con más de dos grupos se aplica en secuencia. Sirve, por ejemplo, para fundir dos brazos experimentales de un mismo ensayo.
 */
export function combinarGrupos(a: GrupoResumen, b: GrupoResumen): GrupoResumen | null {
  if (!esNumero(a.n, a.media, a.sd, b.n, b.media, b.sd) || a.n < 1 || b.n < 1 || a.sd < 0 || b.sd < 0 || a.n + b.n < 2) return null;
  const n = a.n + b.n;
  const media = (a.n * a.media + b.n * b.media) / n;
  const suma = (a.n - 1) * a.sd ** 2 + (b.n - 1) * b.sd ** 2 + ((a.n * b.n) / n) * (a.media ** 2 + b.media ** 2 - 2 * a.media * b.media);
  return { n, media, sd: Math.sqrt(Math.max(0, suma / (n - 1))) };
}

/** DE del cambio desde la línea base imputada con una correlación (6.5.2.8): √(SD_base² + SD_final² − 2·Corr·SD_base·SD_final). */
export function sdCambioDesdeCorrelacion(sdBase: number, sdFinal: number, correlacion: number): number | null {
  if (!esNumero(sdBase, sdFinal, correlacion) || sdBase < 0 || sdFinal < 0 || correlacion < -1 || correlacion > 1) return null;
  return Math.sqrt(Math.max(0, sdBase ** 2 + sdFinal ** 2 - 2 * correlacion * sdBase * sdFinal));
}

/** Correlación base–final que implican tres DE de un mismo grupo (6.5.2.8): (SD_base² + SD_final² − SD_cambio²) / (2·SD_base·SD_final). */
export function correlacionDesdeSdCambio(sdBase: number, sdFinal: number, sdCambio: number): number | null {
  if (!esNumero(sdBase, sdFinal, sdCambio) || sdBase <= 0 || sdFinal <= 0 || sdCambio < 0) return null;
  const r = (sdBase ** 2 + sdFinal ** 2 - sdCambio ** 2) / (2 * sdBase * sdFinal);
  // Un valor fuera de [−1, 1] no es una correlación: los tres datos no son coherentes entre sí.
  return r >= -1 && r <= 1 ? r : null;
}

/** DE aproximada desde el rango intercuartílico cuando la distribución es normal y la muestra grande (6.5.2.5): IQR / 1,35. */
export function sdDesdeIQRNormal(q1: number, q3: number): number | null {
  if (!esNumero(q1, q3) || q3 <= q1) return null;
  return (q3 - q1) / 1.35;
}

export type EscenarioWan = "min-mediana-max" | "cuartiles" | "cinco-numeros";
export interface EntradaWan {
  n: number;
  mediana: number;
  q1?: number;
  q3?: number;
  min?: number;
  max?: number;
}
export interface ResultadoWan {
  escenario: EscenarioWan;
  media: number;
  sd: number;
  advertencias: string[];
}

/** 2·Φ⁻¹((n − 0,375)/(n + 0,25)): ξ(n) del artículo de Wan (ecuación 9), por la aproximación de Blom. */
const xi = (n: number) => 2 * normalCuantil((n - 0.375) / (n + 0.25));
/** 2·Φ⁻¹((0,75n − 0,125)/(n + 0,25)): η(n) del artículo de Wan (ecuación 16), por la aproximación de Blom. */
const eta = (n: number) => 2 * normalCuantil((0.75 * n - 0.125) / (n + 0.25));

/**
 * Media y DE estimadas desde la mediana y el rango y/o los cuartiles (Wan et al. 2014, que el Handbook 6.5.2.5 y 6.5.2.9
 * recomienda). Tres escenarios según lo que reporte el estudio:
 *  · mínimo, mediana, máximo:  media ≈ (a + 2m + b)/4 + (a − 2m + b)/(4n) (Hozo et al., ec. 2 del artículo);  DE ≈ (b − a) / ξ(n).
 *  · cuartiles:                media ≈ (q₁ + m + q₃)/3;  DE ≈ (q₃ − q₁) / η(n).
 *  · cinco números:            media ≈ (a + 2q₁ + 2m + 2q₃ + b)/8;  DE ≈ ½ [ (b − a)/ξ(n) + (q₃ − q₁)/η(n) ].
 * Suponen datos normales; con datos asimétricos son una aproximación (ligeramente sesgada en el artículo) y con tamaños muy
 * pequeños son poco fiables. Si se dan más datos de los necesarios se usa el escenario de los cinco números solo si están
 * todos; ante cualquier otra combinación devuelve `null`. Usa las fórmulas finales del artículo (aproximación de Blom) para
 * cualquier n; `meta::metamean` usa además los valores tabulados del artículo para n ≤ 50, de modo que con n pequeños los
 * resultados pueden diferir unas décimas de por ciento (ver VALIDACION_ESTADISTICA.md).
 */
export function mediaSdWan(e: EntradaWan): ResultadoWan | null {
  const { n, mediana: m, q1, q3, min: a, max: b } = e;
  if (!esNumero(n, m) || n < 5 || !Number.isInteger(n)) return null;
  const tieneRango = a !== undefined && b !== undefined;
  const tieneCuartiles = q1 !== undefined && q3 !== undefined;
  const advertencias: string[] = [];
  if (n < 15) advertencias.push("Con menos de 15 participantes por grupo estas aproximaciones son poco fiables.");

  if (tieneRango && tieneCuartiles) {
    if (![a!, q1!, m, q3!, b!].every(Number.isFinite) || !(a! <= q1! && q1! <= m && m <= q3! && q3! <= b!) || b! <= a!) return null;
    return {
      escenario: "cinco-numeros",
      media: (a! + 2 * q1! + 2 * m + 2 * q3! + b!) / 8,
      sd: 0.5 * ((b! - a!) / xi(n) + (q3! - q1!) / eta(n)),
      advertencias,
    };
  }
  if (tieneRango && q1 === undefined && q3 === undefined) {
    if (![a!, m, b!].every(Number.isFinite) || !(a! <= m && m <= b!) || b! <= a!) return null;
    advertencias.push("El rango es inestable: crece con el tamaño de la muestra y un solo valor extremo lo cambia (Handbook 6.5.2.6). Úsalo solo si no hay cuartiles.");
    return { escenario: "min-mediana-max", media: (a! + 2 * m + b!) / 4 + (a! - 2 * m + b!) / (4 * n), sd: (b! - a!) / xi(n), advertencias };
  }
  if (tieneCuartiles && a === undefined && b === undefined) {
    if (![q1!, m, q3!].every(Number.isFinite) || !(q1! <= m && m <= q3!) || q3! <= q1!) return null;
    return { escenario: "cuartiles", media: (q1! + m + q3!) / 3, sd: (q3! - q1!) / eta(n), advertencias };
  }
  return null;
}

const RAIZ3_ENTRE_PI = Math.sqrt(3) / Math.PI;

/** Odds ratio → diferencia de medias estandarizada (10.6): DME = ln(OR) · √3/π; el EE se multiplica por la misma constante. */
export function dmeDesdeOR(or: number, eeLogOR: number): { dme: number; ee: number } | null {
  if (!esNumero(or, eeLogOR) || or <= 0 || eeLogOR <= 0) return null;
  return { dme: Math.log(or) * RAIZ3_ENTRE_PI, ee: eeLogOR * RAIZ3_ENTRE_PI };
}

/** Diferencia de medias estandarizada → odds ratio (10.6): ln(OR) = DME · π/√3; el EE del log OR se multiplica por la misma constante. */
export function orDesdeDME(dme: number, ee: number): { or: number; eeLogOR: number } | null {
  if (!esNumero(dme, ee) || ee <= 0) return null;
  return { or: Math.exp(dme / RAIZ3_ENTRE_PI), eeLogOR: ee / RAIZ3_ENTRE_PI };
}
