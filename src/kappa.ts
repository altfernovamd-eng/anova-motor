// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
// ===========================================================================
//  κ de Cohen (2 evaluadores) y κ de Fleiss (≥3 evaluadores), sin ponderar — matemática pura, sin dependencias
//  de la app (para poder validarla con Node contra R); solo usa `funciones-especiales.ts`, también autocontenido.
//
//  Fórmulas (todas publicadas; la procedencia está en docs/PROCEDENCIA.md):
//   · Cohen 1960: κ = (po − pe) / (1 − pe), con po = acuerdo observado y pe = acuerdo esperado por azar.
//   · Fleiss 1971: la misma expresión, con po = promedio por sujeto de la proporción de pares de evaluadores que coinciden
//     y pe = Σⱼ pⱼ² (pⱼ = proporción de todas las calificaciones que cayeron en la categoría j).
//   · Error estándar para contrastar H₀: κ = 0 (la «z» que se muestra), desarrollando la varianza bajo H₀ de
//     Fleiss, Cohen y Everitt (1969) para dos evaluadores y de Fleiss, Nee y Landis (1979) para varios:
//        Cohen:  Var₀ = [ pe + pe² − Σᵢ pᵢ. p.ᵢ (pᵢ. + p.ᵢ) ] / [ N (1 − pe)² ]
//        Fleiss: Var₀ = 2 [ pe + pe² − 2 Σⱼ pⱼ³ ] / [ N n (n − 1) (1 − pe)² ]      (n = evaluadores)
//     OJO: ese error estándar vale SOLO para contrastar «no hay acuerdo». No es el de un intervalo de confianza de κ.
//   · Interpretación cualitativa: Landis y Koch 1977.
//  Comprobado contra el paquete `irr` de R (kappa2 y kappam.fleiss): `npm run test:kappa-r` y `npm run benchmark:kappa`;
//  ver VALIDACION_ESTADISTICA.md §3.3.
// ===========================================================================

import { normalDosColas } from "./funciones-especiales.ts";

export type MetodoKappa = "Cohen" | "Fleiss";

export interface ResultadoKappa {
  metodo: MetodoKappa;
  kappa: number;
  /** Acuerdo observado (proporción). Para 2 revisores, el % de acuerdo bruto. */
  acuerdoObservado: number;
  /** Acuerdo esperado por azar. */
  acuerdoEsperado: number;
  /** Error estándar bajo H₀: κ = 0. Sirve para el contraste (z, p), no para un intervalo de confianza. */
  se: number;
  z: number;
  pValor: number;
  nSujetos: number;
  nEvaluadores: number;
  /** Interpretación cualitativa (Landis & Koch 1977). */
  interpretacion: string;
}

/** Interpretación cualitativa de κ según Landis & Koch (1977). */
export function interpretarKappa(k: number): string {
  if (k < 0) return "Pobre (peor que el azar)";
  if (k < 0.2) return "Leve";
  if (k < 0.4) return "Aceptable";
  if (k < 0.6) return "Moderada";
  if (k < 0.8) return "Considerable";
  return "Casi perfecta";
}

/** Etiquetas distintas presentes, en un orden estable. */
function categorias(valores: string[]): string[] {
  return [...new Set(valores)].sort();
}

/** Arma el resultado a partir de κ, el acuerdo y la varianza bajo H₀ (z y p solo si hay varianza). */
function resultado(
  metodo: MetodoKappa,
  po: number,
  pe: number,
  varianzaNula: number,
  nSujetos: number,
  nEvaluadores: number,
): ResultadoKappa {
  const kappa = (po - pe) / (1 - pe);
  const se = Math.sqrt(Math.max(0, varianzaNula));
  const z = se > 0 ? kappa / se : 0;
  return {
    metodo,
    kappa,
    acuerdoObservado: po,
    acuerdoEsperado: pe,
    se,
    z,
    pValor: se > 0 ? normalDosColas(z) : 1,
    nSujetos,
    nEvaluadores,
    interpretacion: interpretarKappa(kappa),
  };
}

/**
 * κ de Cohen para 2 evaluadores (sin ponderar).
 * `pares`: una fila por sujeto con las dos calificaciones [evaluadorA, evaluadorB].
 */
export function cohenKappa(pares: [string, string][]): ResultadoKappa | null {
  const N = pares.length;
  if (N < 1) return null;
  const etiquetas = categorias(pares.flat());
  const posicion = new Map(etiquetas.map((e, i) => [e, i]));
  const k = etiquetas.length;

  // Recuentos: cuántas veces A dijo la categoría i, B la j (coincidencias) y las marginales de cada evaluador.
  const marginalA = new Array<number>(k).fill(0);
  const marginalB = new Array<number>(k).fill(0);
  let coincidencias = 0;
  for (const [a, b] of pares) {
    const i = posicion.get(a)!;
    const j = posicion.get(b)!;
    marginalA[i] += 1;
    marginalB[j] += 1;
    if (i === j) coincidencias += 1;
  }

  const po = coincidencias / N;
  let pe = 0;
  let sumaCruzada = 0; // Σᵢ pᵢ. p.ᵢ (pᵢ. + p.ᵢ)
  for (let i = 0; i < k; i++) {
    const pA = marginalA[i] / N;
    const pB = marginalB[i] / N;
    pe += pA * pB;
    sumaCruzada += pA * pB * (pA + pB);
  }
  if (1 - pe === 0) return null; // todos los sujetos en una sola categoría: κ no está definido

  const varianzaNula = (pe + pe * pe - sumaCruzada) / (N * (1 - pe) ** 2);
  return resultado("Cohen", po, pe, varianzaNula, N, 2);
}

/**
 * κ de Fleiss para n evaluadores (sin ponderar).
 * `sujetos`: una fila por sujeto con las calificaciones de todos los evaluadores (el mismo nº por sujeto).
 */
export function fleissKappa(sujetos: string[][]): ResultadoKappa | null {
  const N = sujetos.length;
  if (N < 1) return null;
  const n = sujetos[0].length;
  if (n < 2 || !sujetos.every((s) => s.length === n)) return null;
  const etiquetas = categorias(sujetos.flat());
  const posicion = new Map(etiquetas.map((e, i) => [e, i]));
  const k = etiquetas.length;

  // Para cada sujeto: cuántos evaluadores eligieron cada categoría; de ahí la proporción de pares que coinciden.
  const total = new Array<number>(k).fill(0);
  let sumaAcuerdoPorSujeto = 0;
  for (const calificaciones of sujetos) {
    const elegidas = new Array<number>(k).fill(0);
    for (const c of calificaciones) elegidas[posicion.get(c)!] += 1;
    let paresQueCoinciden = 0;
    for (let j = 0; j < k; j++) {
      total[j] += elegidas[j];
      paresQueCoinciden += elegidas[j] * (elegidas[j] - 1);
    }
    sumaAcuerdoPorSujeto += paresQueCoinciden / (n * (n - 1));
  }

  const po = sumaAcuerdoPorSujeto / N;
  let pe = 0;
  let sumaCubos = 0; // Σⱼ pⱼ³
  for (let j = 0; j < k; j++) {
    const p = total[j] / (N * n);
    pe += p * p;
    sumaCubos += p * p * p;
  }
  if (1 - pe === 0) return null;

  const varianzaNula = (2 * (pe + pe * pe - 2 * sumaCubos)) / (N * n * (n - 1) * (1 - pe) ** 2);
  return resultado("Fleiss", po, pe, varianzaNula, N, n);
}
