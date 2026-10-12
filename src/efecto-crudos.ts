// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Piezas puras de la extracción (Paso 5) que también necesita el manuscrito (Paso 9): cálculo del efecto desde datos
 * crudos, cero eventos y comparación contra la verificación de un segundo extractor. Viven aquí, y no dentro de la
 * pantalla, para que lo que se calcula en el Paso 5 y lo que se declara en Métodos sea literalmente lo mismo.
 */
import {
  divisorControl,
  efectoDesde2x2,
  efectoDesdeCorrelacion,
  efectoDesdeOEV,
  efectoDesdeRazonDeMedias,
  efectoDesdeTasas,
  efectoDesdeCruzado,
  efectoDesdeMedias,
  efectoDesdeMediasEstandarizada,
  efectoDesdeProporcion,
  type EfectoCalculado,
} from "./statistics";
import type { DatosCrudos, Estudio, Medida } from "./tipos.ts";

export const num = (s: string): number => {
  const n = parseFloat(s);
  return Number.isFinite(n) ? n : NaN;
};

/** Cochrane Handbook cap. 19.5.2: un "0" reportado puede ser un evento realmente ausente (con
 *  vigilancia activa) o solo ausencia de reporte — solo aplica a los tipos con conteo de eventos. */
export function tieneEventoCero(crudos: DatosCrudos): boolean {
  if (crudos.tipo === "2x2") return num(crudos.eventosTrat) === 0 || num(crudos.eventosCtrl) === 0;
  if (crudos.tipo === "proporcion") return num(crudos.eventos) === 0;
  return false;
}
/** ¿Hay un cero reportado (tabla 2×2 o proporción) sin justificación escrita? Handbook 19.5.2: un 0 puede ser ausencia real del evento o ausencia de reporte, y debe justificarse. */
export function faltaJustificarCero(crudos: DatosCrudos): boolean {
  if (!tieneEventoCero(crudos)) return false;
  const justificacion = crudos.tipo === "2x2" || crudos.tipo === "proporcion" ? crudos.justificacionCero : undefined;
  return !justificacion?.trim();
}

/**
 * Calcula la medida de efecto (OR, RR, RAR, PR, DM, DME, Proporción) a partir de datos crudos.
 * @param medida - La medida a calcular (OR/RR/RAR/PR para 2×2; DM/DME para medias; Proporción para un solo brazo)
 * @param crudos - Datos crudos (2x2: eventos/totales; medias: media/DE/n; proporcion: eventos/total)
 * @param dividirControlEn - Si el grupo control se comparte con otras comparaciones del mismo ensayo (Handbook 23.3.4):
 *   en cuántas partes iguales se divide. En una tabla 2×2 se dividen los eventos Y el total del control; en medias solo
 *   el N (la media y la DE no cambian). Sin valor (o 1) no se divide.
 * @returns Efecto calculado (estimación, IC, SE) o null si no se puede calcular
 */
export function calcularEfectoCrudos(medida: Medida, crudos?: DatosCrudos, dividirControlEn?: number): EfectoCalculado | null {
  if (!crudos) return null;
  const k = divisorControl(dividirControlEn);
  if (crudos.tipo === "2x2" && (medida === "OR" || medida === "RR" || medida === "RAR" || medida === "PR")) {
    return efectoDesde2x2(
      num(crudos.eventosTrat),
      num(crudos.totalTrat),
      num(crudos.eventosCtrl) / k,
      num(crudos.totalCtrl) / k,
      medida,
    );
  }
  if (crudos.tipo === "medias" && (medida === "DM" || medida === "DME")) {
    const fn = medida === "DME" ? efectoDesdeMediasEstandarizada : efectoDesdeMedias;
    return fn(
      num(crudos.mediaTrat),
      num(crudos.deTrat),
      num(crudos.nTrat),
      num(crudos.mediaCtrl),
      num(crudos.deCtrl),
      num(crudos.nCtrl) / k,
    );
  }
  if (crudos.tipo === "medias" && medida === "RM") {
    return efectoDesdeRazonDeMedias(
      num(crudos.mediaTrat),
      num(crudos.deTrat),
      num(crudos.nTrat),
      num(crudos.mediaCtrl),
      num(crudos.deCtrl),
      num(crudos.nCtrl) / k,
    );
  }
  if (crudos.tipo === "tasas" && (medida === "RTI" || medida === "DTI")) {
    return efectoDesdeTasas(num(crudos.eventosTrat), num(crudos.tiempoTrat), num(crudos.eventosCtrl), num(crudos.tiempoCtrl), medida);
  }
  if (crudos.tipo === "correlacion" && medida === "r") return efectoDesdeCorrelacion(num(crudos.r), num(crudos.n));
  if (crudos.tipo === "oev" && medida === "HR") return efectoDesdeOEV(num(crudos.oe), num(crudos.v));
  if (crudos.tipo === "proporcion" && medida === "Proporción") {
    return efectoDesdeProporcion(num(crudos.eventos), num(crudos.total));
  }
  if (crudos.tipo === "cruzado" && medida === "DM") {
    return efectoDesdeCruzado(
      num(crudos.mediaTrat),
      num(crudos.mediaCtrl),
      num(crudos.n),
      num(crudos.sdDiff),
      num(crudos.sdTrat),
      num(crudos.sdCtrl),
      num(crudos.correlacion),
    );
  }
  return null;
}

export const EPSILON_VERIFICACION = 1e-6;

/**
 * Compara la extracción principal contra la verificación de un segundo revisor (si existe).
 * @param e - Estudio con `.verificacion` opcional
 * @returns null si no hay verificación cargada; si no, {completa, coincide, discrepancias}
 */
export function compararVerificacion(
  e: Estudio,
): { completa: boolean; coincide: boolean; discrepancias: string[] } | null {
  if (!e.verificacion) return null;
  const campos: { clave: "valor" | "icInferior" | "icSuperior"; etiqueta: string }[] = [
    { clave: "valor", etiqueta: "Estimación" },
    { clave: "icInferior", etiqueta: "IC inferior" },
    { clave: "icSuperior", etiqueta: "IC superior" },
  ];
  const completa = campos.every((c) => e.verificacion![c.clave].trim() !== "");
  if (!completa) return { completa: false, coincide: false, discrepancias: [] };
  const discrepancias: string[] = [];
  for (const c of campos) {
    const a = parseFloat(e[c.clave] as string);
    const b = parseFloat(e.verificacion![c.clave]);
    if (!Number.isFinite(a) || !Number.isFinite(b) || Math.abs(a - b) > EPSILON_VERIFICACION) {
      discrepancias.push(c.etiqueta);
    }
  }
  return { completa: true, coincide: discrepancias.length === 0, discrepancias };
}

