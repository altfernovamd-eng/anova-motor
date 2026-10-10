// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Tipos de datos que consume el motor estadístico.
 * Subconjunto de `src/types/meta-analysis.ts` de la aplicación aNova, generado por `scripts/exportar-motor.mts`
 * (15 declaraciones: las que importa el motor y las que éstas necesitan). No editar a mano.
 */

/** Campo categórico del estudio por el que se puede subagrupar en el Paso 7. */
export type CampoSubgrupo = "diseño" | "calidad" | "evidencia";

export type DisenoEstudio =
  | "ECA"
  | "Cohortes"
  | "Casos y Controles"
  | "Transversal"
  | "Prevalencia";

export type Medida =
  | "RR"
  | "RRR"
  | "RAR"
  | "NNT"
  | "DM"
  | "DME"
  | "HR"
  | "OR"
  | "PR"
  | "Proporción"
  /** Razón de tasas de incidencia (eventos y tiempo-persona; Cochrane Handbook 6.7, 10.8). */
  | "RTI"
  /** Diferencia de tasas de incidencia. */
  | "DTI"
  /** Razón de medias (Handbook 6.5.1.3). */
  | "RM"
  /** Coeficiente de correlación r, combinado con la transformación z de Fisher. */
  | "r";

export type RiesgoSesgo = "Bajo Riesgo" | "Poco Claro" | "Alto Riesgo";

/** Juicio de un dominio individual: el `valor` del nivel definido por cada herramienta (p. ej. "Bajo", "Serio"). */
export type NivelDominio = string;

export type CertezaGrade = "Alta" | "Moderada" | "Baja" | "Muy Baja";

/** Datos crudos capturados del estudio, guardados para poder recalcular y para el test de Peters. */
export type DatosCrudos =
  | {
      tipo: "2x2";
      eventosTrat: string;
      totalTrat: string;
      eventosCtrl: string;
      totalCtrl: string;
      /**
       * Cochrane Handbook cap. 19.5.2: "We generally recommend against extracting data as 'zero'
       * unless it is clearly listed as such in the study report" — un "0" puede ser un evento
       * realmente ausente (con vigilancia activa) o solo ausencia de reporte. Obligatorio solo
       * cuando `eventosTrat` o `eventosCtrl` es 0; registra cómo se determinó, per el Handbook
       * ("review authors should record how they reached a decision of 'zero events'").
       */
      justificacionCero?: string;
    }
  | {
      tipo: "medias";
      mediaTrat: string;
      deTrat: string;
      nTrat: string;
      mediaCtrl: string;
      deCtrl: string;
      nCtrl: string;
    }
  | {
      /** Razón o diferencia de tasas de incidencia: eventos y tiempo-persona de cada brazo (Handbook 6.7). */
      tipo: "tasas";
      eventosTrat: string;
      tiempoTrat: string;
      eventosCtrl: string;
      tiempoCtrl: string;
    }
  | {
      /** Coeficiente de correlación r de un estudio y su tamaño de muestra. */
      tipo: "correlacion";
      r: string;
      n: string;
    }
  | {
      /** Tiempo hasta el evento: O−E (observados menos esperados) y su varianza V, como los da una prueba de rangos logarítmicos (Handbook 6.8, 10.9). */
      tipo: "oev";
      oe: string;
      v: string;
    }
  | {
      /** Un solo brazo (p.ej. prevalencia): eventos sobre el total, sin grupo de comparación. */
      tipo: "proporcion";
      eventos: string;
      total: string;
      /** Ver nota de `justificacionCero` en la variante "2x2" — mismo requisito, aplicado a `eventos`. */
      justificacionCero?: string;
    }
  | {
      /**
       * Ensayo cruzado (crossover), desenlace continuo (solo DM — Cochrane Handbook cap. 23.2.7.1):
       * el MISMO grupo de `n` participantes bajo ambas condiciones, no dos grupos independientes.
       * `sdDiff` si el estudio la reporta directamente; si no, `sdTrat`/`sdCtrl`/`correlacion` la
       * imputan. Ver `efectoDesdeCruzado` en `statistics.ts`.
       */
      tipo: "cruzado";
      n: string;
      mediaTrat: string;
      mediaCtrl: string;
      sdDiff: string;
      sdTrat: string;
      sdCtrl: string;
      correlacion: string;
      /** De dónde sale la correlación cuando se imputa (Handbook 23.2.7.3): la reporta el estudio, se calculó con otro estudio
       *  de la revisión que presenta las tres DE, viene de una fuente externa, o es un supuesto razonado. Ausente = sin indicar. */
      origenCorrelacion?: OrigenCorrelacion;
    };

export type OrigenCorrelacion = "estudio" | "otro_estudio" | "externo" | "supuesto";

/**
 * Variante de ensayo aleatorizado (Cochrane Handbook cap. 23) que requiere un ajuste estadístico
 * distinto del análisis estándar de dos grupos independientes — de lo contrario se comete un "error
 * de unidad de análisis" (cap. 23.1.1): el IC queda artificialmente estrecho y el estudio pesa más de
 * lo que debería en el combinado. Mutuamente excluyente con el análisis estándar (un ensayo
 * clúster-cruzado híbrido, poco común, no se modela aparte — simplificación deliberada).
 */
export type VarianteEnsayo =
  | {
      tipo: "cluster";
      /** Número total de clústeres (ambos brazos). M (tamaño medio del clúster) = `Estudio.n` / esto. */
      numClusters: string;
      /** Coeficiente de correlación intraclase (ICC). Si el estudio no lo reporta, el Handbook
       *  recomienda una estimación externa de la literatura + análisis de sensibilidad. */
      icc: string;
      /** De dónde sale el ICC: el estudio lo reporta, se tomó de una fuente externa o de estudios parecidos, o es un
       *  supuesto razonado. Ausente = sin indicar. Los dos últimos son los que el Handbook (23.1.6) pide contrastar con un
       *  análisis de sensibilidad. */
      origenIcc?: OrigenIcc;
    }
  | { tipo: "cruzado" };

export type OrigenIcc = "estudio" | "externo" | "supuesto";

export interface Estudio {
  id: number;
  autor: string;
  año: string;
  diseño: DisenoEstudio;
  n: string;
  medida: Medida;
  valor: string;
  icInferior: string;
  icSuperior: string;
  /**
   * Variante de ensayo aleatorizado (clúster / cruzado, Cochrane Handbook cap. 23), solo relevante si
   * `diseño === "ECA"`. Ausente = análisis estándar de dos grupos independientes. Ver `VarianteEnsayo`.
   */
  varianteEnsayo?: VarianteEnsayo;
  calidad?: RiesgoSesgo;
  /**
   * Herramienta de riesgo de sesgo con la que se evaluó ESTE estudio (Paso 6). Ausente = la herramienta
   * por defecto del proyecto. Se fija al responder el primer dominio, para que cambiar el valor por
   * defecto después no reinterprete respuestas ya dadas (los ids de dominio se repiten entre
   * herramientas). Una revisión con ensayos y cohortes necesita herramientas distintas: antes había una
   * sola para todo el proyecto y una cohorte recibía las preguntas de aleatorización de RoB 2.0.
   */
  herramientaCalidad?: CalidadTool;
  /**
   * Justificación del riesgo de sesgo de ESTE estudio (Handbook cap. 8: registrar el sustento de cada juicio, p. ej. «sobres
   * opacos sellados, p. 3»). Texto libre y opcional; es un registro interno (no se imprime en el manuscrito) y viaja en el respaldo JSON.
   */
  notaCalidad?: string;
  evidencia?: CertezaGrade;
  origenRegistroId?: number;
  /** Si el efecto se calculó desde datos crudos, aquí quedan guardados los números originales. */
  crudos?: DatosCrudos;
  /**
   * Ensayo con más de dos grupos cuyo grupo CONTROL se comparte con otras filas de la misma revisión (Cochrane Handbook
   * 6.2.9 y 23.3.4): el número de comparaciones que lo comparten. Con k ≥ 2, el efecto se calcula dividiendo el control
   * en k partes iguales (eventos y total en los binarios; solo el N en los continuos), para no contar dos veces a los
   * mismos participantes. Los `crudos` guardan los números ORIGINALES; la división se aplica al calcular.
   */
  dividirControlEn?: number;
  /**
   * Qué dato de este estudio se obtuvo por conversión y cómo (p. ej. «DE imputada desde el IC 95 % de la media»; «media y DE
   * estimadas desde la mediana y los cuartiles, Wan 2014»). Texto libre y opcional: Métodos lo declara, porque un número
   * convertido no es un dato medido (Cochrane Handbook cap. 6).
   */
  conversionDeDatos?: string;
  /** Juicio por dominio de la herramienta de riesgo de sesgo (clave = id del dominio). El veredicto global se deriva de aquí. */
  dominiosCalidad?: Record<string, NivelDominio>;
  /**
   * Segunda extracción independiente (valor/IC) de un segundo revisor, para comparar contra la
   * extracción principal. Práctica recomendada por Cochrane Handbook y AMSTAR-2 ítem 6 (extracción
   * por duplicado) para minimizar errores de trascripción — la causa más común de errores en
   * meta-análisis, más que los errores de fórmula. `revelado` controla el cegado real: mientras es
   * false/ausente, la extracción principal (valor/icInferior/icSuperior de este mismo `Estudio`) se
   * oculta en la tabla del Paso 5 para que el segundo revisor no la vea mientras anota la suya —
   * antes solo se ocultaba de palabra en el texto de ayuda, sin aplicarse de verdad.
   */
  verificacion?: { valor: string; icInferior: string; icSuperior: string; revelado?: boolean };
  /**
   * Valores de los campos de extracción personalizados (Fase 1: plantillas configurables), keyed by
   * `CampoPersonalizado.id`. El set fijo de columnas de arriba (diseño, n, medida, efecto) cubre lo
   * mínimo para el meta-análisis; esto cubre lo demás que una revisión suele extraer (dosis, entorno,
   * umbral diagnóstico, etc.) según el tipo de pregunta — ver `plantillas-extraccion.ts`.
   */
  datosPersonalizados?: Record<string, string>;
  /**
   * Tabla 2×2 de precisión diagnóstica (verdaderos/falsos positivos/negativos), independiente de
   * `medida`/`valor`/`crudos` — un estudio de precisión diagnóstica no tiene un solo efecto, tiene un
   * par sensibilidad/especificidad. Alimenta el análisis SROC (Fase 4) en el Paso 7; ver `statistics.ts`
   * (`calcularSROC`). Vacío/parcial = el estudio no participa del SROC (no es un error).
   */
  datosDiagnosticos?: { vp: string; fp: string; fn: string; vn: string };
}

export type CalidadTool =
  | "RoB 2.0 (Cochrane)"
  | "ROBINS-I"
  | "Newcastle-Ottawa"
  | "QUADAS-2"
  | "AMSTAR 2"
  | "JBI";

/** Medida de efecto soportada por el meta-análisis en red (Fase 4). */
export type MedidaRed = "OR" | "DM";

/**
 * Un brazo (grupo de tratamiento) dentro de un estudio de red. A diferencia del modelo de un solo
 * efecto por estudio (`Estudio`), un estudio de red puede tener ≥2 brazos — cada uno es un tratamiento
 * distinto, no necesariamente "tratamiento vs control". Binario: `eventos`/`total`. Continuo:
 * `media`/`de`/`n`. Solo se llenan los campos de la medida elegida para el estudio.
 */
export interface BrazoRed {
  tratamiento: string;
  eventos?: string;
  total?: string;
  media?: string;
  de?: string;
  n?: string;
}

/**
 * Un estudio dentro de la red (meta-análisis en red, Fase 4). Independiente de `Estudio` — vive en su
 * propia lista (`estudiosRed` en el wizard), porque su forma (número variable de brazos) no encaja en
 * el modelo "un estudio, un efecto, un IC" del resto de la app.
 */
export interface EstudioRed {
  id: number;
  etiqueta: string;
  /**
   * Registro del Tamizaje (incluido en texto completo) del que sale la referencia bibliográfica de este estudio: autores,
   * título, revista y DOI. Sin vínculo, la referencia es solo la etiqueta.
   */
  origenRegistroId?: number;
  medida: MedidaRed;
  brazos: BrazoRed[];
}
