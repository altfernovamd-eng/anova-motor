# aNova motor — motor estadístico de meta-análisis

Motor estadístico en TypeScript (sin dependencias de ejecución) que usa la aplicación aNova Meta Tool. Cada método se compara número a número
con R (`metafor`, `meta`, `netmeta`, `mada`, `irr`, `clubSandwich`) usando R solo como caja negra.

> **English summary.** Dependency-free TypeScript engine for meta-analysis (fixed/random effects, DerSimonian-Laird, Paule-Mandel, REML,
> Hartung-Knapp, prediction intervals, meta-regression, publication-bias tests, trim-and-fill, Mantel-Haenszel, Peto, network meta-analysis,
> bivariate DTA/SROC, multilevel models with robust variance, binomial/hypergeometric GLMMs, kappa). Validated against R packages with
> ≈375,000 randomized comparisons (see `docs/VALIDACION.md`). **Status: pre-release. License: GNU AGPL-3.0-only (commercial license available from the copyright holder).** The code was written by an AI
> assistant (Claude, Anthropic) under the author's direction; see `docs/PROCEDENCIA.md`.

## Contenido

| Módulo (`src/`) | Qué hace |
|---|---|
| `statistics.ts` | Medidas de efecto (OR, RR, RAR, DM, DME, tasas, correlaciones, HR), efectos fijos/aleatorios (DL, PM, REML), Hartung-Knapp, intervalo de predicción, meta-regresión, Egger/Begg/Peters, trim-and-fill, Mantel-Haenszel, Peto, subgrupos, acumulado, influencia, Galbraith, SROC univariado |
| `red.ts` | Meta-análisis en red (método gráfico-teórico), P-scores, node-splitting, descomposición de Q |
| `sroc-bivariado.ts` | Modelo bivariado de precisión diagnóstica (Reitsma) |
| `multinivel.ts` | Efectos dependientes: modelo multinivel con varianza robusta CR2 |
| `glmm.ts` | Modelos mixtos logísticos (hipergeométrico-normal y binomial-normal) |
| `kappa.ts` | κ de Cohen y de Fleiss |
| `conversiones-datos.ts` | Conversiones (DE desde EE/IC, mediana y cuartiles, etc.) |
| `efecto-crudos.ts` | Efecto a partir de datos crudos (2×2, medias, tasas, cruzado…), cero eventos y verificación contra un segundo extractor |
| `funciones-especiales.ts`, `numerica-compartida.ts` | Funciones especiales y numéricas propias |
| `tipos.ts` | Tipos de datos (subconjunto del modelo de la aplicación) |

## Uso

```bash
npm install          # solo herramientas de desarrollo (tsx, typescript)
npm test             # pruebas + validación contra los valores de R guardados
npm run validacion:masiva   # comparación masiva (≈375.000 comprobaciones)
```

```ts
import { calcularMetaAnalisis } from "./src/statistics.ts";

const e = (id: number, valor: string, lo: string, hi: string) =>
  ({ id, autor: `E${id}`, año: "2020", diseño: "ECA", n: "200", medida: "OR", valor, icInferior: lo, icSuperior: hi }) as const;
const r = calcularMetaAnalisis([e(1, "0.80", "0.50", "1.20"), e(2, "0.65", "0.40", "1.05"), e(3, "0.90", "0.60", "1.35")] as never, "OR", "REML");
console.log(r?.aleatorio.estimacion, r?.i2, r?.tau2);
```

La entrada es el registro de estudio de aNova (campos numéricos como texto, tal como salen de un formulario).

## Validación

Los valores de R están guardados en `validacion/referencias/*.json` y se regeneran con los scripts `validacion/generar-referencia*.R`
(requieren R con `metafor`, `meta`, `netmeta`, `mada`, `irr`, `clubSandwich`, `jsonlite`). Resultados y hallazgos: `docs/VALIDACION.md`.
Alcance (qué cubre y qué no): `docs/ALCANCE_MOTOR.md`.

## Contribuciones

Este proyecto **no acepta contribuciones externas de código** (pull requests): el titular conserva todos los derechos para poder ofrecer también una licencia comercial. Los errores y las dudas se pueden reportar como *issues*.

## Procedencia y licencia

Quién escribió el código, de dónde salen los métodos y qué código GPL se consultó durante el desarrollo: `docs/PROCEDENCIA.md`.
Licencia: **GNU AGPL-3.0 solo** (`LICENSE`; identificador SPDX `AGPL-3.0-only`). Copyright © 2026 Fernando Ramírez Nova. Quien incorpore este motor a un producto o servicio debe ofrecer su código completo bajo la misma licencia. Para usarlo dentro de un producto cerrado existe una licencia comercial aparte: altfernovamd@gmail.com. Datos de prueba: resultados publicados, citados en `docs/PROCEDENCIA.md` (sección 5).

*Generado por `scripts/exportar-motor.mts` de la aplicación aNova. No editar a mano: se regenera.*
