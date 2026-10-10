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

## Cómo ejecutarlo (paso a paso, sin experiencia previa)

1. **Instala Node.js 22 o superior** desde https://nodejs.org (es el programa que ejecuta el código; la versión «LTS» sirve).
2. **Descarga el código.** En esta página de GitHub, botón verde **Code → Download ZIP**, y descomprímelo. (Si usas git: `git clone https://github.com/altfernovamd-eng/anova-motor.git`.)
3. **Abre una terminal dentro de esa carpeta** (en Mac: Terminal; en Windows: PowerShell) y escribe estos comandos, uno por uno:

```bash
npm install                 # descarga las herramientas de desarrollo (tsx, typescript); no hay librerías de ejecución
npm test                    # pruebas del motor + comparación con los valores de R guardados (tarda un par de minutos)
npm run ejemplo             # un meta-análisis pequeño de ejemplo (ver abajo)
npm run validacion:masiva   # opcional: comparación masiva (≈375.000 comprobaciones)
```

No hace falta tener R: los resultados de R ya están guardados en `validacion/referencias/`. R solo se necesita para volver a generarlos.

## Ejemplo

`ejemplo.mts` calcula el OR combinado de tres estudios (efectos aleatorios, estimador REML):

```ts
import { calcularMetaAnalisis } from "./src/statistics.ts";
import type { Estudio } from "./src/tipos.ts";

const estudio = (id: number, valor: string, icInferior: string, icSuperior: string): Estudio =>
  ({ id, autor: `Estudio ${id}`, año: "2020", diseño: "ECA", n: "200", medida: "OR", valor, icInferior, icSuperior });

const resultado = calcularMetaAnalisis(
  [estudio(1, "0.80", "0.50", "1.20"), estudio(2, "0.65", "0.40", "1.05"), estudio(3, "0.90", "0.60", "1.35")],
  "OR",
  "REML",
);
console.log(resultado?.aleatorio);      // OR combinado y su IC 95 %
console.log(resultado?.heterogeneidad); // Q, I², τ², ...
```

`npm run ejemplo` muestra el efecto combinado con su intervalo de confianza (`aleatorio`) y las medidas de heterogeneidad (`heterogeneidad`: Q, I², τ²).
Para usar el motor en tu propio proyecto TypeScript, copia la carpeta `src/` e impórtala como en el ejemplo (respetando la licencia AGPL-3.0).
La entrada es el registro de estudio de aNova (campos numéricos como texto, tal como salen de un formulario).

## Ejemplo con datos reales y comparación con R

`ejemplos/ejemplo-bcg.mts` analiza un meta-análisis real de 13 ensayos de la vacuna BCG contra la tuberculosis (datos `dat.bcg` del paquete `metadat`
de R; Colditz et al., 1994) a partir de un CSV con las tablas 2×2 (`ejemplos/bcg.csv`):

```bash
npm run ejemplo:bcg                  # OR combinado ≈ 0,4746, IC 95 % [0,3296; 0,6835], τ² ≈ 0,3378
Rscript ejemplos/comparar-bcg.R      # opcional, si tienes R con metafor y metadat: debe dar los mismos números
```

Para analizar **tus propios datos**, copia `ejemplo-bcg.mts` y cambia el CSV: cada fila es un estudio con eventos y no-eventos de tratados y controles.

## Diferencias de convención con R (documentadas, no son errores)

Si comparas con R y ves una diferencia, probablemente sea una de estas (detalle y comprobaciones en `docs/VALIDACION.md`):

- **I².** Este motor muestra el I² clásico de Higgins, (Q − gl)/Q. `metafor` muestra otro basado en τ², τ²/(τ² + s²); coinciden solo si τ² se estima por DerSimonian-Laird. Para BCG: 92,6 % frente a 92,1 %.
- **Intervalo de predicción.** Aquí se usa t con k − 2 grados de libertad (IntHout et al. 2014; en `meta`, `method.predict = "HTS"`). `predict()` de `metafor` usa por defecto la normal y `meta` ≥ 8 usa t con k − 1; con pocos estudios los intervalos salen distintos (BCG: [0,124; 1,818] aquí, [0,126; 1,794] con el valor por defecto de `meta`).
- **Mantel-Haenszel, diferencia de riesgos.** Se sigue a RevMan 5 y al Cochrane Handbook; `rma.mh` usa otra varianza.
- **Razón de medias.** `escalc("ROM")` añade una corrección de sesgo de primer orden; aquí se usa la fórmula del Handbook.
- **Q de `madauni`.** Usa pesos de efectos aleatorios; aquí se usa la Q de Cochran con pesos de efectos fijos (igual que `metafor` y `meta`).
- **Celdas en cero en la red de OR.** Se sigue la regla de `meta::pairwise` y `netmeta`; los estudios sin ningún contraste informativo se descartan y se declaran.

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
