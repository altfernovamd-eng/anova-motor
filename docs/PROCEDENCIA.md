# Procedencia del motor estadístico

> *Nota: este documento procede del repositorio de la aplicación aNova Meta Tool. Las menciones a «Paso N», pantallas o paneles se refieren a su interfaz; los módulos del motor están en `src/`.*


> Documento técnico. Recoge **hechos comprobados** sobre de dónde sale cada parte del motor. No contiene valoración jurídica.
> El motor se publica bajo AGPL-3.0-only (ver `LICENSE`).

## 1. Resumen

1. **Autoría y asistencia de IA.** El código lo escribió un modelo de lenguaje (Claude, Anthropic) bajo la dirección y revisión del autor, que es
   responsable de su contenido. La verificación es independiente: comparación número a número contra R (`npm run validacion`).
2. **Fuentes de los métodos.** Las fórmulas y los métodos vienen de artículos y manuales publicados (columna «Referencia» de la sección 2).
3. **Código de terceros consultado.** Durante el desarrollo se leyó el código fuente de varios paquetes de R con licencia GPL (netmeta, mada, irr y,
   al investigar diferencias, también meta y metafor) para confirmar detalles de algunos métodos. No consta que se copiara código. Varias partes
   se reescribieron después desde las publicaciones; ninguna de esas reescrituras es independiente en sentido estricto (sección 3).
4. **Sin librerías de terceros en el motor.** Solo importa tipos propios. Las funciones numéricas (gamma, beta, ln Γ) están escritas desde las
   fórmulas matemáticas y comprobadas contra R (sección 4).

## 2. Método por método

«Código de terceros consultado» indica lo que consta en el registro de desarrollo; **«no consta» no prueba que no ocurriera**.

| Módulo / función | Método | Referencia citada en el código | Código de terceros consultado | Comprobación |
|---|---|---|---|---|
| `statistics.ts` · `efectoDesde2x2`, `efectoDesdeMedias`, `efectoDesdeProporcion` | OR, RR, RAR, DM, proporción | Cochrane Handbook 6.4.3.1 (corrección de continuidad) | no consta | contra `metafor::escalc` (`validacion`) |
| `efectoDesdeMediasEstandarizada` | g de Hedges | Hedges 1982 | no consta | contra R |
| `calcularMetaAnalisis` | Efectos fijos, DerSimonian-Laird, Paule-Mandel, REML; Hartung-Knapp; intervalo de predicción | Paule & Mandel 1982; DerSimonian & Kacker 2007; Viechtbauer 2005; Hartung & Knapp 2001; IntHout 2014 | `metafor::rma` leído con `deparse` el 5 oct 2026, para entender diferencias | contra `metafor`/`meta` |
| `confianzaHeterogeneidad` | IC de τ², τ, I² (Q-profile) | Viechtbauer 2007 | `metafor::rma` (ídem) | contra R |
| `ajustarMetaRegresion`, `metaRegresionMulti` | Meta-regresión de efectos mixtos REML | «réplica de `metafor::rma(mods=…)`» | `metafor::rma` (ídem) | 60 conjuntos contra `metafor` |
| `testEgger`, `testBegg`, `testPeters` | Asimetría del embudo | Egger 1997; Begg & Mazumdar 1994; Peters et al. 2006 | `meta::metabias` leído con `deparse` el 5 oct 2026 | contra `meta` |
| `metaAnalisisMH`, `metaAnalisisPeto` | Mantel-Haenszel y Peto | Cochrane Handbook §10.4 (Deeks & Higgins); Robins, Breslow & Greenland 1986; Breslow & Greenland 1986; Rothman et al. 2008 | `metafor::rma.mh` leído con `deparse` el 5 oct 2026 | contra `metafor` (diferencia de riesgos: convención RevMan 5) |
| `trimAndFill` | Trim-and-fill (estimador L0) | Duval & Tweedie 2000 | ejecutado como caja negra (`metafor::trimfill`) | contra `metafor` |
| `baujat`, `analisisSubgrupos`, `metaAnalisisAcumulativo`, `analisisSensibilidad` | Baujat, subgrupos, acumulativo, dejar uno fuera | Baujat et al. 2002 | no consta | contra R |
| `efectoDeDisenoCluster`, `efectoDesdeCruzado` | Ajuste de ensayos clúster y cruzados | Cochrane Handbook cap. 23 (texto oficial leído) | no aplica | contra el Handbook (`test:cluster-cruzado`); sin paquete de R de referencia |
| `moses-littenberg.ts` · `calcularSROC` (parte univariada) | Curva SROC y regresión de Moses-Littenberg, DOR, intervalo de Wilson | Moses, Shapiro & Littenberg 1993; Littenberg & Moses 1993; Glas y cols. 2003; Wilson 1927; Newcombe 1998 | **`mada::madauni` (`madauni.R`) leído el 24 sept 2026**, antes de la reescritura. **Reescrito el 7 oct 2026** desde los artículos (otra organización, regresión en forma cerrada, otros nombres); `mada` solo como caja negra | contra `mada` y 60 conjuntos (`validacion`), `benchmark:sroc`, `test:moses` |
| `sroc-bivariado.ts` | Modelo bivariado de precisión diagnóstica: ML binomial exacta (Gauss-Hermite adaptativa) y REML normal; curva SROC, AUC, regiones | Reitsma y cols. 2005; Chu y Cole 2006; Harbord y cols. 2007; Rutter y Gatsonis 2001; Golub y Welsch 1969; Cochrane Handbook de precisión diagnóstica, cap. 10 | **no se abrió código** de `mada`, `lme4` ni `mvmeta`; usados solo como cajas negras (salidas) | contra implementación independiente en R, fuerza bruta, `lme4`, `mada` (31 conjuntos, `test:sroc-bivariado`) |
| `multinivel.ts` | Modelo CHE (REML) y varianza robusta CR2 con Satterthwaite | Pustejovsky y Tipton 2022; Hedges, Tipton y Johnson 2010; Tipton 2015; Pustejovsky y Tipton 2018; Bell y McCaffrey 2002 | **no se abrió código** de `metafor`, `clubSandwich` ni `robumeta`; usados solo como cajas negras (salidas) | contra `rma.mv` y `clubSandwich` (60 conjuntos, `test:multinivel`) |
| `glmm.ts` | GLMM: hipergeométrico-normal (OR) y binomial-normal (proporciones), máxima verosimilitud con Gauss-Hermite adaptativa | Stijnen, Hamza y Özdemir 2010; Jackson y cols. 2018; Golub y Welsch 1969 | **no se abrió código** de `metafor` ni `lme4`; cajas negras; `BiasedUrn::dFNCHypergeo` como caja negra para las probabilidades | contra implementación independiente en R, fuerza bruta y `rma.glmm` (80 conjuntos, `test:glmm`) |
| `red.ts` · `calcularNMA` | Meta-análisis en red (gráfico-teórico) | Rücker 2012; Rücker & Schwarzer 2014 (Stat Med 33:4353) | **`netmeta`: `netmeta.R`, `nma.ruecker.R`, `multiarm.R` leídos el 24 sept 2026**; `nma_ruecker` y `netmeta` otra vez el 5 oct. **Reescrito el 7 oct 2026** desde los artículos (Laplaciano con un nodo a tierra y Cholesky, pesos multibrazo desde la precisión de los contrastes, τ² por momentos derivado) | contra `netmeta` (Stowe2010 y 40 redes), `test:red` y una prueba diferencial con la versión anterior (792 redes, 58.932 comparaciones, dif. máx. 4e-12) |
| `calcularPScores` | P-scores | Rücker & Schwarzer 2015 (BMC MRM 15:58) | **`netrank.R`, `rankogram.default.R`** (sept 2026). **Reescrito el 7 oct 2026** junto con `red.ts` | contra `netmeta` |
| node-splitting (back-calculation) y SIDDE | Inconsistencia local | König, Krahn & Binder 2013 (Stat Med 32:5414); Efthimiou y cols. 2019 (Stat Med 38:2992); Dias y cols. 2010 | **`netsplit.R`, `sidde.R`**; `netsplit` otra vez el 5 oct. **Reescrito el 7 oct 2026** con `red.ts` (la convención de agrupar la evidencia directa con las varianzas propias, sin la reducción multibrazo, se fijó comparando salidas) | contra `netmeta` |
| descomposición de Q | Q dentro y entre diseños | Krahn, Binder & König 2013 (BMC MRM 13:35, ecs. 3, 4, 8 y 10) | **`decomp.design.R`, `decomp.tau.R`** (sept 2026). **Reescrito el 7 oct 2026** con `red.ts` (mínimos cuadrados generalizados por diseño con Cholesky) | contra `netmeta` |
| `kappa.ts` | κ de Cohen y de Fleiss | Cohen 1960; Fleiss 1971; Fleiss, Cohen & Everitt 1969 y Fleiss, Nee & Landis 1979 (varianza bajo H₀); Landis & Koch 1977 | **`irr`: `kappa2.R` y `kappam.fleiss.R` leídos en sept 2026**. **Reescrito el 6 oct 2026** desde las fórmulas desarrolladas (otra estructura y otros nombres; `irr` solo como caja negra) | contra `irr` (dataset `diagnoses` y 177 conjuntos, `test:kappa-r`) |
| `funciones-especiales.ts` | ln Γ, gamma y beta incompletas | A&S 6.1.40, 6.5.29, 26.5.8; DLMF 8.9.2, 8.17.22; Lentz 1976 | no (reescrito el 6 oct 2026, sección 4) | contra R, 1.109 comprobaciones |
| `statistics.ts` · `efectoDesdeTasas` | Razón y diferencia de tasas de incidencia (RTI, DTI) | Cochrane Handbook 6.7 y 10.8 | no se abrió código; `metafor::escalc` (IRR, IRD) como caja negra (salidas) | contra R (40 conjuntos) |
| `efectoDesdeRazonDeMedias` | Razón de medias (RM) | Handbook 6.5.1.3; Friedrich, Adhikari y Beyene 2008 | no se abrió código | fórmula; `escalc("ROM")` difiere por una corrección de sesgo (convención) |
| `efectoDesdeCorrelacion` | Correlaciones, z de Fisher | Fisher 1915; `escalc("ZCOR")` como caja negra | no se abrió código | contra R (40 conjuntos) |
| `efectoDesdeOEV` | HR desde O−E y V (tiempo hasta el evento) | Handbook 6.8.2 y 10.9; Peto | no se abrió código | identidad exacta con Σ(O−E)/ΣV |
| `diagnosticosInfluencia`, `galbraith` | Diagnóstico de influencia; gráfico radial | Viechtbauer y Cheung 2010; Galbraith 1988 | `influence()` y `radial()` como caja negra, y la **documentación** de `metafor`; las definiciones se fijaron comparando candidatas con su salida | contra R, ≈1e-10 en 835 estudios |
| `conversiones-datos.ts` | DE desde EE, IC, t o P; combinar grupos; DE del cambio; mediana y cuartiles; OR ↔ DME | Handbook cap. 6 y 10.6; Wan y cols. 2014 | `meta::metamean` como caja negra | `test:conversiones` (1.599 comprobaciones) |
| `testBegg` (empates) | τ-b y varianza de Kendall con empates | Kendall (1970) | no se abrió código; `ranktest` como caja negra | contra R |
| `statistics.ts` · `divisorControl` | División del grupo control compartido | Handbook 6.2.9 y 23.3.4 | no aplica | `test:ensayos-repetidos` en la aplicación |
| `statistics.ts` · `evaluarGrade` | Reglas de decisión de GRADE | Guías GRADE y Cochrane Handbook | no aplica | revisión manual; no es comparable con un paquete |

## 3. Código de terceros consultado y reescrituras

- `netmeta`, `mada`, `irr`, `meta` y `metafor` están publicados bajo GPL (≥ 2) o GPL-2 (según sus archivos `DESCRIPTION`).
- **Red (`red.ts`).** La implementación pasó por un intermedio: un programa independiente en Python que reproducía el algoritmo de Rücker «fórmula por
  fórmula desde el artículo», confirmando detalles leyendo el código de `netmeta` (ajuste de estudios multibrazo, P-scores, node-splitting y
  descomposición de Q), y luego se llevó a `red.ts`.
- **`kappa.ts`.** Sus variables coincidían con las de `kappa2.R` de `irr` (se había escrito siguiéndolo de cerca, aunque en otro lenguaje). Se reescribió
  el 6 de octubre de 2026 desde las fórmulas publicadas en su forma desarrollada (una sola suma en vez del doble bucle, otros nombres, otra
  estructura). No es una reescritura limpia en sentido estricto: quien la hizo ya había visto el código de `irr`. `irr` se usó después solo como caja
  negra para compararla (177 conjuntos).
- **`red.ts` y regresión de Moses-Littenberg (7 de octubre de 2026).** Se partió de los artículos (ecuaciones citadas en las cabeceras de los módulos),
  con otra organización del código, otras rutas de cálculo y otros nombres. No es una reescritura limpia en sentido estricto: la hizo la misma sesión de IA
  que había visto el código de `netmeta` y de `mada`, y no se puede descartar que el entrenamiento del modelo incluya ese código. La respalda que el resultado
  coincide con la versión anterior en 792 redes aleatorias (58.932 comparaciones, diferencia máxima relativa 4e-12; comparación hecha una vez con una copia
  temporal de la versión anterior que ya no existe, así que no puede repetirse) y con `netmeta` y `mada` en las validaciones, y que las diferencias de estructura
  se pueden comprobar en el código.
- **Meta-regresión, pruebas de asimetría del embudo, Mantel-Haenszel e IC de τ².** Se escribieron originalmente a partir de las publicaciones; la tabla de la
  sección 2 indica qué función de R se leyó (`rma`, `metabias`, `rma.mh`) para entender diferencias. No se han reescrito después.
- **Piezas añadidas el 6 de octubre de 2026** (tasas, razón de medias, correlaciones, O−E y V, influencia, conversiones) y los módulos bivariado,
  multinivel y GLMM: no se abrió ni se leyó el código fuente de ningún paquete; se usaron los artículos, el Handbook, la documentación (páginas de ayuda) de
  `metafor` y `meta`, y sus resultados como caja negra. No se puede descartar que el modelo de lenguaje conozca esas fórmulas por haber visto código de estos
  paquetes durante su entrenamiento; no hay forma de comprobarlo.
- **La validación contra R trata a R como caja negra** (solo compara salidas); la lectura de código fue anterior o aparte.

## 4. Funciones numéricas reescritas (6 de octubre de 2026)

`statistics.ts`, `red.ts` y `kappa.ts` tenían tres copias de las funciones de gamma y beta incompletas, etiquetadas como «algoritmo de
*Numerical Recipes*» (cuyo código tiene una licencia restrictiva). Se sustituyeron por un único módulo, `funciones-especiales.ts`,
escrito desde las fórmulas matemáticas (serie de Stirling con números de Bernoulli, serie de potencias y fracciones continuas
evaluadas con el método de Lentz modificado; fuentes en la cabecera del módulo).

- Comprobado contra R (`pchisq`, `pt`, `qt`, `pnorm`, `lgamma`): 378 puntos de χ², 240 de la t, 18 de la normal, 22 de ln Γ y 15 cuantiles,
  desde p ≈ 1 hasta p ≈ 1e-290. Peor error relativo: 3,9e-12. `npm run test:funciones-especiales`.
- **Mejora de paso:** la versión anterior devolvía 0 en colas extremas (p < 1e-15). No cambiaba ninguna conclusión (esos p se muestran como
  «<0,001»), pero ahora son exactos.
- Sin regresiones: `npm run validacion` dio entonces 68.443 comparaciones y 0 fuera de tolerancia; al ampliar el alcance (§ 2 de `docs/ALCANCE_MOTOR.md`) son 237.089 comparaciones en 995 grupos, también con 0 fuera.

## 5. Datos usados en las pruebas

| Datos | Origen | Licencia del paquete que los distribuye |
|---|---|---|
| `dat.bcg` (13 estudios, vacuna BCG) | paquete `metadat` | GPL (≥ 2) |
| `ThirdWave` | paquete `dmetar` | no verificada |
| `AuditC` (14 estudios, Kriston et al. 2008) | paquete `mada` | GPL-2 |
| `Stowe2010` (29 estudios, Parkinson) | paquete `netmeta` | GPL (≥ 2) |
| `diagnoses` (30 sujetos, 6 evaluadores) | paquete `irr` | GPL (≥ 2) |
| 705 conjuntos sintéticos | generados con semilla fija por los scripts propios de `validacion/*.R` | propios |

Los datos son resultados publicados de estudios (hechos), y el repositorio los incluye como cifras en los scripts, no los paquetes.
La licencia individual de cada conjunto de datos no se verificó; lo que se comprobó es la licencia del paquete.

## 6. Declaración de uso de IA (para el artículo)

> El motor estadístico se escribió en TypeScript con la asistencia de un modelo de lenguaje (Claude, Anthropic) bajo la dirección y revisión del autor, que es
> responsable de su contenido. Las fórmulas se tomaron de las publicaciones citadas. Para confirmar detalles de algunos métodos se consultó el código fuente de
> paquetes de R (netmeta, mada, irr, meta, metafor), y varias partes se reescribieron después desde las publicaciones. La corrección se verificó comparando los
> resultados con esos paquetes en cientos de conjuntos de datos y redes sintéticos y en conjuntos publicados (más de 370.000 comparaciones numéricas en la
> validación masiva; 255.661 en la validación estándar).
