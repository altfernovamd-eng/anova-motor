# Validación estadística de aNova

> *Nota: este documento procede del repositorio de la aplicación aNova Meta Tool. Las menciones a «Paso N», pantallas o paneles se refieren a su interfaz; los módulos del motor están en `src/`.*


**¿Los números que produce aNova son confiables?** Esta es la pregunta más importante que puede
hacerse quien va a publicar con una herramienta de meta-análisis, porque un solo resultado mal
calculado invalida un artículo y la revista lo rechaza. Este documento responde con evidencia, no
con promesas: aNova reproduce **exactamente** la salida del estándar de referencia mundial.

> **Resumen en una frase:** para el meta-análisis estándar, aNova da los **mismos números, hasta el
> último decimal**, que el paquete `meta` de R — la herramienta que usan y verifican los
> investigadores y revisores de todo el mundo.

---

## 1. ¿Contra qué se compara y por qué?

El paquete **`meta`** (y su hermano `metafor`) para el lenguaje **R** es el *gold standard* de facto
del meta-análisis: es de código abierto, está auditado por la comunidad estadística desde hace más de
una década, y es lo que un revisor de revista usará para recalcular tus resultados. Si aNova coincide
con `meta`, entonces un artículo hecho con aNova es **defendible ante cualquier revisor**.

No basta con "revisamos que las fórmulas estén bien escritas" (eso es revisar la receta). La prueba de
verdad es **hornear el pastel y comprobar que sabe igual que el del maestro**: correr un caso conocido
por ambas herramientas y comparar cifra a cifra.

## 2. ¿Cómo se hizo la prueba?

- **Dataset canónico:** `ThirdWave` (18 estudios), del libro de referencia *Doing Meta-Analysis in R*
  (Harrer, Cuijpers, Furukawa & Ebert). Los datos se tomaron del archivo original del paquete
  `dmetar`, con precisión numérica completa — **no se transcribieron cifras a mano**.
- **Motor evaluado:** el **código de producción real** de aNova (`src/statistics.ts`),
  el mismo que corre dentro de la aplicación. No es una reimplementación hecha para la prueba.
- **Referencia:** la salida **publicada** de `meta` con
  `metagen(TE, seTE, method.tau = "REML", hakn = TRUE)` — es decir, el modelo de efectos aleatorios
  con el estimador **REML** de la heterogeneidad y el ajuste **Knapp-Hartung** del intervalo de
  confianza, que son precisamente las dos partes "difíciles" donde una implementación floja se
  desviaría.

## 3. Resultado: coincidencia exacta

| Cantidad | aNova | R `meta` (publicado) | ¿Coincide? |
|---|---|---|:---:|
| Efecto combinado aleatorio (SMD) | 0.5771 | 0.5771 | ✓ |
| IC 95 % (Knapp-Hartung) | [0.3782; 0.7760] | [0.3782; 0.7760] | ✓ |
| Estadístico *t* (gl = 17) | 6.12 | 6.12 | ✓ |
| τ² — varianza entre estudios (REML) | 0.0820 | 0.0820 | ✓ |
| τ (REML) | 0.2863 | 0.2863 | ✓ |
| I² — heterogeneidad | 62.6 % | 62.6 % | ✓ |
| Q de Cochran | 45.50 | 45.50 | ✓ |
| Grados de libertad | 17 | 17 | ✓ |
| p-valor de Q | 0.0002 | 0.0002 | ✓ |

**Las 10 cantidades coinciden a todos los decimales publicados.** Esto incluye el cálculo iterativo de
REML y el ajuste con distribución *t* de Knapp-Hartung, no solo la aritmética básica.

### 3.2. Segundo método validado: Mantel-Haenszel (datos binarios)

El método de **Mantel-Haenszel** (recomendado cuando los eventos son raros) se validó contra
**`metafor::rma.mh`** sobre el dataset canónico **`dat.bcg`** (13 estudios de la vacuna BCG), para riesgo
relativo (RR) y odds ratio (OR):

| Cantidad | aNova (RR) | metafor (RR) | aNova (OR) | metafor (OR) | ¿Coincide? |
|---|---|---|---|---|:---:|
| log-estimación | −0.4537 | −0.4537 | −0.4734 | −0.4734 | ✓ |
| Error estándar (log) | 0.0393 | 0.0393 | 0.0410 | 0.0410 | ✓ |
| z | −11.5338 | −11.5338 | −11.5444 | −11.5444 | ✓ |
| Efecto (escala natural) | 0.6353 | 0.6353 | 0.6229 | 0.6229 | ✓ |
| IC 95 % inferior | 0.5881 | 0.5881 | 0.5748 | 0.5748 | ✓ |
| IC 95 % superior | 0.6862 | 0.6862 | 0.6750 | 0.6750 | ✓ |
| Q de heterogeneidad | 152.5676 | 152.5676 | 163.9426 | 163.9426 | ✓ |
| I² | 92.13 % | 92.13 % | 92.68 % | 92.68 % | ✓ |

Coincidencia exacta en las **20 cantidades**, incluyendo la varianza de Robins-Breslow-Greenland y la
convención exacta de metafor para la Q (centrada en el estimador de Mantel-Haenszel).

### 3.3. Concordancia entre revisores: kappa (κ)

La fiabilidad inter-observador del tamizaje (**κ de Cohen** para 2 revisores, **κ de Fleiss** para ≥3) se
validó contra el paquete **`irr`** de R sobre su dataset `diagnoses` (30 sujetos, 6 evaluadores):

| Coeficiente | Cantidad | aNova | R `irr` (publicado) | ¿Coincide? |
|---|---|---|---|:---:|
| **Cohen** (`kappa2`, 2 evaluadores) | κ | 0.631 | 0.631 | ✓ |
| | z | 7.56 | 7.56 | ✓ |
| **Fleiss** (`kappam.fleiss`, 6 evaluadores) | κ | 0.43 | 0.43 | ✓ |
| | z | 17.7 | 17.7 | ✓ |

Coinciden el coeficiente κ y el estadístico z (que valida la fórmula de la varianza). El p-valor de estos
casos es del orden de 10⁻¹⁴ (prácticamente cero) en ambos.

**Prueba masiva (6 de octubre de 2026, tras reescribir `kappa.ts` desde las fórmulas publicadas).** 87 conjuntos de Cohen
(2 evaluadores, 2 a 6 categorías, 8 a 250 sujetos) y 90 de Fleiss (3 a 8 evaluadores) sintéticos con semilla fija, más tres casos
límite (acuerdo perfecto, peor que el azar), comparados con `irr` 0.85 usado como caja negra (`npm run test:kappa-r`, referencia
generada por `validacion/generar-referencia-kappa.R`). Peor diferencia: κ 2,5e-15 (absoluta), z 3,7e-14 y p 3,4e-12
(relativas). Además se comprueban a mano un caso 2×2 de Cohen (κ = 0,4; SE = √0,0192) y uno de Fleiss (κ = 1/3), y propiedades
(acuerdo perfecto, una sola categoría, intercambiar evaluadores). Verificada con 5 mutaciones de las fórmulas.

Dos observaciones. **(1)** El p que imprime `irr` pierde cifras cuando es muy pequeño (p < 1e-8: error relativo de hasta 2e-2 frente a
`2·pnorm(−|z|)` de R); el de aNova coincide con esa cola exacta, así que la comparación se hace contra ella. **(2)** El error estándar
que usa este cálculo es el que vale **bajo H₀: κ = 0** (Fleiss, Cohen y Everitt 1969; Fleiss, Nee y Landis 1979): sirve para el contraste
(z y p), **no** para construir un intervalo de confianza de κ (para eso hay otra varianza, p. ej. la de Gwet por linealización). aNova solo
muestra κ, su interpretación y el % de acuerdo; no muestra ni usa un intervalo.

### 3.4. Método de Peto (OR de un paso)

El método de **Peto** (odds ratio de un paso, para eventos raros con grupos balanceados) se validó contra
**`metafor::rma.peto`** sobre `dat.bcg`:

| Cantidad | aNova | metafor | ¿Coincide? |
|---|---|---|:---:|
| log-estimación | −0.4744 | −0.4744 | ✓ |
| Error estándar (log) | 0.0407 | 0.0407 | ✓ |
| z | −11.6689 | −11.6689 | ✓ |
| OR | 0.6222 | 0.6222 | ✓ |
| IC 95 % | [0.5746; 0.6738] | [0.5746; 0.6738] | ✓ |
| Q de heterogeneidad | 167.7302 | 167.7302 | ✓ |
| I² | 92.85 % | 92.85 % | ✓ |

Coincidencia exacta en las **10 cantidades** (diferencia observado−esperado y varianza hipergeométrica).

### 3.5. Diferencia de medias estandarizada: g de Hedges

El cálculo de la **g de Hedges** (SMD, con la corrección de sesgo de muestra pequeña) se validó reproduciendo
la fórmula exacta de **`metafor::escalc(measure="SMD")`** —tomada del código fuente de metafor— sobre el
dataset `dat.normand1999` (9 estudios). Coincide el efecto **g** y su varianza **v** de los 9 estudios; una
muestra:

| Estudio | g (aNova) | g (metafor) | v (aNova) | v (metafor) |
|---|---|---|---|---|
| Edinburgh | −0.3552 | −0.3552 | 0.0131 | 0.0131 |
| Orpington-Moderate | −2.3176 | −2.3176 | 0.0458 | 0.0458 |
| Umeå | −0.4246 | −0.4246 | 0.0149 | 0.0149 |
| … (9/9) | ✓ | ✓ | ✓ | ✓ |

Incluye el factor de corrección exacto `J = Γ(m/2)/(√(m/2)·Γ((m−1)/2))` y la varianza de Hedges (1982, ec. 8).
La referencia se computó con una implementación independiente de la función `lgamma` (no la de aNova), de modo
que la coincidencia comprueba de verdad la aritmética del motor.

### 3.6. Meta-análisis acumulativo

El meta-análisis acumulativo (ordenar por año y recalcular el efecto añadiendo un estudio a la vez) se validó
contra **`metafor::cumul`** sobre `dat.bcg` (riesgo relativo, modelo aleatorio REML, 13 pasos). **La estimación
combinada y su intervalo de confianza coinciden en los 13 pasos** (log-escala, 3 decimales); una muestra:

| Paso (k) | Estudio añadido | log-RR aNova | log-RR metafor | IC 95 % (aNova) |
|---|---|---|---|---|
| 1 | Aronson, 1948 | −0.889 | −0.889 | [−2.008; 0.229] |
| 3 | Stein & Aronson, 1953 | −0.972 | −0.972 | [−1.453; −0.491] |
| 13 | TPT Madras, 1980 | −0.715 | −0.715 | [−1.067; −0.362] |

> **Nota honesta sobre I²:** aNova reporta el I² **clásico de Higgins** `(Q−gl)/Q` (independiente del estimador),
> mientras que metafor reporta un I² basado en el estimador de τ² (`τ²/(τ²+s²)`). Coinciden cuando se usa
> DerSimonian-Laird, pero divergen un poco con REML. No es un error: son **dos definiciones válidas** del mismo
> concepto. Lo que importa del acumulativo —la evolución del efecto y su IC— coincide exactamente.

### 3.7. Test de Begg (correlación de rangos)

El test de rangos de **Begg & Mazumdar** para asimetría del funnel se validó contra **`metafor::ranktest`**
sobre `dat.bcg` (RR):

| Cantidad | aNova | metafor `ranktest` | ¿Coincide? |
|---|---|---|:---:|
| τ de Kendall | 0.0256 | 0.0256 | ✓ |
| p-valor (exacta) | 0.9524 | 0.9524 | ✓ |

Incluye el **p-valor exacto** de Kendall (distribución de inversiones, como `cor.test` de R para k pequeño sin
empates); para k>18 aNova usa la aproximación normal.

### 3.8. Intervalos de confianza de τ², τ e I² (método Q-profile)

Los IC del 95% de la heterogeneidad por el **método Q-profile** (Viechtbauer 2007) se validaron sobre ThirdWave
contra la salida del paquete **`meta`** (que usa Q-profile por defecto):

| Cantidad | aNova | R `meta` | ¿Coincide? |
|---|---|---|:---:|
| IC de τ² | [0.0295; 0.3533] | [0.0295; 0.3533] | ✓ |
| IC de τ | [0.1717; 0.5944] | [0.1717; 0.5944] | ✓ |

El IC de τ² (y por tanto el de τ) **coincide exactamente**.

> **Nota honesta sobre el IC de I²:** a partir del IC de τ², aNova deriva el de I² con la transformación
> Q-profile `I²=τ²/(τ²+s²)` — que es la que usa **`metafor::confint`**. El paquete `meta` (metagen), en cambio,
> calcula el IC de I² por el método test-based de Higgins-Thompson (vía el estadístico H), y da un intervalo algo
> más estrecho. Son **dos métodos válidos y publicados** para el mismo concepto; aNova usa el coherente con el
> Q-profile que pediste. El IC de τ², que es el que importa y no tiene ambigüedad, coincide con ambos.

### 3.9. Meta-regresión (moderadores múltiples y categóricos)

La meta-regresión de efectos mixtos (τ² por REML) se validó contra **`metafor::rma(mods=…, method="REML")`** sobre
`dat.bcg` (RR), con un moderador **continuo** (latitud) y uno **categórico** (método de asignación, 3 niveles):

| Modelo | Cantidad | aNova | metafor |
|---|---|---|---|
| ~ latitud | τ² residual | 0.0763 | 0.0763 ✓ |
| | pendiente (SE) | −0.0291 (0.0072) | −0.0291 (0.0072) ✓ |
| | QM (test moderador) | 16.3582 | 16.3582 ✓ |
| | R² | 75.63 % | 75.63 % ✓ |
| ~ asignación (categórico) | τ² residual | 0.3615 | 0.3615 ✓ |
| | coef. random / systematic | −0.4478 / 0.0890 | −0.4478 / 0.0890 ✓ |
| | QM (gl=2), p | 1.7675, 0.4132 | 1.7675, 0.4132 ✓ |

Coincidencia exacta en las 11 cantidades, incluyendo el álgebra matricial (diseño p×p), la codificación dummy de
la variable categórica y el test ómnibus QM.

> **Sobre la referencia:** la salida de metafor se reprodujo con una implementación independiente en Python/numpy
> del mismo algoritmo REML, **anclada al valor canónico** de `dat.bcg` (τ²=0.3132, I²=92.22 %, que este REML
> reproduce exacto) — así la comparación es contra el REML de metafor, no contra el propio motor de aNova.

### 3.10. SROC (precisión diagnóstica: sensibilidad/especificidad combinadas)

`calcularSROC` se validó con el dataset **AuditC** (14 estudios; Kriston, Hölzel, Weiser, Berner & Härter 2008),
incluido en el paquete de referencia **`mada`** de R para meta-análisis diagnóstico. Los conteos exactos (VP/FP/
FN/VN, con 2 estudios de celda cero) se extrajeron parseando el `.rda` oficial del paquete.

**Intervalo de Wilson** (sensibilidad/especificidad por estudio) — validado contra el ejemplo de texto clásico
(Newcombe 1998), x=81, n=263:

| Cantidad | aNova | Referencia |
|---|---|---|
| IC 95% inferior | 0.2553 | 0.2553 ✓ |
| IC 95% superior | 0.3662 | 0.3662 ✓ |

**Regresión de Moses-Littenberg** (D~S, WLS; Moses, Shapiro & Littenberg 1993) — validada contra una
implementación independiente en Python/numpy (ecuaciones normales) sobre AuditC:

| Cantidad | aNova | numpy |
|---|---|---|
| Intercepto (a) | 2.9159 | 2.9159 ✓ |
| Pendiente (b) | 0.3863 | 0.3863 ✓ |
| SE(a) / SE(b) | 0.0832 / 0.0846 | 0.0832 / 0.0846 ✓ |
| z(pendiente) | 4.565 | 4.565 ✓ |

**DOR combinado (DerSimonian-Laird, aleatorio)** — validado contra una implementación independiente en Python
que replica la misma fórmula (Glas et al. 2003; DerSimonian & Laird 1986):

| Cantidad | aNova | Referencia |
|---|---|---|
| DOR (IC 95%) | 27.944 (18.738–41.671) | 27.944 (18.738–41.671) ✓ |
| τ² | 0.339 | 0.339 ✓ |

> **Por qué no se compara literalmente contra `summary(madauni(AuditC))` de mada** (que publica DOR=26.337,
> τ²=0.311): se leyó el código fuente de `mada` (R/madauni.R) y aplica la corrección de continuidad (+0.5) a
> **las 4 celdas de todos los estudios siempre**, no solo a los que tienen una celda en cero. aNova sigue aquí,
> como en el resto del motor (`efectoDesde2x2`), la recomendación del Cochrane Handbook (6.4.3.1): corregir
> solo el estudio que de verdad lo necesita — la práctica moderna recomendada. Se confirmó reproduciendo el
> algoritmo de mada con su misma convención en Python: con "corregir siempre" sí se obtiene 26.337/0.311 exacto
> — la fórmula de *pooling* es idéntica a mada; la única diferencia real es el alcance de la corrección.
> Detalle completo en el encabezado de `scripts/benchmark-sroc.mts`.

### 3.11. Meta-análisis en red (Rücker, múltiples tratamientos)

`calcularNMA` (`src/red.ts`) implementa el método gráfico-teórico de Rücker (2012)
— el mismo algoritmo del paquete de referencia del campo, **`netmeta`** de R — incluido el ajuste de
varianza para estudios multibrazo (Rücker & Schwarzer 2014). Validado con el dataset **Stowe2010**
(29 estudios, enfermedad de Parkinson; 28 de 2 brazos + 1 de 3 brazos — LARGO), incluido en `netmeta`
(obtenido parseando su `.rda` oficial vía Python, sin R en este entorno). Referencia:
`netmeta(TE, seTE, treat1, treat2, studlab, common=FALSE, ref="Placebo")`, vignette de JSS 2023
(Balduzzi, Rücker, Nikolakopoulou et al.), sección 5.1.

| Cantidad | aNova | netmeta |
|---|---|---|
| k (estudios) / m (comparaciones) / n (tratamientos) | 29 / 31 / 4 | 29 / 31 / 4 ✓ |
| gl | 27 | 27 ✓ |
| Q | 48.15 | 48.15 ✓ |
| τ² / τ | 0.1518 / 0.3896 | 0.1518 / 0.3896 ✓ |
| I² | 43.9% | 43.9% ✓ |
| COMTI vs Placebo (DM, IC 95%) | −0.95 [−1.27, −0.63] | −0.95 [−1.27, −0.63] ✓ |
| Dopamine Agonist vs Placebo | −1.49 [−1.86, −1.12] | −1.49 [−1.86, −1.12] ✓ |
| MAOBI vs Placebo | −1.03 [−1.56, −0.50] | −1.03 [−1.56, −0.50] ✓ |

Coincidencia exacta en las 17 cantidades — incluida la parte más delicada del algoritmo, el ajuste de
varianza del estudio de 3 brazos (LARGO), cuya fórmula exacta (Lᵗ = −½·H·R·H con doble centrado,
pseudo-inversa, W = diag(diag(L))−L) se confirmó leyendo el código fuente real de `multiarm.R` en el
repositorio de `netmeta` en GitHub, no adivinada. El camino binario (OR) y el modelo de efectos fijos
se verificaron por separado con una prueba de consistencia interna (transitividad en una red sin ciclos:
ln(A/C) = ln(A/B)+ln(B/C) exacto a precisión de máquina).

**Ranking de tratamientos (P-scores)** — `calcularPScores`, fórmula exacta de Rücker & Schwarzer (2015),
confirmada leyendo `R/netrank.R` del código fuente de netmeta: P(i mejor que j) = Φ(−z) (o Φ(z) si un
valor mayor es mejor), P-score = promedio contra todos los demás tratamientos. Validado contra
`netrank(net1)` sobre el mismo Stowe2010:

| Tratamiento | aNova | netmeta |
|---|---|---|
| Dopamine Agonist | 0.9677 | 0.9677 ✓ |
| MAOBI | 0.5664 | 0.5664 ✓ |
| COMTI | 0.4659 | 0.4659 ✓ |
| Placebo | 0.0000 | 0.0000 ✓ |

Coincidencia exacta a 4 decimales en los 4 tratamientos.

**Inconsistencia local (node-splitting)** — `ResultadoNMA.nodeSplitting`, método de back-calculation de
König, Krahn & Binder (2013) (el método por defecto de `netsplit()` en netmeta), confirmado leyendo
`R/netmeta.R` (la fórmula exacta `TE.indirecto=(TE.red−P·TE.directo)/(1−P)`, `seTE.indirecto=
√(seTE.red²/(1−P))`, con P=proporción de información directa). La vignette no publica un ejemplo
numérico de `netsplit()`, así que se validó por **consistencia algebraica**: para COMTI vs MAOBI (con
evidencia directa real vía el estudio de 3 brazos LARGO), la fórmula garantiza que
`P·directo+(1−P)·indirecto` debe reproducir *exacto* el estimado de red ya validado arriba —
coincidencia exacta (4 decimales).

**Node-splitting por SIDDE** — el método alternativo de `netsplit()` (Separate Indirect from Direct
Design-based Evidence): en vez de despejar el indirecto algebraicamente, **descarta por completo** los
estudios con evidencia directa de un par (el estudio entero, no solo esa arista — confirmado leyendo
`R/sidde.R`) y **reajusta la red desde cero** con lo que queda; la estimación de red de ese reajuste es
el indirecto. Implementado reutilizando `calcularNMA` recursivamente (sin duplicar el álgebra de la red).
Tampoco tiene un ejemplo numérico publicado en la vignette, así que se validó con una implementación
Python/numpy independiente que quita el estudio LARGO (el único con comparación directa COMTI-MAOBI) y
reajusta la red completa desde cero:

| Cantidad | aNova | Python |
|---|---|---|
| COMTI vs MAOBI, indirecto (SIDDE) | 0.2507 | 0.2507 ✓ |

(Nota: SIDDE y back-calculation dan valores *distintos* para este par —0.2507 vs 0.1327— porque son
métodos genuinamente diferentes que solo coinciden exacto en redes sin redundancia, como se confirmó por
separado con una red triangular sintética de prueba donde ambos métodos SÍ dan resultados idénticos.)

**Descomposición de Q (heterogeneidad-dentro-de-diseño vs. inconsistencia-entre-diseños)** —
`calcularDecomposicionQ`, agrupa los estudios por "diseño" (conjunto exacto de tratamientos comparados;
{A,B} y {A,B,C} son diseños distintos aunque compartan el par A-B) y combina por mínimos cuadrados
generalizados (GLS) los que comparten diseño (fórmula confirmada leyendo `R/decomp.tau.R`; para un
diseño de 2 brazos se reduce exacto a la Q de Cochran ordinaria). Q-entre se obtiene por la identidad
algebraica de Krahn et al. (2013): Q-entre = Q_total − Q-dentro. Validado contra
`decomp.design(net1)` sobre Stowe2010:

| Cantidad | aNova | netmeta |
|---|---|---|
| Q-dentro (within designs) | 46.91 (25 gl, p=0.0050) | 46.91 (25 gl, p=0.0050) ✓ |
| Q-entre (between designs) | 1.24 (2 gl, p=0.5376) | 1.24 (2 gl, p=0.5376) ✓ |

Coincidencia exacta, incluido el desglose por diseño (15/14, 11/10, 2/1 estudios/gl para los 3 diseños de
2 brazos con ≥2 estudios; el diseño de 3 brazos LARGO, único en su tipo, no aporta —correctamente— al no
tener con qué compararse).

> **Alcance actual** (declarado también en la UI del Paso 7): el modelo NMA estándar, su descomposición
> de Q, el ranking por P-score, y la inconsistencia local por AMBOS métodos (back-calculation y SIDDE)
> están completos y validados. No incluye SUCRA (la versión por simulación del ranking — el P-score
> frecuentista da la misma conclusión sin necesitar simulación) ni gráficos especializados adicionales
> (net heat plot, embudo ajustado por comparación). Detalle completo en el encabezado de
> `scripts/benchmark-red.mts`.

### 3.12. Validación diferencial masiva contra R (octubre 2026)

Las secciones anteriores comparan el motor con **un dataset publicado** por método. Esta añade una prueba de otro
tipo: **705 conjuntos de datos sintéticos** (semilla fija; 534 + 31 de §3.16 + 60 de §3.17 + 80 de §3.18) se pasan por R (`metafor` 5.2.1, `meta` 8.5.0, `netmeta` 3.7.0) y
por el motor real de aNova, y se comparan **≈255.700 números** en 1.057 grupos de comparación. (Nota de cuentas: la primera versión
de este apartado, con 214 conjuntos, comparó ≈52.200 números; con meta-regresión, SROC y red fueron 374 conjuntos y ≈68.400 números
en 414 grupos. Antes de octubre de 2026 la documentación escribía «414 conjuntos», pero 414 eran los grupos de comparación, no los
conjuntos de datos; las medidas nuevas de §3.15 elevan las cifras a las de aquí.) Cubre tamaños de 2 a 100 estudios,
heterogeneidad nula y alta, eventos raros, celdas en cero y dobles-ceros, y redes con estudios de tres brazos.

| Familia | Conjuntos | Qué entra |
|---|:---:|---|
| Efecto genérico (efecto + SE) | 64 | modelos fijo/aleatorio, DL, REML, Paule-Mandel, Knapp-Hartung, intervalo de predicción, Q, I², IC de τ² e I² (Q-profile), Egger, Begg, Baujat, trim-and-fill, dejar-uno-fuera, acumulado, subgrupos |
| Tablas 2×2 | 70 | OR, RR y diferencia de riesgos por estudio; todo lo anterior sobre esos efectos; Mantel-Haenszel; Peto; Peters |
| Medias (DM, g de Hedges) | 40 | efecto y error estándar por estudio, y el pooling completo |
| Proporciones (logit) | 40 | ídem |
| **Tasas de incidencia** (RTI, DTI) | 40 | razón y diferencia de tasas con eventos y tiempo-persona (`escalc` IRR e IRD), con sus pruebas de sesgo, subgrupos, influencia, etc. |
| **Razón de medias** | 40 | ln(RM) y su EE de primer orden (fórmula del Handbook; ver la convención de §3.15) y el pooling completo |
| **Correlaciones** (z de Fisher) | 40 | z = atanh r con EE 1/√(n−3) (`escalc` ZCOR) y el pooling completo |
| **O−E y V** (HR) | 40 | ln HR = (O−E)/V con EE 1/√V, y que el efecto fijo sea exactamente Σ(O−E)/ΣV (Peto) |
| **Meta-regresión** (REML, pruebas z) | 60 | 1 a 3 moderadores (año, N, diseño categórico); coeficientes, EE, z, p, IC, τ² residual, QM y su p, R² |
| **SROC / precisión diagnóstica** | 60 | DOR, sensibilidad y especificidad combinadas (fijo, DL, PM, REML) con su IC, Q, I², τ², regresión de Moses-Littenberg (a, b, EE, z, p) e IC de Wilson por estudio |
| **SROC bivariado** (§3.16) | 31 | modelo binomial exacto y normal REML, curva, AUC y elipses contra implementación independiente, fuerza bruta, `lme4`, `mada` |
| **Efectos dependientes** (§3.17) | 60 | modelo CHE por REML y varianza robusta CR2 con Satterthwaite, contra `metafor::rma.mv` y `clubSandwich` |
| **Modelos mixtos logísticos** (§3.18) | 80 | proporciones y OR por máxima verosimilitud, contra una implementación independiente, fuerza bruta y `metafor::rma.glmm` |
| **Meta-análisis en red** | 40 | redes de 3 a 6 tratamientos y 8 a 25 estudios (≈12 % de tres brazos), diferencia de medias y OR, modelo común y aleatorio: efecto, EE e IC de cada par, Q, gl, I², τ², P-scores, descomposición de Q (dentro/entre diseños) y node-splitting por back-calculation y SIDDE |

**Resultado: 0 grupos fuera de tolerancia.** Efecto por estudio, modelo fijo, DL, REML, Paule-Mandel, Knapp-Hartung,
intervalo de predicción, Q, p de Q, I², IC de τ² e I², Baujat, Begg, dejar-uno-fuera, acumulado, subgrupos,
Mantel-Haenszel (OR, RR, RAR), Peto, Egger, Peters, trim-and-fill, meta-regresión, SROC y red coinciden con R entre
1e-8 y 1e-15 (la meta-regresión, el SROC y la red, a 1e-10 o mejor).

**Errores reales que esta prueba encontró en aNova y que ya están corregidos** (cada uno tiene su test de regresión en
`npm run test:regresion-r`; los cuatro primeros salieron de las familias originales, el 5 y el 6 de la ampliación a meta-regresión, SROC y red, y el 7 y el 8 de validar las medidas nuevas de §3.15, con su prueba en `npm run test:medidas-nuevas`):

1. **El p-valor de Egger y el de Peters salía DUPLICADO** (se multiplicaba por 2 una cola que ya era bilateral). Un
   sesgo de publicación "significativo" (p = 0,03) se mostraba como p = 0,06, y al revés. Corregido.
2. **El efecto calculado desde datos crudos se guardaba con 4 decimales**, y el motor reconstruye el error estándar a
   partir de los límites del IC; el resultado combinado se apartaba de R hasta un 1,6 %. Ahora se guarda con 9 cifras
   significativas (`formatoParaGuardar`).
3. **Mantel-Haenszel de la diferencia de riesgos excluía las tablas doble-cero** (R y el Cochrane Handbook las incluyen: son
   estudios con diferencia 0 y peso propio) y su Q no aplicaba la corrección de celdas cero. Corregido.
4. **Trim-and-fill con efectos aleatorios usaba el centro de efectos fijos** para decidir cuántos estudios imputar; ahora usa
   el del modelo elegido, como `metafor`. En un caso de la prueba pasaba de imputar 10 estudios a imputar 2.

5. **Node-splitting en redes con estudios de tres brazos.** aNova agrupaba la evidencia directa de un par con las varianzas
   ya ajustadas por multibrazo (las que pesan en la red completa); `netmeta` usa las varianzas originales de cada contraste
   (más τ² en efectos aleatorios) y define la proporción directa como Var(red)/Var(directa). Sin estudios multibrazo ambas
   coinciden (50 de 50 estimaciones exactas), pero con ellos fallaban 110 de 239: la estimación directa, la indirecta, la
   proporción directa, la z y el p. La red completa (efectos, Q, τ², P-scores, descomposición de Q) y la estimación
   indirecta de SIDDE no estaban afectadas. Corregido en `red.ts`.
6. **La cola normal era aproximada.** Los p-valores de z (coeficientes de la meta-regresión, pendiente de Moses-Littenberg,
   node-splitting y κ) usaban la aproximación de Abramowitz-Stegun (error absoluto hasta 1,5e-7; poco fiable para p muy
   pequeños). Ahora se calcula exacta (a ~1e-15) con la misma rutina de la χ² que ya coincidía con R.

7. **El REML se detenía antes de converger con efectos de magnitud pequeña.** Su criterio de parada era una tolerancia absoluta
   (1e-12); con efectos del orden de 0,001 (como una diferencia de tasas) τ² es del orden de 1e-8 y esa tolerancia deja un error
   relativo de ≈1e-4 en τ² y en lo que depende de él (intervalos, intervalo de predicción, Q entre subgrupos). La tolerancia es
   ahora relativa a la escala de las varianzas. El óptimo exacto de la verosimilitud restringida (maximizada con `optimize` sobre los
   datos re-escalados) confirma que aNova ahora es exacto; el REML de `metafor` tiene la misma limitación (su umbral también es
   absoluto), así que estas cantidades se comparan en una escala donde R es fiable (la diferencia de tasas «por 1000 personas-tiempo»).
8. **Begg con empates.** La τ de Kendall se calculaba sin corregir empates (τ-a) y su p exacta no es válida con empates; `metafor`
   usa τ-b y la varianza de Kendall con empates. Los empates son habituales (estudios con el mismo tamaño de muestra tienen la misma
   varianza) y, además, las varianzas se reconstruyen desde el IC y pueden diferir en el último dígito: ahora se detectan con una
   tolerancia de 1e-9. Corregido: τ-b, varianza con empates y aproximación normal cuando los hay.

Además se alineó **Peters** con su definición original y con `meta::metabias`: el +0,5 corrige solo el ln(OR) y el peso usa
los conteos originales (antes aNova ponderaba con los conteos ya corregidos).

**Diferencias de convención que se conservan, a propósito y documentadas** (no son errores; en cada una se comprobó que
aNova coincide con la fuente que dice seguir):

- **Mantel-Haenszel, diferencia de riesgos (RAR).** aNova sigue a **RevMan 5 y al Cochrane Handbook** (varianza de
  Greenland-Robins 1985, `SE = √J/K`; Q con las tablas doble-cero y +0,5 en celdas cero). `metafor::rma.mh` usa otro
  estimador de varianza (Sato-Greenland-Robins 1989) y deja las doble-cero fuera del Q. La prueba compara aNova con las
  fórmulas literales de RevMan 5 (transcritas en R, no con salida de un paquete) y mide la distancia a `metafor`: el error
  estándar se aparta hasta ≈2 % y el Q difiere en 14 de 70 conjuntos. La estimación puntual coincide con las dos.
  RevMan Web (la versión actual) no ofrece Mantel-Haenszel, así que la referencia de Cochrane sigue siendo RevMan 5.
- **Trim-and-fill con REML.** `metafor::trimfill()` reajusta el modelo con la tolerancia por defecto de REML (1e-5), por lo
  que su efecto ajustado solo es fiable a ≈1e-6. Se reajusta en R con control estricto sobre los mismos datos aumentados y
  aNova coincide con ese valor.
- **Referencias exactas de Paule-Mandel y del IC Q-profile.** La búsqueda de raíz de `metafor` se detiene antes de converger
  (diferencias en la 4.ª–5.ª cifra); la referencia se calcula en R con `uniroot(tol = 1e-15)` y es contra esa raíz exacta que
  se compara.

**Convenciones de las familias nuevas.**

- **SROC.** aNova corrige con +0,5 las cuatro celdas de un estudio si alguna es cero (el mismo estudio para sensibilidad,
  especificidad y DOR). La referencia aplica esa misma regla antes de llamar a R y usa `add = 0`: se compara el método, no la
  regla de corrección. Combina sensibilidad y especificidad por separado (univariado), no el modelo bivariado
  de Reitsma; esa limitación está declarada en la propia pantalla del Paso 7 y no se compara aquí.
- **Regresión de Moses-Littenberg.** El error estándar es el de efectos fijos (`rma(method = "FE", mods = ~S)`), sin reescalar
  por la varianza residual, que es lo que implementa aNova.
- **Red con OR.** La referencia solo usa tablas sin celdas en cero: aNova corrige cada brazo por separado y
  `netmeta::pairwise` corrige cada comparación, así que con ceros son reglas de corrección distintas.
- **Meta-regresión.** El R² solo se compara cuando el τ² sin moderadores no es un cero numérico en R (un cociente de dos ceros
  no significa nada); en ese caso se exige que aNova dé 0.

**Qué demuestra y qué no esta prueba.** Los datos son simulados: prueban que aNova calcula **lo mismo que R** sobre
datos de forma realista, incluidos los casos límite, no que los datos de un estudio real estén bien. No sustituye a los
benchmarks con datasets publicados de las secciones anteriores; los complementa. Quedan **fuera** de esta prueba masiva
las reglas de GRADE y el modelo bivariado (que aNova no implementa); los ensayos clúster y cruzados se comprobaron aparte contra el Handbook (§3.13).

### 3.13. Ensayos clúster-aleatorizados y cruzados: comprobación contra el Cochrane Handbook, cap. 23 (octubre 2026)

Para estos dos diseños **no existe un paquete de R que sirva de referencia** (son ajustes de las fórmulas, no un método
con software estándar), así que la comprobación es distinta de las anteriores: se leyeron las fórmulas del capítulo 23 del
*Cochrane Handbook for Systematic Reviews of Interventions* (versión 6.5, 2024; capítulo actualizado por última vez en octubre de 2019, página
<https://www.cochrane.org/authors/handbooks-and-manuals/handbook/current/chapter-23>) y se compararon una a una con el código. En
la página web las fórmulas son imágenes; se descargaron y se leyeron directamente, no de memoria ni de un resumen. Además se
reproducen **los ejemplos resueltos del propio capítulo** (23.1.4.1, 23.1.5 y 23.2.7.4) y las **equivalencias** que el
capítulo afirma. Todo está en `scripts/test-cluster-cruzado.mts` (`npm run test:cluster-cruzado`, 38 comprobaciones).

| Fórmula o regla | Handbook | Dónde está en aNova | Comprobación | Resultado |
|---|---|---|---|:---:|
| Efecto de diseño: DEFF = 1 + (M − 1) × ICC, con M = tamaño medio del clúster | 23.1.4 | `efectoDeDisenoCluster` (`statistics.ts`) | M = participantes totales / clústeres totales, como en el ejemplo del capítulo (295 + 330) / (10 + 11) | ✓ |
| Ejemplo: 10 + 11 aulas, 625 niños, ICC 0,02 → DEFF 1,576 | 23.1.4.1 | ídem | aNova da 1,5752: el libro redondea M a 29,8 antes de calcular (la aritmética del libro, comprobada) | ✓ |
| Tamaño efectivo: dividir participantes (y eventos, si es binario) entre el DEFF; con datos continuos solo N | 23.1.4 | no se usa (ver siguiente fila) | Ejemplo 63/295 y 84/330 → 40,0/187,2 y 53,3/209,4 | ✓ |
| Inflar el error estándar: EE × √DEFF, «equivalente a calcular tamaños efectivos» | 23.1.5 | `extraerEstudiosValidos` (`statistics.ts`) | Equivalencia demostrada numéricamente para OR, RR, diferencia de riesgos y diferencia de medias (celdas sin redondear: EE idéntico a 1e-12) | ✓ |
| Ejemplo: OR 0,795 (0,548 a 1,154), ln OR −0,23, EE 0,19 → 0,19 × √1,576 = 0,24 | 23.1.5 | `efectoDesde2x2` + ajuste | aNova: OR 0,795 (0,548 a 1,154), EE 0,190, EE inflado 0,239 | ✓ |
| Con datos de clúster incompletos no se ajusta «a ojo» | 23.1.1 (error de unidad de análisis) | el estudio se excluye del combinado | Probado: sin nº de clústeres o con ICC inválido el estudio queda fuera | ✓ |
| Diferencia de medias cruzada: MD = ME − MC | 23.2.7.1 | `efectoDesdeCruzado` | 7,0 − 6,5 = 0,5 | ✓ |
| EE(MD) = SDdiff / √N | 23.2.7.1 | ídem | Ejemplo: 2/√10 = 0,632 | ✓ |
| SDdiff = √(SDE² + SDC² − 2·Corr·SDE·SDC) | 23.2.7.1 | ídem | Ejemplo con Corr 0,68: SDdiff 1,8426, EE 0,583 | ✓ |
| Corr = (SDE² + SDC² − SDdiff²) / (2·SDE·SDC) | 23.2.7.3 | no se implementa (aNova no estima la correlación desde otro estudio) | Es la inversa de la anterior: devuelve 0,68 | ✓ |
| Con correlación 0 no hay ventaja del cruzado; con correlación positiva el análisis de grupos independientes da un IC demasiado ancho | 23.2.6, 23.2.7.1 | — | Con Corr = 0 el pareado coincide con el independiente; con 0,68 el independiente es 1,03 y el pareado 0,58 | ✓ |
| Desenlaces binarios cruzados requieren un estadístico | 23.2.5 | `calcularEfectoCrudos` solo calcula cruzado con DM | Con OR o DME devuelve «no calculable» | ✓ |

**Una errata del propio Handbook.** En el ejemplo 23.2.7.4 el libro escribe «SDdiff = 1.846» y en la línea siguiente usa 1,8426. El
valor correcto es 1,8426 (la raíz de 2,38² + 2,21² − 2 × 0,68 × 2,38 × 2,21); aNova calcula 1,8426.

**Lo que esta revisión corrigió en aNova** (las fórmulas estaban bien; lo que rodea a las fórmulas, no del todo):

1. **ICC mayor que 1.** El motor lo aceptaba; un coeficiente de correlación no pasa de 1. Ahora el estudio se excluye como cualquier dato de clúster inválido.
2. **Riesgo de ajustar dos veces.** El Handbook (23.1.3) dice que lo ideal es el efecto de un análisis que ya tenga en cuenta el clúster. El panel de clúster ahora lo dice: solo se marca esa variante si el efecto e IC capturados vienen de un análisis que ignora el clúster.
3. **Correlación imputada en cruzados.** El manuscrito decía «derivada con una correlación de X»; ahora dice «imputada de X, no reportada por el estudio». La pantalla recuerda que el Handbook (23.2.7.3) pide repetir el análisis con otros valores.
4. **Diferencia de medias estandarizada con clúster.** El Handbook (23.1.4) remite a consideraciones especiales (White y Thomas, 2005); aNova aplica el mismo ajuste del error estándar y ahora avisa en pantalla de que conviene consultar a un estadístico.

**Qué no cubre esta comprobación** (se deja dicho para no sugerir más de lo que hay):

- **Diferencia de medias estandarizada cruzada** (23.2.7.2) no está implementada; aNova solo calcula cruzados con DM.
- **Cruzados con desenlace binario, cruzados por clúster (23.2.9), ensayos escalonados (23.1.7) y ensayos individuales con clústeres (23.1.8)** no se calculan; el caso binario se bloquea con un aviso.
- **Tamaños de clúster desiguales:** aNova usa el tamaño medio (total/clústeres). El Handbook admite «otras formas más apropiadas» de estimar M cuando varían mucho (Eldridge 2006) y no se implementan.
- **Origen del ICC y de la correlación, y análisis de sensibilidad (resuelto el 5 de octubre de 2026).** El Handbook (23.1.6,
  23.2.7.3 y 23.2.8) pide declarar cómo se trataron estos datos y probar otros valores. aNova ahora guarda el origen del ICC
  (lo reporta el estudio, fuente externa o supuesto) y el de la correlación de los cruzados (la reporta el estudio, calculada
  con otro estudio de la revisión, fuente externa o supuesto), y el Paso 7 ofrece dos **sensibilidades**: al **ICC** (un mismo
  ICC de 0, 0,01, 0,02, 0,05, 0,10 o el que se quiera, aplicado a todos los ensayos por clúster) y a la **correlación** (0, 0,25,
  0,5, 0,75, 0,9 o la que se quiera, aplicada a todos los cruzados cuya DE de las diferencias se imputó). Ambas dicen si la
  significancia cambia y pueden llevarse al manuscrito como análisis exploratorios. La de la correlación calcula además la
  que implican los cruzados que sí reportan las tres DE (Corr = (SDE² + SDC² − SDdiff²) / (2·SDE·SDC), 23.2.7.3), y solo
  reajusta los cruzados cuyo efecto guardado coincide con sus datos crudos. Se comprueban contra cálculos manuales
  independientes en `npm run test:sensibilidad-icc` y `npm run test:sensibilidad-correlacion`.
- **Estudios con más de dos grupos y factoriales (23.3)** no se revisaron en esta pasada.

### 3.14. Funciones especiales de probabilidad contra R (octubre 2026)

De ellas salen todos los p de χ² (heterogeneidad, Q-profile, meta-regresión), los p y cuantiles de la t (Hartung-Knapp, intervalo de
predicción, Egger, Peters) y la cola de la normal. Se reescribieron en `src/funciones-especiales.ts` (fórmulas en su
cabecera; procedencia en `docs/PROCEDENCIA.md`) y se compararon con las funciones base de R (`pchisq`, `pt`, `qt`, `pnorm`, `lgamma`;
`validacion/generar-referencia-funciones.R`) con `npm run test:funciones-especiales`:

| Función | Puntos | Rango | Peor error relativo |
|---|---|---|---|
| χ², cola superior | 378 | gl de 0,5 a 1.000; p de ≈1 a ≈1e-290 | 3,4e-13 |
| t de Student, dos colas | 240 | gl de 0,5 a 10.000; t de 1e-6 a 1.000 | 3,9e-12 |
| Normal, dos colas | 18 | z de 0 a 30 | 3,3e-14 |
| ln Γ | 22 | x de 0,01 a 5.000 | 4,3e-14 |
| Cuantil 97,5 % de la t | 15 | gl de 1 a 1.000 | 8,9e-14 |

Además se comprueban identidades matemáticas que no dependen de R (P + Q = 1, simetría de la beta, χ² con 2 gl = e^(−q/2), t con 1 gl = Cauchy,
Γ(n) = (n−1)!). La prueba se verificó con mutaciones (cambiar un coeficiente de Stirling o un término de la fracción de la beta la hace fallar).
La versión anterior era exacta hasta p ≈ 1e-8 y devolvía 0 con p < 1e-15; ninguna conclusión cambia, pero ahora las colas extremas son exactas.

### 3.15. Medidas nuevas, influencia, conversiones y reglas de datos (octubre 2026)

Ampliación del alcance del motor (ver `docs/ALCANCE_MOTOR.md`). Todo lo siguiente salió de las fórmulas publicadas (Handbook,
artículos citados) y se comprobó contra R usando sus paquetes como **caja negra** (solo sus resultados; ver `docs/PROCEDENCIA.md`).

**Medidas nuevas** (4 familias × 40 conjuntos en `npm run validacion`, más `npm run test:medidas-nuevas`):

| Medida | Cálculo | Contra R |
|---|---|---|
| Razón de tasas (RTI) | ln[(e₁/t₁)/(e₀/t₀)], EE √(1/e₁ + 1/e₀); +0,5 a los dos conteos si alguno es cero | `escalc("IRR", add = 1/2, to = "only0")` |
| Diferencia de tasas (DTI) | e₁/t₁ − e₀/t₀, EE √(e₁/t₁² + e₀/t₀²) | `escalc("IRD")` |
| Razón de medias (RM) | ln(m₁/m₂), varianza de primer orden (Handbook 6.5.1.3) | fórmula; `escalc("ROM")` difiere por una corrección de sesgo (ver abajo) |
| Correlación (r) | z = atanh r, EE 1/√(n−3); se combina en z y vuelve a r con tanh | `escalc("ZCOR")` |
| HR desde O−E y V | ln HR = (O−E)/V, EE 1/√V; el modelo fijo da Σ(O−E)/ΣV (Peto) | identidad exacta |

Convenciones documentadas (no son errores): **(a)** `metafor::escalc("ROM")` suma al logaritmo una corrección de sesgo de primer orden
(Hedges y cols. 1999); aNova usa la fórmula del Handbook (Friedrich 2008) sin ella; la diferencia máxima en los 380 estudios de la
prueba es 1,5e-2 en el logaritmo (O(1/n)). **(b)** La diferencia de tasas se combina «por 1000 personas-tiempo» en la comparación contra R
por el motivo del error 7. **(c)** Un estudio sin eventos en ningún brazo se excluye de la razón de tasas, igual que del OR y el RR.

**Estudios sin eventos en ningún brazo (Handbook 10.4.4.2).** Las tablas 2×2 sin eventos en ningún brazo (o con el evento en todos)
no informan el OR ni el RR y se excluyen de la combinación (el `drop00` de `metafor`); en la diferencia de riesgos sí entran. Antes
aNova las dejaba dentro con +0,5 solo en el cálculo principal. La referencia de R se regeneró con la misma regla (28 de los 70 conjuntos
de tablas 2×2 tienen alguna tabla excluida); `npm run test:sin-eventos`.

**Diagnóstico de influencia y gráfico de Galbraith** (Viechtbauer y Cheung 2010; diez cantidades por estudio, más las dos coordenadas
del gráfico, con los modelos fijo, DL y REML en cada familia): residuo estandarizado externamente, DFFITS, distancia de Cook, razón de covarianzas, τ² y Q sin el estudio, hat, peso, DFBETAS y la marca de
influyente, y las coordenadas de `radial()`. Coinciden con `influence()` y `radial()` a ≈1e-10 en 835 estudios
por cantidad (OR) y en todas las demás familias, y la marca de «influyente» coincide en todos. Las definiciones se fijaron comparando
candidatas con la salida de `influence()` (por ejemplo, el error estándar de DFBETAS sale de todos los estudios con el τ² de la
exclusión). Los umbrales son los convencionales de `metafor` y se presentan como orientativos. `npm run test:influencia`.

**Conversiones de datos del Handbook** (cap. 6 y 10.6; `npm run test:conversiones`, 1.599 comprobaciones): DE desde el EE, desde el IC
de una media o de una diferencia, desde la t o desde el P; EE e IC desde el P y desde el IC; combinar grupos; DE del cambio y
correlación implícita; OR ↔ DME; media y DE desde mediana y cuartiles o rango (Wan y cols. 2014). Las identidades se prueban contra datos
simulados con la verdad conocida (R calcula con `t.test` el IC, la t y el p; aNova debe recuperar la DE verdadera, a ≈1e-13) y contra el
cálculo directo (combinar grupos, DE del cambio). Wan se compara con `meta::metamean` como caja negra: la media coincide siempre; la DE
a ≈6e-4 como máximo con n > 50 y a ≈1e-9 con n ≥ 400. Con n ≤ 50, `meta` usa las tablas del artículo, y para tamaños que no son de la forma 4Q+1
las aplica por tramos (un sondeo con n de 6 a 8 dio diferencias de hasta ≈9 % respecto a la fórmula del propio artículo); aNova usa las
fórmulas finales del artículo (Blom) para cualquier n y lo declara.

**Ensayos con más de dos grupos (Handbook 6.2.9, 23.3.4).** La división del grupo control compartido (eventos y total en binarios; solo
el N en continuos), la detección de estudios repetidos y que Mantel-Haenszel, Peto, Peters y el riesgo basal usen el mismo control
dividido tienen su prueba (`npm run test:ensayos-repetidos`, 38 comprobaciones). No existe un paquete de R que haga esta división,
así que se comprueba contra el cálculo manual del Handbook; es una aproximación (las comparaciones se tratan como independientes
aunque comparten el control) y Métodos lo declara.

**κ de Cohen y Fleiss (reescrito).** `npm run test:kappa-r`: 177 conjuntos contra `irr` 0.85 (caja negra); ver §3.3.

### 3.16. Modelo bivariado de precisión diagnóstica contra R (octubre 2026)

El modelo recomendado por el Handbook de precisión diagnóstica (cap. 10) es el **bivariado de efectos aleatorios** (Reitsma 2005; Chu y Cole
2006): el logit de la sensibilidad y el de la especificidad de cada estudio siguen una normal bivariada entre estudios (σ_sens, σ_esp, ρ).
`sroc-bivariado.ts` lo ajusta de dos maneras: **(a) principal, verosimilitud binomial exacta** dentro de cada estudio (máxima verosimilitud,
cuadratura de Gauss-Hermite adaptativa de 21×21 puntos con nodos de Golub-Welsch; no necesita corregir celdas en cero); **(b) aproximación
normal de Reitsma (REML)**, que es la de `mada::reitsma`, incluida como comparación y como respaldo declarado. Resultados: sensibilidad y
especificidad resumen con IC de Wald en logit, curva SROC de Rutter-Gatsonis dibujada solo en el rango de FPR observado, área bajo la curva y
regiones de confianza y de predicción del 95 %.

**Cómo se validó** (`validacion/generar-referencia-bivariado.R` → `referencia-bivariado.json` → `comparar-bivariado.mts`, dentro de
`npm run validacion`): 30 conjuntos simulados (k = 5 a 35; heterogeneidad, correlación y tamaños variados; estudios pequeños con celdas en
cero) más `mada::AuditC`. Cada pieza se compara con la fuente que le corresponde:

| Pieza | Contra | Resultado (dif. máx.) |
|---|---|---|
| Log-verosimilitud binomial en el óptimo | implementación **independiente** en R (otra cuadratura con `eigen`, `optim`, `numDeriv`) | 6e-6 (31 conjuntos) |
| μ, σ, ρ del modelo binomial | ídem | 2e-5, 3e-5 y 2e-3 (ρ en zonas planas) |
| Covarianza de las medias | ídem (Hessiano de `numDeriv`) | 3e-6 |
| Cuadratura 21×21 | **fuerza bruta** (`integrate` anidado por tramos), 5 puntos, 6 estudios | 6e-9 |
| Modelo binomial | `lme4::glmer` (Laplace), solo como comprobación de que el modelo es el mismo | μ 9e-3, σ 5e-2 (es otra aproximación) |
| Aproximación normal REML | óptimo **exacto** de la misma función objetivo en R | σ 4e-5, ρ 9e-4 |
| Aproximación normal REML | `mada::reitsma` (caja negra) | μ 5e-4, var(μ) 1e-5, sens. e IC 4e-5 |
| Curva SROC | `mada::sroc` | 1e-14 (248 valores) |
| AUC | `mada::AUC` (rejilla fina) | 5e-5 |
| Elipses de confianza | forma cuadrática d′V⁻¹d = χ²₂(0,95) para los puntos de aNova y los de `ROCellipse` | 1e-12; extensión 8e-5 |

**Resultado: 0 fuera de tolerancia** (≈19.000 comparaciones nuevas). Hallazgos de la validación, ya corregidos: con ρ o una varianza en el límite,
el error estándar de las medias se calculaba condicionando a todos los parámetros y salía hasta un 22 % más bajo que el de R; ahora se
invierte el Hessiano de los parámetros libres con el del límite fijado (coincide con R a 3e-6). Una referencia de fuerza bruta inicial con
`integrate` en un solo intervalo subestimaba la verosimilitud con muchos pacientes (pierde el pico estrecho): se cambió a tramos y a un conjunto
pequeño; la cuadratura, estable de 21 a 81 nodos, era la correcta.

**Dónde mada difiere de aNova (documentado):** `mvmeta`, que usa `mada`, se detiene antes de converger; con una varianza casi nula (< 0,05) su σ
puede quedar a ~1e-2 del óptimo exacto, que aNova sí alcanza (se comprueba contra el óptimo exacto). Por eso la comparación directa con `mada`
usa una tolerancia de 2e-2 en σ solo en esa zona.

**Límites declarados al usuario:** con pocos estudios (< 10) σ y ρ se estiman con poca precisión; si ρ = ±1 o σ ≈ 0 el ajuste lo dice; las
estimaciones son un promedio sobre los umbrales (efecto de umbral: Handbook 10.4.1); el modelo no incluye covariables ni compara pruebas.
Pruebas sin R: `npm run test:sroc-bivariado` (cuadratura, simetrías, casos límite, geometría, caché; 7 mutaciones atrapadas).

### 3.17. Efectos dependientes: modelo multinivel y varianza robusta contra R (octubre 2026)

Cuando un estudio aporta varias filas de la misma medida (varios desenlaces, seguimientos o comparaciones), tratarlas como independientes
sobreestima la precisión. `multinivel.ts` ajusta el modelo **CHE** (efectos correlacionados y jerárquicos; Pustejovsky y Tipton 2022):
yᵢⱼ = x′β + uᵢ + wᵢⱼ + eᵢⱼ, con σ²ᵤ entre estudios y σ²ᵥ entre efectos del mismo estudio (REML), error de muestreo conocido y una correlación ρ
supuesta entre los de un mismo estudio. La inferencia usa el **error estándar robusto CR2** (Bell y McCaffrey, con la corrección de
Pustejovsky y Tipton 2018) y los **grados de libertad de Satterthwaite** para un coeficiente (Tipton 2015); con menos de 4 df la app avisa (la
aproximación se considera válida con df > 4: simulaciones de Tipton 2015).

| Pieza | Contra | Resultado (dif. máx.) |
|---|---|---|
| Estimador sándwich CR2, df de Satterthwaite y coeficientes **aislados** (con los σ² de metafor) | `clubSandwich::coef_test(vcov = "CR2")` (caja negra) | 1e-15, 1e-13 y 6e-15 (99 coeficientes en 60 conjuntos) |
| σ² entre estudios y entre efectos | `metafor::rma.mv` (REML, `random = ~ 1 \| estudio/efecto`, V de `vcalc`) | 4e-6 (metafor se detiene antes de converger) |
| Log-verosimilitud restringida (con constantes) | `rma.mv` | 1,4e-6 (aNova queda por encima) |
| Coeficientes y EE del modelo | `rma.mv` | 2e-6 y 4e-6 |
| EE robusto, df, valor p e IC del ajuste completo | `clubSandwich` | 5e-7, 4e-6, 7e-7 y 2e-6 |

60 conjuntos simulados con estudios de 1 a 12 efectos, ρ de 0 a 0,8, σ² entre estudios y entre efectos desde 0, y moderadores a nivel de efecto o de
estudio. **Resultado: 0 fuera de tolerancia.** Pruebas sin R (`npm run test:multinivel`): solución analítica del CR2 con observaciones
independientes (HC2 ponderada), df exactos = m − 1 en un diseño balanceado, equivalencia con el REML de efectos aleatorios cuando cada estudio
aporta una fila, invariancias (orden, nombres, escala, desplazar un moderador) y 5 mutaciones atrapadas.
**Dónde difiere `metafor`:** su optimizador de σ² se detiene a ~4e-6 del óptimo; por eso el ajuste completo se compara con esa tolerancia y el sándwich,
aislado, a precisión de máquina. **Límites declarados:** la ρ es un supuesto (la app muestra la sensibilidad a ρ = 0 … 0,8); la inferencia robusta
vale para un coeficiente a la vez (no hay contrastes de varios coeficientes: Tipton y Pustejovsky 2015); sin covariables en la interfaz (el motor las admite y
están validadas); no hay Hartung-Knapp aquí.

### 3.18. Modelos mixtos logísticos (GLMM) contra R (octubre 2026)

`glmm.ts` modela los recuentos directamente: **proporciones** con el modelo binomial-normal (logit(pᵢ) = μ + uᵢ) y **odds ratio** con el modelo
hipergeométrico-normal (verosimilitud condicional exacta; «CM.EL» de `metafor::rma.glmm`; Stijnen y cols. 2010), por máxima verosimilitud con
cuadratura de Gauss-Hermite adaptativa de 21 puntos. No usa corrección de continuidad, conserva los estudios con un brazo en cero y declara, sin
corregirlos, los que no aportan información (sin eventos en ningún brazo o con evento en todos). El error estándar sale de la inversa del Hessiano
completo (μ, ln τ).

| Pieza | Contra | Resultado (dif. máx.) |
|---|---|---|
| Log-verosimilitud, μ, τ y EE de μ (40 proporciones y 40 tablas 2×2) | implementación **independiente** en R (Gauss-Hermite adaptativa propia, `optim`, `numDeriv`) | 5e-5, 1,4e-4, 6e-7 y 3e-5 (rel. 4,6e-3 en proporciones con τ grande) |
| Verosimilitud marginal (5 puntos de cada familia) | **fuerza bruta** (`integrate` por tramos); probabilidades de la hipergeométrica no central de `BiasedUrn::dFNCHypergeo` | 1,2e-7 |
| Modelo | `metafor::rma.glmm` (PLO y CM.EL; solo donde converge: 23 de 40 tablas) | μ 2e-6 (OR) y 5e-3 (prop.); EE hasta 1,1 % |

**Resultado: 0 fuera de tolerancia.** Dónde difiere: (1) el EE de `metafor` para proporciones (vía `lme4`) está condicionado a τ; el de aNova usa el
Hessiano completo, y la diferencia es de hasta el 1,1 % en estos conjuntos (con τ grande llegaría al 15 %: está fijado con dos casos en
`test:glmm`); (2) `rma.glmm(CM.EL)` no converge en 17 de 40 tablas, así que la referencia principal es la implementación independiente; (3) con
τ ≈ 0 la verosimilitud es casi plana y el BFGS de R se detiene en τ ≈ 0,003–0,007 con una verosimilitud algo menor que la de aNova, que llega al
límite (se acepta si la de aNova no es peor). **Hallazgo:** la cuadratura adaptativa se corrompe con más de ~100 nodos (con integrando estrecho, los
nodos de las colas multiplican un cociente integrando/normal que crece); se comprobó que 21 a 61 coinciden a 1e-14 y `nodosGaussHermite` ahora
rechaza más de 64. **Límites declarados:** solo OR y proporciones (no RR ni tasas); intervalos de Wald sin ajuste para pocos estudios (aviso con
< 5); sin intervalo de predicción; los ensayos por clúster o cruzados y los de control dividido (recuentos no enteros) se excluyen y se declaran;
τ en el tope de 30 se declara como heterogeneidad no estimable. Pruebas sin R: `npm run test:glmm` (hipergeométrica y binomial calculadas a mano,
simetrías, ceros, dobles-ceros, integración con el manuscrito y 8 mutaciones atrapadas).

### 3.19. Validación masiva de la red y del SROC univariado contra `netmeta`, `meta`, `metafor` y `mada` (octubre 2026)

Tras reescribir `red.ts` y la regresión de Moses-Littenberg desde los artículos (docs/PROCEDENCIA.md), se amplió la validación con **una referencia masiva**
de R (`AVANZADA_MASIVA=1 Rscript validacion/generar-referencia-avanzada.R`, semilla 20261010; `npm run validacion:masiva`; la generación tarda ~50 min
por el node-splitting SIDDE de `netmeta`): **300 redes** de 3 a 7 tratamientos y 6 a 30 estudios (DM y OR; estudios de **2, 3 y 4 brazos**: 4.211, 612 y 221)
y **400 conjuntos de SROC** de 3 a 30 estudios, comparados con `netmeta` 3.7.0, `meta`, `metafor` y `mada` 0.5.12 (`madauni`, DerSimonian-Laird, corrección de
continuidad «single»).

**Resultado: 328.574 comparaciones en 1.050 grupos, 0 fuera de tolerancia (346.745 tras alinear las celdas en cero: ver el hallazgo 1).** En las 225 redes sin celdas en cero: efectos de red, errores estándar e IC de todos
los pares, Q, gl, I², τ², descomposición de Q, P-scores y node-splitting por back-calculation y por SIDDE coinciden con `netmeta` entre 1e-13 y 1e-15 (la
descomposición de Q, a 5e-13). En los 400 SROC: DOR, IC, τ² y gl coinciden con `mada::madauni` a 6e-12 o mejor; la regresión de Moses-Littenberg (328 conjuntos con ≥ 4
estudios) con `metafor` a 2e-13; el intervalo de Wilson por estudio (8.638 intervalos) con `prop.test` a 8e-11.

**Hallazgo 1 — celdas en cero en la red (convención distinta: encontrada, cuantificada y ALINEADA).** Las otras 75 redes de OR tenían brazos con 0 eventos en 37 de ellas, y
las 37 diferían de `netmeta` (efecto de red, log OR: mediana 5,5e-3, máxima 5,5e-2; error estándar hasta ~3 %). Causa: aNova sumaba 0,5 a las celdas del BRAZO con el cero
y conservaba los contrastes de doble cero; `netmeta` (vía `meta::pairwise`) aplica la corrección **por estudio** (si alguna celda de algún brazo es cero, +0,5 a las
cuatro celdas de TODAS las comparaciones del estudio, incluso las que no tienen ceros: se comprobó con la salida de `pairwise`) y descarta los contrastes de doble cero (el
Handbook 10.4.4.2 excluye de un OR los estudios sin eventos en ambos brazos). No lo había introducido la reescritura (la versión anterior daba lo mismo). **El 8 de octubre de 2026 `red.ts` se alineó**:
los contrastes se calculan por pares de brazos con la corrección por estudio, los estudios de dos brazos con doble cero (sin eventos o con evento en todos) se
descartan y se declaran, y la reducción de pesos de los estudios multibrazo parte de las varianzas de contraste (Σ_ab = (R_0a + R_0b − R_ab)/2; con varianzas aditivas es la de antes, y
coincide con el `seTE.adj` de `netmeta`). Un estudio multibrazo con un contraste de doble cero se conserva con su corrección (equivale a `pairwise(allstudies = TRUE)`: `netmeta` no admite
ese caso por defecto); uno sin ningún contraste informativo se descarta y se declara. **Resultado: las 75 redes con ceros coinciden con `netmeta` a ~2e-14 en efectos de red y a 1e-15 en τ²
(346.745 comparaciones, 0 fuera de tolerancia con la referencia masiva) y las 225 sin ceros siguen igual.** Seis casos de la red base (3 brazos con un cero, doble cero entre brazos,
dobles «todos con evento», esquinas opuestas) se fijaron con los valores de `netmeta` en `test:red`.

**Segunda tanda de ceros (120 redes de OR, `referencia-red-ceros.json`; `AVANZADA_MASIVA=ceros`).** 60 redes `redm` (estudios multibrazo con contrastes de doble cero, contra `pairwise(allstudies = TRUE)`) y 60 `redh` (tasas altas: brazos con eventos = n y dobles «todos con evento»). Total con la masiva: **374.147 comparaciones en 1.050 grupos, 0 fuera de tolerancia.** Las `redh` y 51 de las 60 `redm` coincidieron a la primera. Las otras 9 `redm` (exactamente las que tenían un estudio multibrazo con todos los brazos en cero) mostraron una **diferencia de regla**: `allstudies = TRUE` conserva ese estudio con TE = 0 y varianza enorme (k, Q, gl y τ² cambian), mientras que aNova lo descarta y lo declara por no informar ningún contraste (mismo criterio que el doble cero de dos brazos, Handbook 10.4.4.2). Se mantiene la regla de aNova; a R se le entregan los mismos datos sin esos estudios y a aNova los completos (así se comprueba que los descarta). Con eso las 9 coinciden en todo.

**Hallazgo 2 — la «Q de Cochran» de `madauni`.** Con el método DSL, `madauni` devuelve la Q calculada con los pesos de EFECTOS ALEATORIOS (1/(v + τ²)), no la de Cochran con pesos de
efectos fijos que da el I². La de aNova es la de efectos fijos y coincide con la de `metafor` y con el cálculo a mano (p. ej. AuditC recortado: 33,319 frente a 15,651 de `madauni`);
el τ² y el DOR coinciden. Por eso esa Q no se compara con `mada` (sí con `metafor` y `meta`, §3.12).

## 4. Alcance honesto: qué prueba esto y qué no

Ser transparente sobre los límites es parte de ser confiable.

**Lo que esta validación demuestra (el núcleo que usa la gran mayoría de las revisiones publicadas):**

- Combinación de efectos por inverso-varianza (modelos fijo y aleatorio).
- Estimador **REML** de la heterogeneidad (también hay DerSimonian-Laird y Paule-Mandel disponibles).
- Ajuste **Knapp-Hartung** del intervalo de confianza.
- Heterogeneidad: **Q, I², τ², τ** y su valor de p.
- Intervalo de predicción.
- **Mantel-Haenszel** para datos binarios (OR, RR y diferencia de riesgos), con la varianza de
  Robins-Breslow-Greenland — validado contra `metafor::rma.mh` (§3.2).
- **Concordancia entre revisores**: κ de Cohen y de Fleiss — validado contra `irr` (§3.3).
- **Método de Peto** (OR de un paso, eventos raros) — validado contra `metafor::rma.peto` (§3.4).
- **g de Hedges** (diferencia de medias estandarizada) — validado contra `metafor::escalc(measure="SMD")` (§3.5).
- **Meta-análisis acumulativo** (efecto + IC por paso) — validado contra `metafor::cumul` (§3.6).
- **Test de Begg** (correlación de rangos, τ de Kendall + p exacta) — validado contra `metafor::ranktest` (§3.7).
- **IC de τ²/τ/I² por Q-profile** — IC de τ²/τ validados contra `meta` (§3.8).
- **Meta-regresión** (moderadores múltiples y categóricos, REML) — validada contra `metafor::rma(mods=…)` (§3.9).
- **SROC** (precisión diagnóstica: sensibilidad/especificidad combinadas, DOR, regresión de Moses-Littenberg,
  intervalo de Wilson) — validado contra el dataset AuditC de `mada` e implementaciones independientes en
  Python/numpy (§3.10). Nota de alcance: la sensibilidad/especificidad se combinan por separado (univariado),
  no con el modelo bivariado completo (Reitsma/HSROC) — ver el detalle en §3.10 y en la propia UI del Paso 7.
- **Meta-análisis en red** (múltiples tratamientos, método de Rücker con ajuste multibrazo) — validado
  contra el dataset Stowe2010 de `netmeta` (§3.11): ranking por P-score (`netrank`), descomposición de Q
  dentro/entre diseños (`decomp.design`), e inconsistencia local por node-splitting (back-calculation y
  SIDDE) — ver §3.11.

- **Validación diferencial masiva** contra R sobre 705 conjuntos sintéticos y ≈255.700 comparaciones, que incluye
  subgrupos, dejar-uno-fuera, trim-and-fill, Egger, Peters, las medidas de efecto por estudio, la meta-regresión, el SROC
  y el meta-análisis en red (§3.12), y las medidas nuevas, la influencia y Galbraith (§3.15). Encontró y corrigió ocho errores reales.

**Lo que aún NO cubre este benchmark** (son diferencias de *amplitud*, no de *exactitud*):

- Los ensayos clúster y cruzados no tienen paquete de R de referencia: se comprobaron contra las fórmulas y los
  ejemplos del Handbook cap. 23 (§3.13, con sus límites). GRADE son reglas de decisión, no un cálculo comparable.

En otras palabras: en lo que aNova hace, lo hace con exactitud de gold standard; lo que le falta son
métodos adicionales, no precisión.

## 5. Cómo reproducir la prueba tú mismo

La validación es una **prueba de regresión permanente** dentro del proyecto. Cualquiera puede correrla:

```bash
npm run benchmark        # corre las cuatro comprobaciones
# o por separado:
npm run benchmark:meta    # meta-análisis (ThirdWave + dat.bcg) vs meta/metafor
npm run benchmark:kappa   # concordancia (diagnoses) vs irr
npm run benchmark:sroc    # SROC (AuditC) vs mada / numpy
npm run benchmark:red     # meta-análisis en red (Stowe2010) vs netmeta

npm run validacion:masiva     # además, 400 SROC y 300 redes contra netmeta, metafor, meta y mada (§3.19); NO necesita R
npm run validacion            # 705 conjuntos sintéticos vs R (≈255.700 comparaciones); NO necesita R
npm run test:regresion-r      # los errores corregidos, uno por uno
npm run test:cluster-cruzado  # clúster y cruzados contra los ejemplos del Handbook, cap. 23 (§3.13)
npm run test:sensibilidad-icc # sensibilidad al ICC contra un cálculo manual independiente
npm run test:sensibilidad-correlacion # sensibilidad a la correlación de cruzados, ídem
npm run validacion:regenerar  # vuelve a generar la referencia con R y repite la comparación (SÍ necesita R)
```

La referencia de R está guardada en `validacion/referencias/referencia.json` y `referencia-avanzada.json`, así que `npm run validacion`
corre en cualquier equipo con Node. Regenerarla exige R con los paquetes `metafor`, `meta`, `mada`, `netmeta` y `jsonlite`
(en macOS: `brew install r cmake pkgconf`, y luego instalar los paquetes desde R). `validacion/generar-referencia.R`
y `generar-referencia-avanzada.R` fijan la semilla y registran las versiones de R y de los paquetes dentro del propio JSON.

Los scripts (`scripts/benchmark-meta.mts`, `scripts/benchmark-kappa.mts`, `scripts/benchmark-sroc.mts` y
`scripts/benchmark-red.mts`) alimentan datasets canónicos al
motor real de aNova, comparan cada cantidad con la cifra publicada de R, imprimen las tablas anteriores y
**fallan (código de salida 1) si alguna cifra deja de coincidir**. Así, si en el futuro alguien modifica el
motor estadístico y rompe la equivalencia sin querer, la prueba lo detecta de inmediato.

> Detalle técnico: aNova combina los efectos por inverso-varianza sobre el efecto transformado, y esa
> matemática es idéntica para cualquier medida (RR, OR, DM, etc.). Por eso el benchmark alimenta los
> efectos ya calculados (TE) y sus errores estándar (seTE) del dataset, y compara el *pooling* — que es
> exactamente lo que se está validando.

## 6. Un documento vivo

Esta validación **crece con la herramienta**. La política de aNova es que **cada nueva metodología
estadística** que se incorpore se valide con el mismo procedimiento cabeza-a-cabeza contra `meta`/`metafor`
sobre un dataset canónico publicado, se **añada a este documento** (su fila o sección de resultados) y se
**incluya en la prueba de regresión automática** (`npm run benchmark`). Así, la lista de métodos
verificados solo puede crecer, y ninguno entra a producción sin su comprobación numérica.

## 7. Fuentes

- Harrer M, Cuijpers P, Furukawa TA, Ebert DD. *Doing Meta-Analysis in R: A Hands-On Guide.* (dataset
  ThirdWave y salida de `metagen`): <https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/>
- Balduzzi S, Rücker G, Schwarzer G. *How to perform a meta-analysis with R: a practical tutorial.*
  Evid Based Ment Health. 2019 (paquete `meta`).
- Verificación complementaria del motor contra las fórmulas oficiales de Cochrane/RevMan (Deeks &
  Higgins) — ver el documento interno de estado de la aplicación.
- Higgins JP, Eldridge S, Li T. *Chapter 23: Including variants on randomized trials* (last updated October 2019). En: Higgins JP,
  Thomas J, Chandler J, et al. (eds.) *Cochrane Handbook for Systematic Reviews of Interventions*, versión 6.5. Cochrane, 2024: <https://www.cochrane.org/authors/handbooks-and-manuals/handbook/current/chapter-23>
  (fórmulas de §3.13; consultado el 5 de octubre de 2026).
- Doebler P, Holling H, Sousa-Pinto B. *Meta-Analysis of Diagnostic Accuracy with mada* (vignette del
  paquete `mada` de R; dataset AuditC): <https://cran.r-project.org/web/packages/mada/vignettes/mada.pdf>
- Moses LE, Shapiro D, Littenberg B. *Combining independent studies of a diagnostic test into a summary
  ROC curve.* Stat Med. 1993 (regresión de Moses-Littenberg).
- Newcombe RG. *Two-sided confidence intervals for the single proportion.* Stat Med. 1998 (intervalo de
  Wilson, ejemplo de referencia x=81/n=263).
- Rücker G. *Network meta-analysis, electrical networks and graph theory.* Res Synth Methods. 2012
  (método gráfico-teórico, base de `calcularNMA`).
- Rücker G, Schwarzer G. *Reduce dimension or reduce weight? Comparing two approaches to multi-arm
  studies in network meta-analysis.* Stat Med. 2014 (ajuste de varianza multibrazo).
- Balduzzi S, Rücker G, Nikolakopoulou A, Papakonstantinou T, Salanti G, Efthimiou O, Schwarzer G.
  *netmeta: An R Package for Network Meta-Analysis Using Frequentist Methods.* J Stat Softw. 2023
  (vignette; dataset Stowe2010): <https://cran.r-project.org/web/packages/netmeta/vignettes/netmeta.pdf>
- Stowe RL, Ives NJ, Clarke C, et al. *Meta-analysis of the comparative efficacy and safety of adjuvant
  treatment to levodopa in later Parkinson's disease.* Cochrane review (dataset Stowe2010, vía `netmeta`).

---

*Última actualización: 5 de octubre de 2026 (validación diferencial masiva, §3.12, y comprobación de clúster y cruzados contra el Handbook, §3.13).*
