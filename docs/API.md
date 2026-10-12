# Referencia de la API

Generada automáticamente por `scripts/exportar-motor.mts` desde los comentarios del código: 107 funciones exportadas, con su firma y la primera
descripción de su comentario. Los tipos de datos (`Estudio`, `ResultadoMetaAnalisis`, …) están en `src/tipos.ts` y en cada módulo. Las funciones
que reciben `Estudio` esperan los números como texto, tal como salen de un formulario (ver el ejemplo del README).

## `src/funciones-especiales.ts`

### `lnGamma`

```ts
export function lnGamma(x: number): number
```

ln Γ(x) para x > 0 (y para x no entero negativo, por reflexión).

### `gammaIncompletaQ`

```ts
export function gammaIncompletaQ(a: number, x: number): number
```

Q(a, x) = Γ(a, x) / Γ(a): cola superior de la gamma. Para a > 0 y x ≥ 0. Exacta también en colas extremas.

### `gammaIncompletaP`

```ts
export function gammaIncompletaP(a: number, x: number): number
```

P(a, x) = γ(a, x) / Γ(a): cola inferior de la gamma.

### `betaIncompleta`

```ts
export function betaIncompleta(x: number, a: number, b: number): number
```

I_x(a, b): beta incompleta regularizada, para 0 ≤ x ≤ 1, a > 0, b > 0.

### `chiCuadradoCola`

```ts
export function chiCuadradoCola(q: number, gl: number): number
```

p de la cola superior de χ² con `gl` grados de libertad: P(χ² ≥ q).

### `tStudentDosColas`

```ts
export function tStudentDosColas(t: number, gl: number): number
```

P(|T| > t) para una t de Student con `gl` grados de libertad: la cola de AMBOS lados (ya es el p bilateral).

### `normalDosColas`

```ts
export function normalDosColas(z: number): number
```

Cola de dos lados de la normal estándar: P(|Z| > z) = Q(½, z²/2).

### `tStudentCuantilDosColas`

```ts
export function tStudentCuantilDosColas(alfa: number, gl: number): number
```

Valor t > 0 tal que P(|T| > t) = `alfa` para una t de Student con `gl` grados de libertad (cuantil bilateral). Por bisección sobre `tStudentDosColas`, que baja de forma estricta al subir t. `alfa` en (0, 1).

### `normalCuantilDosColas`

```ts
export function normalCuantilDosColas(alfa: number): number
```

Valor z > 0 tal que P(|Z| > z) = `alfa` (cuantil bilateral de la normal estándar). `alfa` en (0, 1).

### `normalCuantil`

```ts
export function normalCuantil(p: number): number
```

Cuantil de la normal estándar Φ⁻¹(p), para p en (0, 1).

## `src/numerica-compartida.ts`

### `jacobiSimetrica`

```ts
export function jacobiSimetrica(a: number[][]):
```

Autovalores y autovectores de una matriz simétrica por el método de Jacobi cíclico. Devuelve los autovalores y los vectores en columnas.

### `nodosGaussHermite`

```ts
export function nodosGaussHermite(q: number): NodosGH
```

Nodos y pesos de Gauss-Hermite para la normal estándar con `q` puntos (Golub-Welsch: valores propios de la matriz de Jacobi de los polinomios de Hermite). Integran exactamente los polinomios hasta grado 2q − 1 contra la densidad normal estándar.

### `nelderMead`

```ts
export function nelderMead(f: (x: number[]) => number, x0: number[], opc: { paso?: number; tolX?: number; tolF?: number; maxEval?: number; reinicios?: number } = {}):
```

Minimiza `f` por Nelder-Mead con reinicios; devuelve el mejor punto y su valor.

### `hessianoNumerico`

```ts
export function hessianoNumerico(f: (x: number[]) => number, x: number[], h = 1e-3): number[][]
```

Hessiano numérico (diferencias centrales) de `f` en `x`, con paso relativo `h`.

### `invertir`

```ts
export function invertir(a: number[][]): number[][] | null
```

Inversa de una matriz cuadrada por Gauss-Jordan con pivote parcial; `null` si es singular.

### `cholesky`

```ts
export function cholesky(A: Mat): Mat | null
```

Factor de Cholesky inferior L (A = L L'); `null` si A no es definida positiva.

### `inversaSPD`

```ts
export function inversaSPD(A: Mat):
```

Inversa y log-determinante de una matriz simétrica definida positiva por Cholesky.

### `raizInversaSimetrica`

```ts
export function raizInversaSimetrica(B: Mat): Mat
```

Raíz cuadrada inversa simétrica (con inversa generalizada: los valores propios casi nulos se tratan como 0).

## `src/statistics.ts`

### `chiSquarePValue`

```ts
export function chiSquarePValue(q: number, gl: number): number
```

p de la cola superior de χ² con `gl` grados de libertad.

### `studentTQuantile975`

```ts
export function studentTQuantile975(gl: number): number
```

Percentil 97,5 de la t de Student con `gl` grados de libertad (cuantil bilateral del 95 %), por bisección.

### `transformar`

```ts
export function transformar(valor: number, escala: EscalaAnalisis): number
```

Lleva un valor a la escala de análisis: log (OR, RR, HR), logit (proporciones), z de Fisher (correlaciones) o el mismo valor (escala lineal).

### `destransformar`

```ts
export function destransformar(valor: number, escala: EscalaAnalisis): number
```

Inversa de `transformar`: devuelve un valor de la escala de análisis a su escala natural.

### `efectoDesde2x2`

```ts
export function efectoDesde2x2(eventosTrat: number, totalTrat: number, eventosCtrl: number, totalCtrl: number, medida: "OR" | "RR" | "RAR" | "PR"): EfectoCalculado | null
```

Calcula la medida de efecto (OR, RR, RAR) a partir de una tabla 2×2.

### `efectoDesdeProporcion`

```ts
export function efectoDesdeProporcion(eventos: number, total: number): EfectoCalculado | null
```

Calcula una proporción de un solo brazo (p.ej. prevalencia) a partir de eventos/total — sin grupo de comparación. Transforma en escala logit (igual que `escalaPorMedida.Proporción`) para el error estándar, y aplica corrección de continuidad cuando el resultado es 0% o 100% (donde logit(0) y logit(1) no están definidos — sin esto, esos estudios, perfectamente legítimos en la práctica, quedarían excluidos del meta-análisis o forzarían a calcular el IC a mano).

### `efectoDesdeMedias`

```ts
export function efectoDesdeMedias(mediaTrat: number, deTrat: number, nTrat: number, mediaCtrl: number, deCtrl: number, nCtrl: number): EfectoCalculado | null
```

Calcula la diferencia de medias (DM) a partir de medias, DE y n por grupo.

### `efectoDesdeMediasEstandarizada`

```ts
export function efectoDesdeMediasEstandarizada(mediaTrat: number, deTrat: number, nTrat: number, mediaCtrl: number, deCtrl: number, nCtrl: number): EfectoCalculado | null
```

Calcula la g de Hedges (diferencia de medias estandarizada, DME) a partir de medias, DE y n.

### `efectoDesdeTasas`

```ts
export function efectoDesdeTasas(eventosTrat: number, tiempoTrat: number, eventosCtrl: number, tiempoCtrl: number, medida: "RTI" | "DTI"): EfectoCalculado | null
```

Tasas de incidencia (eventos y tiempo-persona de cada brazo) → **razón de tasas** (`RTI`) o **diferencia de tasas** (`DTI`) (Cochrane Handbook 6.7, 10.8). Réplica de `metafor::escalc(measure = "IRR" / "IRD", add = 1/2, to = "only0")`: · RTI: ln(RTI) = ln[(e₁/t₁)/(e₀/t₀)], con error estándar √(1/e₁ + 1/e₀) del logaritmo (eventos de Poisson). · DTI: e₁/t₁ − e₀/t₀, con error estándar √(e₁/t₁² + e₀/t₀²). Si algún brazo tiene cero eventos se suma 0,5 a los dos conteos (no al tiempo), igual que en la tabla 2×2. Un estudio sin eventos en ningún brazo no informa la razón y queda fuera del combinado (ver `excluidoPorSinEventos`). `null` si los datos no son válidos (tiempos ≤ 0, eventos negativos o no finitos).

### `efectoDesdeRazonDeMedias`

```ts
export function efectoDesdeRazonDeMedias(mediaTrat: number, deTrat: number, nTrat: number, mediaCtrl: number, deCtrl: number, nCtrl: number): EfectoCalculado | null
```

Medias (2 grupos) → **razón de medias** (Cochrane Handbook 6.5.1.3; Friedrich, Adhikari y Beyene 2008): ln(RM) = ln(m₁/m₂), con varianza (de primer orden) sd₁²/(n₁·m₁²) + sd₂²/(n₂·m₂²). Solo tiene sentido con medidas que no pueden ser negativas (las dos medias deben ser positivas). No incluye la corrección de sesgo de muestra pequeña de `metafor::escalc("ROM")` (Hedges y cols. 1999), que cambia el logaritmo en una fracción de 1/n: la del Handbook es la fórmula simple. `null` si no son válidos.

### `efectoDesdeCorrelacion`

```ts
export function efectoDesdeCorrelacion(r: number, n: number): EfectoCalculado | null
```

Coeficiente de correlación r de un estudio → su efecto en la **escala z de Fisher**: z = atanh(r), con error estándar 1/√(n − 3) (`metafor::escalc("ZCOR")`). El IC se calcula en z y se devuelve a la escala de r con tanh. Necesita n > 3 y −1 < r < 1. `se` está en la escala z, como el de las demás medidas en su escala de análisis.

### `efectoDesdeOEV`

```ts
export function efectoDesdeOEV(oMenosE: number, varianza: number): EfectoCalculado | null
```

Tiempo hasta el evento: **O−E y V** (observados menos esperados y su varianza, como los da una prueba de rangos logarítmicos) → razón de riesgos instantáneos HR (Cochrane Handbook 6.8.2 y 10.9; Peto): ln(HR) = (O−E)/V, con error estándar 1/√V. Combinar estos ln(HR) con la varianza inversa en un modelo fijo equivale exactamente al método de Peto, Σ(O−E)/ΣV. `null` si V ≤ 0.

### `efectoDesdeCruzado`

```ts
export function efectoDesdeCruzado(mediaTrat: number, mediaCtrl: number, n: number, sdDiff: number, sdTrat: number, sdCtrl: number, correlacion: number): EfectoCalculado | null
```

Ensayo **cruzado (crossover)**, desenlace continuo → diferencia de medias PAREADA (Cochrane Handbook cap. 23.2.7.1). A diferencia de `efectoDesdeMedias` (dos grupos INDEPENDIENTES), un cruzado mide al MISMO grupo de `n` participantes bajo ambas condiciones — usar la fórmula de grupos independientes aquí ignora la correlación entre periodos y calcula un error estándar incorrecto (el "error de unidad de análisis" del cap. 23), sea de más o de menos según el signo de la correlación real.

### `efectoDeDisenoCluster`

```ts
export function efectoDeDisenoCluster(n: number, numClusters: number, icc: number): number | null
```

**Efecto de diseño** de un ensayo clúster-aleatorizado (Cochrane Handbook cap. 23.1): infla el error estándar por `√(efecto de diseño)` en vez de dividir N/eventos por él ("effective sample size" method) — el propio Handbook los llama equivalentes (cap. 23.1.5, "a slightly more flexible approach"), y este evita redondear N a un entero, funcionando con la MISMA fórmula tanto para desenlaces binarios como continuos (a diferencia del método de N efectiva, que divide también los eventos SOLO para binarios).

### `esTablaSinInformacion`

```ts
export function esTablaSinInformacion(crudos: DatosCrudos | undefined): boolean
```

Tabla 2×2 que no informa NINGÚN efecto relativo: sin eventos en los dos brazos, o con el evento en todos los participantes de los dos brazos (sin «no eventos»). Es el criterio de `drop00` de metafor y el que ya usa Mantel-Haenszel en aNova. `false` si no hay tabla 2×2 completa y numérica.

### `excluidoPorSinEventos`

```ts
export function excluidoPorSinEventos(e: Estudio): boolean
```

¿Queda este estudio fuera del OR/RR/PR por no tener eventos en ningún brazo? Cochrane Handbook 10.4.4.2: «Studies with no events in either arm … do not provide any information about the direction of the effect and are excluded from odds ratio or risk ratio analyses». No aplica a la diferencia de riesgos, donde sí aportan (RD = 0). Solo se detecta con datos crudos 2×2; un estudio capturado solo con su efecto e IC no se puede reconocer.

### `divisorControl`

```ts
export function divisorControl(dividirControlEn: number | undefined): number
```

Divisor válido del grupo control compartido (Handbook 23.3.4): un entero ≥ 2, o 1 (sin dividir) si no hay división indicada o el valor no es válido. Lo usan el cálculo del efecto y todas las pruebas que releen los datos crudos (Mantel-Haenszel, Peto, Peters, riesgo basal) para que el mismo control no cuente dos veces en ninguna.

### `confianzaHeterogeneidad`

```ts
export function confianzaHeterogeneidad(resultado: ResultadoMetaAnalisis): ICHeterogeneidad | null
```

**Intervalos de confianza del 95% para τ², τ e I² por el método Q-profile** (Viechtbauer 2007). Invierte el estadístico Q generalizado `Q_gen(τ²)=Σ(yᵢ−ȳ)²/(vᵢ+τ²)` (que sigue χ²_{k−1}): los límites de τ² son los τ² que igualan Q_gen a los cuantiles χ² 0.975 y 0.025. I² se deriva con la varianza intra-estudio típica `s²=(k−1)Σw/((Σw)²−Σw²)`. Réplica de `metafor::confint` / `meta` (método por defecto). `null` si <2 estudios.

### `ajustarMetaRegresion`

```ts
export function ajustarMetaRegresion(puntos: readonly PuntoEfecto[], X: number[][], nombres: string[]): ResultadoMetaRegresionMulti | null
```

**Meta-regresión de efectos mixtos** con una matriz de diseño arbitraria `X` (primera columna = intercepto), τ² residual por **REML** (réplica de `metafor::rma(mods=…, method="REML")`). Soporta moderadores múltiples y categóricos (codificados como columnas dummy por quien construye `X`). Devuelve coeficientes con IC, el test ómnibus QM y R². `null` si hay menos estudios que parámetros o el diseño es singular.

### `metaRegresionMulti`

```ts
export function metaRegresionMulti(estudios: Estudio[], medida: Medida, moderadores: EspecModerador[]): ResultadoMetaRegresionMulti | null
```

Meta-regresión desde los campos de los estudios: construye la matriz de diseño (intercepto + moderadores continuos centrados + categóricos codificados como dummies, referencia = primer nivel alfabético) y llama a `ajustarMetaRegresion` (REML). Excluye estudios sin todos los moderadores seleccionados. `null` si no hay datos suficientes (menos estudios que parámetros).

### `calcularMetaAnalisis`

```ts
export function calcularMetaAnalisis(estudios: Estudio[], medida: Medida, metodoTau2: MetodoTau2 = "DL"): ResultadoMetaAnalisis | null
```

**Motor del meta-análisis.** Calcula efecto combinado, heterogeneidad, IC, etc. para una medida.

### `agruparPorMedida`

```ts
export function agruparPorMedida(estudios: Estudio[], metodoTau2: MetodoTau2 = "DL"): ResultadoMetaAnalisis[]
```

**Orquestador del meta-análisis.** Calcula meta-análisis para TODAS las medidas presentes en los estudios.

### `metaAnalisisAcumulativo`

```ts
export function metaAnalisisAcumulativo(estudios: Estudio[], medida: Medida, modelo: "fijo" | "aleatorio", metodoTau2: MetodoTau2 = "DL"): FilaAcumulativa[] | null
```

**Meta-análisis acumulativo**: ordena los estudios (por defecto por año, ascendente) y recalcula el efecto combinado añadiendo un estudio a la vez, mostrando cómo se ha ido acumulando la evidencia. Reusa el motor ya validado (`calcularMetaAnalisis`) sobre cada prefijo; el primer paso (k=1) es el estudio solo. Replica `metafor::cumul`. `null` si hay <2 estudios válidos con esta medida.

### `analisisSensibilidad`

```ts
export function analisisSensibilidad(estudios: Estudio[], medida: Medida, modelo: "fijo" | "aleatorio", metodoTau2: MetodoTau2 = "DL"): FilaSensibilidad[] | null
```

Análisis de sensibilidad «dejar uno fuera»: repite el meta-análisis omitiendo cada estudio por turno, para ver si un solo estudio determina el resultado global. Requiere ≥ 3 estudios válidos (así quedan ≥ 2 al quitar uno). Usa el mismo motor y el modelo fijo o aleatorio elegido.

### `testEgger`

```ts
export function testEgger(resultado: ResultadoMetaAnalisis): ResultadoEgger | null
```

**Test de Egger:** regresión lineal para detectar sesgo de publicación en funnel plots.

### `testBegg`

```ts
export function testBegg(resultado: ResultadoMetaAnalisis): ResultadoBegg | null
```

**Test de Begg:** correlación de Kendall entre efecto y varianza.

### `testPeters`

```ts
export function testPeters(estudios: Estudio[], medida: Medida): ResultadoPeters | null
```

**Test de Peters:** regresión robusta para sesgo de publicación (alternativa a Egger).

### `metaAnalisisMH`

```ts
export function metaAnalisisMH(estudios: Estudio[], medida: Medida): ResultadoMH | null
```

Meta-análisis por el método de **Mantel-Haenszel** (efecto fijo) para datos binarios 2×2. Preferible al inverso-varianza cuando hay **eventos raros o muestras pequeñas**, porque no depende de varianzas intra-estudio de muestra grande y maneja bien las celdas cero (sin corrección de continuidad). - OR: `OR_MH = Σ(aᵢdᵢ/Nᵢ) / Σ(bᵢcᵢ/Nᵢ)`; varianza de ln(OR) por Robins-Breslow-Greenland. - RR: `RR_MH = Σ(aᵢn2ᵢ/Nᵢ) / Σ(cᵢn1ᵢ/Nᵢ)`; varianza de ln(RR) por Greenland-Robins. - RD (RAR): `RD_MH = Σwᵢ(RDᵢ)/Σwᵢ`, `wᵢ=n1ᵢn2ᵢ/Nᵢ`; varianza por Greenland-Robins. Fuentes: Cochrane Handbook §10.4 (Deeks & Higgins); Robins, Breslow & Greenland (1986); verificado contra `metafor::rma.mh` (Rothman et al. 2008). `null` si <2 estudios binarios con tabla 2×2 válida.

### `metaAnalisisPeto`

```ts
export function metaAnalisisPeto(estudios: Estudio[]): ResultadoPeto | null
```

Meta-análisis por el método de **Peto** (efecto fijo, solo OR) para datos binarios 2×2. Recomendado con **eventos raros y grupos de tamaño similar** (Cochrane Handbook §10.4.2). Por estudio usa la diferencia observado−esperado y la varianza hipergeométrica: Eᵢ = n1ᵢ(aᵢ+cᵢ)/Nᵢ ; Vᵢ = n1ᵢ·n2ᵢ·(aᵢ+cᵢ)·(bᵢ+dᵢ) / (Nᵢ²·(Nᵢ−1)) ; (O−E)ᵢ = aᵢ − Eᵢ ln(OR_Peto) = Σ(O−E)/ΣV ; var = 1/ΣV ; Q = Σ(O−E)²/V − (Σ(O−E))²/ΣV No necesita corrección de continuidad (maneja celdas cero). Verificado contra `metafor::rma.peto`. `null` si <2 estudios OR con tabla 2×2 utilizable. Advertencia: Peto se sesga con OR lejos de 1 o grupos desbalanceados; en ese caso preferir Mantel-Haenszel o inverso-varianza.

### `trimAndFill`

```ts
export function trimAndFill(resultado: ResultadoMetaAnalisis, modelo: "fijo" | "aleatorio", metodoTau2: MetodoTau2 = "DL"): ResultadoTrimFill | null
```

**Trim-and-Fill:** método de Duval & Tweedie para ajustar sesgo de publicación en funnel plots.

### `analisisSubgrupos`

```ts
export function analisisSubgrupos(estudios: Estudio[], medida: Medida, campo: CampoSubgrupo, modelo: "fijo" | "aleatorio", metodoTau2: MetodoTau2 = "DL"): ResultadoSubgrupos | null
```

Análisis de subgrupos: combina los estudios dentro de cada nivel de un moderador categórico y prueba si los subgrupos difieren con la Q entre subgrupos `Q_bet = Σ Wg(θg − θ̄)²` (χ², gl = G−1), como RevMan. Requiere ≥2 subgrupos con datos válidos.

### `baujat`

```ts
export function baujat(resultado: ResultadoMetaAnalisis): PuntoBaujat[] | null
```

Diagrama de Baujat (Baujat et al. 2002): por cada estudio, X = su contribución a la Q de heterogeneidad (`w(y−θF)²`) e Y = su influencia en el efecto combinado (`(θF − θF₍₋ᵢ₎)² · Σw₍₋ᵢ₎`), en efectos fijos. Los estudios arriba-derecha son los que más conviene revisar. Requiere ≥3 estudios.

### `diagnosticosInfluencia`

```ts
export function diagnosticosInfluencia(resultado: ResultadoMetaAnalisis, modelo: "fijo" | "aleatorio", metodoTau2: MetodoTau2 = "DL"): DiagnosticosInfluencia | null
```

**Diagnóstico de influencia** (Viechtbauer y Cheung 2010, como `metafor::influence`): para cada estudio, cuánto cambian el efecto combinado y la heterogeneidad al quitarlo, y cuán atípico es respecto al resto. Todo se calcula en la escala de análisis (log para razones). Con modelo aleatorio el τ² se vuelve a estimar sin cada estudio con el mismo método. Los umbrales para marcar un estudio como influyente son los convencionales del propio metafor (son orientativos, no una prueba): |DFFITS| > 3·√(1/(k−1)), Cook > mediana de χ²₁ (0,455), hat > 3/k, |DFBETAS| > 1. `null` con menos de 3 estudios.

### `galbraith`

```ts
export function galbraith(resultado: ResultadoMetaAnalisis, modelo: "fijo" | "aleatorio"): GraficoGalbraith | null
```

**Gráfico de Galbraith (radial)**: cada estudio se dibuja como (precisión, efecto estandarizado) = (1/√(vᵢ+τ²), yᵢ/√(vᵢ+τ²)). Con el efecto común μ̂ como pendiente, un estudio compatible con él cae a menos de 1,96 unidades verticales de la recta por el origen. Igual que `metafor::radial`: con el modelo aleatorio se usa el τ² del ajuste. `null` con menos de 2 estudios.

### `evaluarGrade`

```ts
export function evaluarGrade(estudios: Estudio[], medida: Medida, señales: SeñalesGrade): ResultadoGrade | null
```

Sugiere la certeza GRADE de un desenlace: arranca en Alta (mayoría ECA) o Baja (observacional) y ajusta por los 5 dominios de descenso (riesgo de sesgo, inconsistencia, imprecisión, evidencia indirecta, sesgo de publicación) y, si arrancó en Baja, los factores de aumento (efecto grande automático; dosis-respuesta y confusión residual quedan siempre a tu juicio). Es una sugerencia editable, con el motivo de cada dominio — no un veredicto final.

### `calcularEfectoAbsoluto`

```ts
export function calcularEfectoAbsoluto(estudios: Estudio[], medida: Medida, estimacionPooled: number): EfectoAbsoluto | null
```

Convierte un efecto relativo (OR o RR) en efecto absoluto — lo que GRADE llama "anticipated absolute effects" y que Cochrane Handbook exige en toda tabla Summary of Findings. Un OR de 0.66 no dice si el beneficio es grande o trivial; "80 eventos menos por cada 1,000 pacientes" sí.

### `valorNulo`

```ts
export function valorNulo(escala: EscalaAnalisis): number | null
```

Valor nulo (sin efecto) de una escala: 1 para la escala log (OR, RR, HR), 0 para la lineal y la z de Fisher; `null` si la escala no tiene un nulo con sentido (p. ej. proporciones).

### `esSignificativo`

```ts
export function esSignificativo(icInferior: number, icSuperior: number, nulo: number | null): boolean
```

¿El intervalo de confianza excluye el valor nulo? (es decir, ¿el efecto combinado es estadísticamente significativo?).

### `validarEfectoIC`

```ts
export function validarEfectoIC(estudio: Pick<Estudio, "medida" | "valor" | "icInferior" | "icSuperior">): string | null
```

Validación en vivo de una fila de la tabla de extracción: devuelve `null` mientras la fila esté sin tocar (para no molestar en una fila nueva en blanco) y, si no, un mensaje de error específico.

### `calcularSROC`

```ts
export function calcularSROC(estudios: Estudio[], modelo: "fijo" | "aleatorio" = "aleatorio", metodoTau2: MetodoTau2 = "DL", opciones: { bivariado?: boolean } = {}): ResultadoSROC | null
```

**SROC (Summary Receiver Operating Characteristic)** — precisión diagnóstica combinada. Lee `Estudio.datosDiagnosticos` (VP/FP/FN/VN), independiente de `medida`/`valor`.

## `src/moses-littenberg.ts`

### `intervaloWilson`

```ts
export function intervaloWilson(x: number, n: number):
```

Intervalo de Wilson de una proporción x/n (más estable que el de Wald cuando x = 0 o x = n). Wilson 1927; Newcombe 1998.

### `puntoDiagnostico`

```ts
export function puntoDiagnostico(vp: number, fp: number, fn: number, vn: number): PuntoDiagnostico | null
```

Tabla de un estudio diagnóstico (VP, FP, FN, VN) → ln DOR con su error estándar y S = logit(sensibilidad) + logit(FPR), las dos variables de la regresión de Moses-Littenberg; suma 0,5 a las cuatro celdas si alguna es cero. `null` si los conteos no son válidos (negativos, no finitos, o sin enfermos o sin sanos).

### `regresionMosesLittenberg`

```ts
export function regresionMosesLittenberg(puntos: readonly PuntoDiagnostico[]): RegresionMosesLittenberg | null
```

Mínimos cuadrados ponderados de D sobre S con peso 1/Var(D). Con la matriz de información [[Σw, Σws], [Σws, Σws²]], su inversa se escribe con el determinante Δ = Σw·Σws² − (Σws)²: Var(α) = Σws²/Δ, Var(β) = Σw/Δ. `null` si Δ ≈ 0 (todos los S iguales).

### `curvaMosesLittenberg`

```ts
export function curvaMosesLittenberg(reg: { a: number; b: number }, n = 61):
```

Curva SROC en `n` puntos de FPR entre 0,001 y 0,999 (escala logit lineal en L_F); `null` si |β| ≥ 1 (la curva no está definida).

## `src/sroc-bivariado.ts`

### `logVerosimilitudBinomial`

```ts
export function logVerosimilitudBinomial(tablas: TablaDiagnostica[], p: ParametrosBivariado, nodosPorDimension = 15): number
```

Log-verosimilitud marginal del modelo binomial-normal (con los coeficientes binomiales).

### `sensibilidadEnCurva`

```ts
export function sensibilidadEnCurva(p: ParametrosBivariado, fpr: number): number
```

Curva SROC implícita en el modelo: en el espacio logit-logit es la recta que pasa por el punto resumen con pendiente σA/σB (la «generalización de Rutter-Gatsonis», coincide con el HSROC cuando hay equivalencia con el bivariado: Harbord 2007): logit(sensibilidad) = μA + (σA/σB)·(logit(FPR) − logit(FPR̄)), con logit(FPR̄) = −μB.

### `curvaSROC`

```ts
export function curvaSROC(p: ParametrosBivariado, desde: number, hasta: number, n = 80): PuntoROC[]
```

Curva SROC en `n` valores de FPR entre `desde` y `hasta`, equiespaciados en la escala logit (más puntos donde la curva se dobla).

### `areaBajoCurva`

```ts
export function areaBajoCurva(p: ParametrosBivariado): number
```

Área bajo la curva SROC, ∫₀¹ TPR(FPR) dFPR, por la regla del trapecio en la escala logit (FPR = σ(u), dFPR = σ(u)(1−σ(u)) du) entre u = −40 y 40 con 8000 intervalos: error ≈ 1e-10. Incluye la extrapolación fuera del rango de FPR observado.

### `elipseROC`

```ts
export function elipseROC(centro: [number, number], cov: [[number, number], [number, number]], nivel = 0.95, n = 100): PuntoROC[]
```

Elipse de nivel `nivel` en el plano logit (centrada en `centro`, con matriz de covarianza `cov`) pasada al espacio ROC (FPR, sensibilidad): d' cov⁻¹ d = χ²₂(nivel) = −2 ln(1 − nivel). Los puntos están en el orden del ángulo. Coordenadas logit: (logit sens, logit esp); FPR = σ(−logit esp).

### `regionesROC`

```ts
export function regionesROC(r: ResultadoBivariado, nivel = 0.95):
```

Región de confianza (para el punto resumen) y de predicción (para un estudio nuevo) en el espacio ROC.

### `tablasUtilizables`

```ts
export function tablasUtilizables(tablas: TablaDiagnostica[]): TablaDiagnostica[]
```

Valida las tablas: conteos finitos, no negativos y con enfermos y sanos en cada estudio. Devuelve las utilizables.

### `ajustarBivariado`

```ts
export function ajustarBivariado(tablasEntrada: TablaDiagnostica[], opc: OpcionesAjuste = {}): ResultadoBivariado | null
```

Ajusta el modelo bivariado. `null` con menos de 4 estudios utilizables o si el ajuste no se puede calcular (matriz singular). Con «binomial», si la optimización no converge cae en la aproximación normal e indica el motivo en `advertencias`.

### `ajustarBivariadoConCache`

```ts
export function ajustarBivariadoConCache(tablas: TablaDiagnostica[], metodo: MetodoBivariado = "binomial"): ResultadoBivariado | null
```

Igual que `ajustarBivariado`, pero recuerda los últimos ajustes por (método + tablas): el Paso 7, el manuscrito y la cobertura piden el mismo SROC varias veces y el ajuste binomial tarda ~1 s con 40 estudios.

## `src/multinivel.ts`

### `inferenciaRobusta`

```ts
export function inferenciaRobusta(datos: EfectoDependiente[], sigma2Estudio: number, sigma2Efecto: number, opc: OpcionesMultinivel = {}):
```

Resultado de la parte de varianza robusta con los componentes de varianza dados (los que ajusta `ajustarMultinivel`, o los de otro programa para comprobarla). Expuesta para validar el estimador sándwich de forma aislada.

### `ajustarMultinivel`

```ts
export function ajustarMultinivel(datos: EfectoDependiente[], opc: OpcionesMultinivel = {}): ResultadoMultinivel | null
```

Ajusta el modelo CHE por REML y calcula la inferencia robusta. `null` con menos de 3 efectos, con menos de 2 estudios o si el ajuste no se puede calcular. La escala es la de los `y` (log para OR/RR/HR…); la conversión a la escala natural la hace quien llama.

### `sensibilidadRho`

```ts
export function sensibilidadRho(datos: EfectoDependiente[], opc: OpcionesMultinivel = {}, rhos: readonly number[] = RHO_SENSIBILIDAD):
```

Repite el ajuste multinivel con varias correlaciones supuestas entre efectos (por defecto 0; 0,2; 0,4; 0,6 y 0,8) para ver cuánto depende el resultado de esa suposición.

## `src/glmm.ts`

### `ajustarGLMMProporcion`

```ts
export function ajustarGLMMProporcion(datos: DatosProporcion[]): ResultadoGLMM | null
```

Proporción combinada con el modelo binomial-normal (logit). `null` con menos de 2 estudios válidos.

### `ajustarGLMMOdds`

```ts
export function ajustarGLMMOdds(datos: DatosTabla[]): ResultadoGLMM | null
```

Odds ratio combinado con el modelo hipergeométrico-normal (verosimilitud condicional exacta). `null` con menos de 2 estudios con información.

### `logVerosimilitudProporcion`

```ts
export function logVerosimilitudProporcion(datos: DatosProporcion[], mu: number, tau: number, nodosPorEstudio = NODOS): number
```

Log-verosimilitud marginal (con los coeficientes binomiales) de las proporciones en (μ, τ); expuesta para validar la cuadratura.

### `logVerosimilitudOdds`

```ts
export function logVerosimilitudOdds(datos: DatosTabla[], mu: number, tau: number, nodosPorEstudio = NODOS): number
```

Log-verosimilitud marginal de las tablas 2×2 (hipergeométrica no central) en (ln OR, τ); expuesta para validar la cuadratura.

## `src/red.ts`

### `calcularNMA`

```ts
export function calcularNMA(estudiosRed: EstudioRed[], medida: MedidaRed, modelo: "fijo" | "aleatorio" = "aleatorio"): ResultadoNMA | null
```

Calcula el meta-análisis en red para la medida indicada.

### `calcularPScores`

```ts
export function calcularPScores(resultado: ResultadoNMA, menorEsMejor = true): PScoreTratamiento[]
```

Ranking por P-score (Rücker y Schwarzer 2015): para cada par, la probabilidad de que i sea mejor que j es Φ(−z) si un valor menor es mejor (Φ(z) si es mayor), con z = diferencia/EE en la escala de análisis; el P-score de i es su promedio contra los demás tratamientos (1 = seguramente el mejor; 0 = seguramente el peor). Es la versión frecuentista de SUCRA, sin remuestreo.

### `calcularNodeSplittingSIDDE`

```ts
export function calcularNodeSplittingSIDDE(estudiosRed: EstudioRed[], medida: MedidaRed, modelo: "fijo" | "aleatorio" = "aleatorio"): NodeSplitEstimacion[]
```

Node-splitting por SIDDE (Efthimiou y cols. 2019): para cada comparación con estimado directo, se quitan los estudios completos que comparan directamente el par (un estudio multibrazo que lo incluya también sale) y se reajusta la red con el resto; la estimación de red de ese reajuste del par es el estimado indirecto. El estimado directo es el del back-calculation, que no depende del método.

## `src/kappa.ts`

### `interpretarKappa`

```ts
export function interpretarKappa(k: number): string
```

Interpretación cualitativa de κ según Landis & Koch (1977).

### `cohenKappa`

```ts
export function cohenKappa(pares: [string, string][]): ResultadoKappa | null
```

κ de Cohen para 2 evaluadores (sin ponderar). `pares`: una fila por sujeto con las dos calificaciones [evaluadorA, evaluadorB].

### `fleissKappa`

```ts
export function fleissKappa(sujetos: string[][]): ResultadoKappa | null
```

κ de Fleiss para n evaluadores (sin ponderar). `sujetos`: una fila por sujeto con las calificaciones de todos los evaluadores (el mismo nº por sujeto).

## `src/conversiones-datos.ts`

### `sdDesdeEE`

```ts
export function sdDesdeEE(ee: number, n: number): number | null
```

DE de un grupo desde el error estándar de su media: SD = SE · √N (6.5.2.2).

### `sdDesdeICMedia`

```ts
export function sdDesdeICMedia(inf: number, sup: number, n: number, opc: OpcionesIC = {}):
```

DE de un grupo desde el IC de su media (6.5.2.2): SE = (límite superior − inferior) / (2·q), SD = SE·√N, con q el cuantil de la t con N − 1 grados de libertad (o de la normal; 1,96 al 95 %: «3,92 de ancho»). El Handbook avisa de que con menos de unos 60 participantes por grupo el IC se calcula con la t.

### `sdDesdeEEDiferencia`

```ts
export function sdDesdeEEDiferencia(se: number, nE: number, nC: number): number | null
```

DE intragrupo desde el error estándar de la diferencia de medias: SD = SE / √(1/N₁ + 1/N₂) (6.5.2.3).

### `eeDiferenciaDesdeIC`

```ts
export function eeDiferenciaDesdeIC(inf: number, sup: number, nE: number, nC: number, opc: OpcionesIC = {}): number | null
```

Error estándar de la diferencia de medias desde su IC (6.5.2.3), con t de N₁ + N₂ − 2 grados de libertad (o normal).

### `eeDiferenciaDesdeT`

```ts
export function eeDiferenciaDesdeT(diferencia: number, t: number): number | null
```

Error estándar de la diferencia de medias desde el estadístico t: SE = |DM| / |t| (6.5.2.3).

### `eeDiferenciaDesdeP`

```ts
export function eeDiferenciaDesdeP(diferencia: number, p: number, nE: number, nC: number): number | null
```

Error estándar de la diferencia de medias desde el valor de P de dos colas (6.5.2.3): primero la t que corresponde a ese P con N₁ + N₂ − 2 grados de libertad, luego SE = |DM| / t. Un P aproximado («p < 0,05», «p = 0,01») da un SE aproximado. `null` si P no está en (0, 1) o la diferencia es 0 (con DM = 0 el SE no se puede recuperar).

### `eeDesdeICAbsoluto`

```ts
export function eeDesdeICAbsoluto(inf: number, sup: number, nivel = 0.95): number | null
```

Error estándar desde el IC de una medida ABSOLUTA con distribución normal (DM, DME, diferencia de riesgos), 6.3.1: SE = (sup − inf) / 3,92 al 95 %, 3,29 al 90 %, 5,15 al 99 % (es decir, 2·z).

### `eeLogDesdeICRazon`

```ts
export function eeLogDesdeICRazon(inf: number, sup: number, nivel = 0.95): number | null
```

Error estándar del logaritmo de una razón (RR, OR, HR) desde su IC, 6.3.2: (ln sup − ln inf) / 3,92 al 95 %.

### `eeDesdePAbsoluto`

```ts
export function eeDesdePAbsoluto(estimacion: number, p: number): number | null
```

Error estándar desde el valor de P de dos colas de una medida ABSOLUTA (6.3.1): z desde P y SE = |estimación| / z.

### `icDesdeEEAbsoluto`

```ts
export function icDesdeEEAbsoluto(estimacion: number, ee: number, nivel = 0.95):
```

IC (inferior, superior) de una medida ABSOLUTA desde su estimación y su error estándar, con la normal al nivel dado.

### `icDesdeEELogRazon`

```ts
export function icDesdeEELogRazon(razon: number, eeLog: number, nivel = 0.95):
```

IC (inferior, superior) de una RAZÓN desde la razón y el error estándar de su logaritmo: exp(ln(razón) ± z · EE).

### `eeLogDesdePRazon`

```ts
export function eeLogDesdePRazon(razon: number, p: number): number | null
```

Error estándar del log de una razón desde su P de dos colas (6.3.2): SE = |ln(razón)| / z.

### `combinarGrupos`

```ts
export function combinarGrupos(a: GrupoResumen, b: GrupoResumen): GrupoResumen | null
```

Combina dos grupos en uno solo (Handbook Tabla 6.5.a), como si nunca se hubieran dividido: N = N₁ + N₂, media ponderada y SD² = [ (N₁−1)·SD₁² + (N₂−1)·SD₂² + N₁N₂/(N₁+N₂) · (M₁² + M₂² − 2·M₁·M₂) ] / (N₁ + N₂ − 1). Con más de dos grupos se aplica en secuencia. Sirve, por ejemplo, para fundir dos brazos experimentales de un mismo ensayo.

### `sdCambioDesdeCorrelacion`

```ts
export function sdCambioDesdeCorrelacion(sdBase: number, sdFinal: number, correlacion: number): number | null
```

DE del cambio desde la línea base imputada con una correlación (6.5.2.8): √(SD_base² + SD_final² − 2·Corr·SD_base·SD_final).

### `correlacionDesdeSdCambio`

```ts
export function correlacionDesdeSdCambio(sdBase: number, sdFinal: number, sdCambio: number): number | null
```

Correlación base–final que implican tres DE de un mismo grupo (6.5.2.8): (SD_base² + SD_final² − SD_cambio²) / (2·SD_base·SD_final).

### `sdDesdeIQRNormal`

```ts
export function sdDesdeIQRNormal(q1: number, q3: number): number | null
```

DE aproximada desde el rango intercuartílico cuando la distribución es normal y la muestra grande (6.5.2.5): IQR / 1,35.

### `mediaSdWan`

```ts
export function mediaSdWan(e: EntradaWan): ResultadoWan | null
```

Media y DE estimadas desde la mediana y el rango y/o los cuartiles (Wan et al. 2014, que el Handbook 6.5.2.5 y 6.5.2.9 recomienda). Tres escenarios según lo que reporte el estudio: · mínimo, mediana, máximo: media ≈ (a + 2m + b)/4 + (a − 2m + b)/(4n) (Hozo et al., ec. 2 del artículo); DE ≈ (b − a) / ξ(n). · cuartiles: media ≈ (q₁ + m + q₃)/3; DE ≈ (q₃ − q₁) / η(n). · cinco números: media ≈ (a + 2q₁ + 2m + 2q₃ + b)/8; DE ≈ ½ [ (b − a)/ξ(n) + (q₃ − q₁)/η(n) ]. Suponen datos normales; con datos asimétricos son una aproximación (ligeramente sesgada en el artículo) y con tamaños muy pequeños son poco fiables. Si se dan más datos de los necesarios se usa el escenario de los cinco números solo si están todos; ante cualquier otra combinación devuelve `null`. Usa las fórmulas finales del artículo (aproximación de Blom) para cualquier n; `meta::metamean` usa además los valores tabulados del artículo para n ≤ 50, de modo que con n pequeños los resultados pueden diferir unas décimas de por ciento (ver VALIDACION_ESTADISTICA.md).

### `dmeDesdeOR`

```ts
export function dmeDesdeOR(or: number, eeLogOR: number):
```

Odds ratio → diferencia de medias estandarizada (10.6): DME = ln(OR) · √3/π; el EE se multiplica por la misma constante.

### `orDesdeDME`

```ts
export function orDesdeDME(dme: number, ee: number):
```

Diferencia de medias estandarizada → odds ratio (10.6): ln(OR) = DME · π/√3; el EE del log OR se multiplica por la misma constante.

## `src/efecto-crudos.ts`

### `tieneEventoCero`

```ts
export function tieneEventoCero(crudos: DatosCrudos): boolean
```

Cochrane Handbook cap. 19.5.2: un "0" reportado puede ser un evento realmente ausente (con vigilancia activa) o solo ausencia de reporte — solo aplica a los tipos con conteo de eventos.

### `faltaJustificarCero`

```ts
export function faltaJustificarCero(crudos: DatosCrudos): boolean
```

¿Hay un cero reportado (tabla 2×2 o proporción) sin justificación escrita? Handbook 19.5.2: un 0 puede ser ausencia real del evento o ausencia de reporte, y debe justificarse.

### `calcularEfectoCrudos`

```ts
export function calcularEfectoCrudos(medida: Medida, crudos?: DatosCrudos, dividirControlEn?: number): EfectoCalculado | null
```

Calcula la medida de efecto (OR, RR, RAR, PR, DM, DME, Proporción) a partir de datos crudos.

### `compararVerificacion`

```ts
export function compararVerificacion(e: Estudio):
```

Compara la extracción principal contra la verificación de un segundo revisor (si existe).
