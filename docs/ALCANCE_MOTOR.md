# Alcance del motor estadístico: qué hay y qué falta

> *Nota: este documento procede del repositorio de la aplicación aNova Meta Tool. Las menciones a «Paso N», pantallas o paneles se refieren a su interfaz; los módulos del motor están en `src/`.*


> **Estado:** inventario del 6 de octubre de 2026 y su cierre el mismo día. Se hizo comparando el código con los índices de los
> capítulos 6 y 10 del Cochrane Handbook v6.5 (leídos en cochrane.org) y con lo que ofrecen RevMan, metafor y meta. **Alcance
> aprobado por el autor: A, B y C completos, solo la influencia y el gráfico de Galbraith de D, el modelo bivariado de precisión diagnóstica (E) y, después, los efectos dependientes con varianza robusta (F) y los modelos mixtos logísticos (G).** El detalle de cómo
> se validó cada pieza está en `VALIDACION_ESTADISTICA.md` §3.12–3.18 y de dónde sale cada una en `PROCEDENCIA.md`.
> Criterio de «completo»: lo que el Handbook recomienda para revisiones de intervenciones (cap. 6, 10, 23) más lo que ya ofrecía la
> aplicación para precisión diagnóstica y red. Todo lo demás se declara fuera de alcance, por escrito (sección 3).

## 1. Lo que hay

| Área | Qué hay | Handbook |
|---|---|---|
| Medidas binarias | RR, OR, PR, diferencia de riesgos (RAR), proporción de un brazo; Mantel-Haenszel y Peto; **estudios sin eventos en ningún brazo excluidos del OR/RR** | 10.4, 6.4 |
| Tasas | **Razón (RTI) y diferencia (DTI) de tasas de incidencia** con eventos y tiempo-persona | 6.7, 10.8 |
| Medidas continuas | DM, DME (g de Hedges), **razón de medias (RM)**, cruzados (DM pareada), clúster (efecto de diseño) | 10.5, 6.5.1.3, 23 |
| Correlaciones | **r, combinada con la z de Fisher** | — |
| Tiempo hasta el evento | HR por varianza inversa (log HR y su EE) y **desde O−E y V (equivalente al método de Peto)** | 10.9, 6.8 |
| Combinación | Efectos fijos y aleatorios; **modelos mixtos logísticos (hipergeométrico-normal para el OR, binomial-normal para proporciones)**; **efectos dependientes: modelo multinivel CHE con varianza robusta CR2**; τ² por DL, Paule-Mandel y REML (**con tolerancia relativa a la escala**); Hartung-Knapp; intervalo de predicción; IC de τ² e I² (Q-profile) | 10.3, 10.10 |
| Pocos estudios | Con 2–3 estudios el aviso **compara Wald con Knapp-Hartung** en vez de recomendar el segundo | 10.10.4.5 |
| Unidad de análisis | **Detección de ensayos repetidos y división del grupo control compartido** (también en MH, Peto, Peters y riesgo basal) | 6.2.4, 6.2.9, 23.3.4 |
| Conversiones de datos | **DE desde EE, IC, t o P; EE e IC desde el P; combinar grupos; DE del cambio y su correlación; media y DE desde mediana y cuartiles o rango (Wan); OR ↔ DME** | 6.3, 6.5.2, 10.6 |
| Exploración | Subgrupos, meta-regresión multivariable, sensibilidad (dejar uno fuera), acumulado, Baujat, **diagnóstico de influencia y gráfico de Galbraith** | 10.11, 10.14 |
| Sesgo de publicación | Egger, Begg (**con empates**), Peters (OR), trim-and-fill (L0) | cap. 13 |
| Precisión diagnóstica | **Modelo bivariado de Reitsma** (verosimilitud binomial exacta, con la aproximación normal REML como comparación), curva SROC, AUC y regiones de confianza y de predicción; Moses-Littenberg y agrupación por separado como exploratorios | DTA 10 |
| Otros métodos | meta-análisis en red (Rücker, P-scores, node-splitting, SIDDE); κ de Cohen y Fleiss (reescrito) | — |

## 2. Lo que se pidió y cómo quedó

### A. Riesgo metodológico: hecho

| # | Qué | Estado |
|---|---|---|
| A1 | Ensayos con más de dos grupos en un análisis por pares (el control contado dos veces) | ✅ Detección + decisión por fila (dividir el control en k, o «no lo comparte») + «por definir» en el manuscrito si falta decidir. Combinar los grupos experimentales es la otra salida (conversor). `test:ensayos-repetidos` |
| A2 | Estudios sin eventos en ningún brazo en OR/RR | ✅ Excluidos (también de la razón de tasas), declarado en Métodos, cobertura y avisos; la diferencia de riesgos los conserva. Referencia de R regenerada con `drop00`. `test:sin-eventos` |
| A3 | Pocos estudios (2–3) | ✅ El aviso compara Wald y Knapp-Hartung con los números reales |
| A4 | Observaciones repetidas del mismo estudio | ✅ Misma detección que A1 (el aviso distingue varios brazos de varios seguimientos; la decisión es de la persona) |

### B. Conversiones de datos: hecho

| # | Qué | Estado |
|---|---|---|
| B1 | DE desde EE, IC, t o P | ✅ `test:conversiones`, contra R (datos simulados con verdad conocida, ≈1e-13) |
| B2 | Media y DE desde mediana, cuartiles y rango | ✅ Wan y cols. 2014 (fórmulas finales del artículo, aproximación de Blom). Contra `meta::metamean`: media exacta, DE a ≈6e-4 con n > 50. Con n pequeño difiere de las tablas que usa `meta` (documentado). No se implementaron Luo 2018 ni Shi 2020 (mejoras posteriores; el Handbook cita Wan) |
| B3 | DE del cambio desde la línea base, con tabla de sensibilidad a la correlación | ✅ |
| B4 | Combinar grupos | ✅ Tabla 6.5.a, contra los datos completos |
| B5 | OR ↔ DME | ✅ Handbook 10.6 |
| — | Anotar qué se convirtió | ✅ Campo por estudio que Métodos declara |

### C. Medidas: hecho, con una decisión de no hacer

| # | Qué | Estado |
|---|---|---|
| C1 | Razón y diferencia de tasas | ✅ 40 conjuntos contra `escalc` IRR e IRD |
| C2 | Razón de medias | ✅ Fórmula del Handbook; difiere de `escalc("ROM")` en una corrección de sesgo O(1/n), documentada |
| C3 | Tiempo hasta el evento con O−E y V | ✅ Identidad exacta con Σ(O−E)/ΣV |
| C4 | **Doble arcoseno (Freeman-Tukey) para proporciones** | ❌ **No se implementa.** Schwarzer, Chemaitelly, Abu-Raddad y Rücker (2019, *Research Synthesis Methods*) mostraron que su transformación inversa puede dar resultados «seriamente engañosos» y proponen los modelos mixtos como alternativa (una réplica de 2025 la defiende: sigue siendo discutido). Se mantiene el logit con corrección de continuidad y se declara el límite |
| C5 | Correlaciones (z de Fisher) | ✅ 40 conjuntos contra `escalc` ZCOR |

### D. Análisis adicionales: solo lo acordado

| # | Qué | Estado |
|---|---|---|
| D1 | Diagnóstico de influencia | ✅ Diez cantidades por estudio contra `influence()` (≈1e-10); umbrales de su documentación |
| D2 | Gráfico de Galbraith | ✅ Contra `radial()` |
| D2 | Gráfico de L'Abbé | ⏸ no hecho |
| D3 | Pruebas de asimetría de Harbord y Rücker; embudo con contornos | ⏸ no hecho |
| D4 | Mejor y peor caso para datos faltantes | ⏸ no hecho |
| D5 | Otros estimadores de τ² (ML, SJ, HE…) | ⏸ no hecho (poca urgencia) |

### Arreglos del motor que salieron de validar lo nuevo

1. **REML con tolerancia absoluta** (se detenía antes de converger con efectos de magnitud pequeña) → relativa a la escala.
2. **Begg con empates** (τ-a, p exacta no válida) → τ-b y varianza de Kendall con empates.

### E. Precisión diagnóstica: modelo bivariado (hecho)

| # | Qué | Estado |
|---|---|---|
| E1 | Modelo bivariado de efectos aleatorios (Reitsma 2005; Chu y Cole 2006) con verosimilitud binomial exacta dentro de cada estudio, máxima verosimilitud por cuadratura de Gauss-Hermite adaptativa | ✅ Contra una implementación independiente en R y contra `lme4::glmer` (comprobación del modelo). Es el que recomienda el Handbook de precisión diagnóstica |
| E2 | Aproximación normal con REML (la de `mada::reitsma`) | ✅ Incluida como comparación y como respaldo declarado si el modelo exacto no se puede ajustar. Contra el óptimo exacto de R y contra `mada` |
| E3 | Curva SROC (Rutter-Gatsonis), área bajo la curva, regiones de confianza y de predicción | ✅ Contra `mada::sroc`, `AUC` y `ROCellipse`; la curva solo se dibuja en el rango de FPR observado |
| E4 | Moses-Littenberg, sensibilidad y especificidad por separado, DOR | ✅ Se conservan como análisis **exploratorio** (el Handbook no los recomienda para los intervalos); el manuscrito lo dice así |
| E5 | Parámetros en el límite (varianza ≈ 0, ρ = ±1) | ✅ Declarados, y los errores estándar se calculan con el parámetro fijado en el límite (coincide con R) |

### F. Efectos dependientes: multinivel con varianza robusta (hecho)

| # | Qué | Estado |
|---|---|---|
| F1 | Modelo CHE (σ² entre estudios y entre efectos, REML; correlación supuesta ρ entre errores de muestreo) | ✅ Contra `metafor::rma.mv` y `vcalc` (60 conjuntos). `test:multinivel` |
| F2 | Error estándar robusto CR2 con df de Satterthwaite (un coeficiente) | ✅ Contra `clubSandwich` a precisión de máquina; solución analítica y df exactos en diseño balanceado |
| F3 | Sensibilidad a ρ (0 a 0,8) en pantalla y manuscrito; la ρ elegida se guarda y se declara | ✅ |
| F4 | Contrastes de varios coeficientes (Hotelling), covariables en la interfaz, Hartung-Knapp aquí | ⏸ no hecho (el motor admite covariables y están validadas) |

### G. Modelos mixtos logísticos (hecho)

| # | Qué | Estado |
|---|---|---|
| G1 | OR: hipergeométrico-normal (verosimilitud condicional exacta, sin corrección de continuidad) | ✅ Contra implementación independiente en R, fuerza bruta (con `BiasedUrn`) y `metafor::rma.glmm` (CM.EL, donde converge) |
| G2 | Proporciones: binomial-normal | ✅ Contra implementación independiente, fuerza bruta y `rma.glmm` (PLO) |
| G3 | Estudios sin información (sin eventos en ningún brazo, o todos con evento), ensayos por clúster o cruzados, control dividido | ✅ Se declaran en pantalla y manuscrito; no se corrigen ni se redondean |
| G4 | Binomial con efectos fijos por estudio (UM.FS), RR, tasas, intervalo de predicción, ajuste para pocos estudios | ⏸ no hecho |

## 3. Fuera de alcance (declarado)

- Del modelo bivariado de precisión diagnóstica (hecho el 6 de octubre de 2026, §2 E): covariables (meta-regresión diagnóstica), comparación
  directa entre pruebas, el HSROC parametrizado por separado (su curva coincide con la del bivariado: Harbord 2007), ajuste bayesiano y valores
  predictivos con prevalencia.
- Modelos bayesianos (Handbook 10.13), multivariantes (varios desenlaces con su correlación estimada), dosis-respuesta (10.11.8), meta-regresión en
  red, CINeMA y matriz de contribución. (Los GLMM y el multinivel con varianza robusta ya están: §2 F y G.)
- Datos ordinales con odds proporcionales (6.6, 10.7).
- Doble arcoseno de Freeman-Tukey (ver C4).
- Los pendientes de D de la tabla anterior (L'Abbé, Harbord/Rücker, mejor y peor caso, otros estimadores de τ²).
- Diseños sin cálculo (ya declarados): cruzado con DME o binarios, clúster-cruzado, escalonados, aleatorizados individuales con agrupamiento.
- Luo (2018) y Shi (2020) para la media y la DE desde cuantiles (el Handbook cita Wan).

## 4. Cómo se hizo

Cada pieza se construyó **junto con su validación** (referencia de R en `validacion/`, prueba y mutaciones), escribiendo desde las fórmulas
publicadas y no leyendo código de otros programas (ver `PROCEDENCIA.md`). Resultado a esta fecha: **255.661 comparaciones contra R en 1.057 grupos
(705 conjuntos de datos), 0 fuera de tolerancia**; 45 suites, 24 pruebas E2E, 14 escenarios de auditoría sin hallazgos.
