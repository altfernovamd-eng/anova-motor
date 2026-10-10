// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * **Meta-análisis en red frecuentista** (método gráfico-teórico de Rücker): combina evidencia directa e indirecta entre ≥ 3 tratamientos.
 * Reescrito el 7 de octubre de 2026 desde las fórmulas publicadas (ver docs/PROCEDENCIA.md): la versión anterior se había escrito
 * confirmando detalles contra el código de `netmeta` (GPL); esta parte de las ecuaciones de los artículos, con otra organización, otras
 * rutas de cálculo (Laplaciano con un nodo a tierra, Cholesky) y otros nombres, y usa `netmeta` solo como caja negra de comparación.
 *
 * El modelo, en lo esencial:
 *  · Cada estudio con k brazos aporta contrastes entre sus brazos. Un contraste entre los brazos a y b tiene efecto yₐ − y_b y varianza
 *    vₐ + v_b (+ τ² en efectos aleatorios: cada brazo lleva un efecto aleatorio de varianza τ²/2).
 *  · Red como circuito eléctrico (Rücker 2012): tratamientos = nodos, contrastes = cables con conductancia (peso) w = 1/varianza.
 *    Matriz Laplaciana ponderada L (L_ii = Σ pesos incidentes, L_ij = −peso entre i y j). Los efectos de red son diferencias de los
 *    «potenciales» θ = L⁺ b, con b_i = Σ (± w·y) el flujo neto en el nodo i y L⁺ la pseudoinversa de Moore-Penrose.
 *    Var(θ_i − θ_j) = L⁺_ii + L⁺_jj − 2 L⁺_ij.
 *  · Estudios multibrazo (Rücker y Schwarzer 2014, Stat Med 33:4353): sus k brazos tienen k − 1 contrastes independientes, pero se
 *    reparten en k(k−1)/2 cables. Se «reducen los pesos» para que el conjunto reproduzca la covarianza real del estudio: con la matriz
 *    de covarianza Σ de los k − 1 contrastes respecto al primer brazo, la matriz de precisión Σ⁻¹ es el Laplaciano del estudio con el
 *    primer brazo a tierra, de donde w_ab = −(Σ⁻¹)_ab y w_0a = Σ_b (Σ⁻¹)_ab. Los autores demuestran que es equivalente al modelo
 *    estándar de «reducir la dimensión» (mínimos cuadrados generalizados con la covarianza de cada estudio).
 *  · Heterogeneidad e inconsistencia: Q = Σ w (y − ŷ)², con gl = Σ(kₛ − 1) − (n − 1). τ² por momentos (generalización de
 *    DerSimonian-Laird): como E[Q] = gl + τ²·tr[(I − H)·U·W] (H = B L⁺ Bᵀ W, U la estructura de covarianza de los efectos aleatorios,
 *    U_ef = ½·(B_e·B_f) entre contrastes del mismo estudio), τ² = max(0, (Q − gl)/tr[(I − H) U W]).
 *  · Q se descompone en heterogeneidad dentro de diseños e inconsistencia entre diseños (Krahn, Binder y König 2013, BMC MRM 13:35,
 *    ecs. 3, 4, 8 y 10): un «diseño» es el conjunto exacto de tratamientos de un estudio.
 *  · P-scores (Rücker y Schwarzer 2015, BMC MRM 15:58): promedio de las probabilidades Φ(±z) de que un tratamiento supere a cada otro.
 *  · Inconsistencia local: node-splitting por back-calculation (König, Krahn y Binder 2013, Stat Med 32:5414: proporción de evidencia
 *    directa P = Var(red)/Var(directa); indirecto = (red − P·directo)/(1 − P)) y por SIDDE (Efthimiou y cols. 2019, Stat Med 38:2992:
 *    se quitan los estudios que comparan directamente el par y se reajusta la red; el estimado de ese reajuste es el indirecto).
 *
 * Validado contra `netmeta` (dataset Stowe2010 y 40 redes sintéticas con estudios de tres brazos, DM y OR; `npm run validacion` y
 * `scripts/benchmark-red.mts`). Lógica pura, sin dependencias de la app salvo `funciones-especiales` y `numerica-compartida`.
 */
import type { BrazoRed, EstudioRed, MedidaRed } from "./tipos.ts";
import { chiCuadradoCola, normalDosColas } from "./funciones-especiales.ts";
import { cholesky, invertir, inversaSPD, matriz, type Mat } from "./numerica-compartida.ts";

const Z95 = 1.959963985;

// ---------------------------------------------------------------------------------------------
// Datos de cada estudio
// ---------------------------------------------------------------------------------------------

interface ContrasteBruto {
  /** Efecto del contraste brazo(a) − brazo(b) en la escala de análisis y su varianza propia. */
  y: number;
  v: number;
  /** Doble cero: sin eventos en los dos brazos, o con evento en todos los participantes de los dos (el OR no está definido). */
  dobleCero: boolean;
}

/**
 * Contraste entre dos brazos binarios: ln OR con su varianza 1/a + 1/b + 1/c + 1/d. La corrección de continuidad (Handbook 10.4.4.1) se
 * aplica a las cuatro celdas de TODOS los contrastes del estudio si alguna celda de algún brazo es cero (`incremento` = 0,5): es la regla
 * de `meta::pairwise` («only0» por estudio), que para un estudio de dos brazos coincide con corregir su único contraste y para uno de tres o
 * más corrige también los contrastes sin ceros. Un contraste con doble cero no informa.
 */
function contrasteBinario(ea: number, na: number, eb: number, nb: number, incremento: number): ContrasteBruto {
  const [a, b, c, d] = [ea, na - ea, eb, nb - eb].map((x) => x + incremento);
  return { y: Math.log(a) + Math.log(d) - Math.log(b) - Math.log(c), v: 1 / a + 1 / b + 1 / c + 1 / d, dobleCero: (ea === 0 && eb === 0) || (ea === na && eb === nb) };
}

interface BrazoPreparado {
  trat: string;
  /** Binario: eventos y total. Continuo: media y varianza de la media. */
  eventos?: number;
  total?: number;
  media?: number;
  varMedia?: number;
}

function brazoValido(brazo: BrazoRed, medida: MedidaRed): BrazoPreparado | null {
  const trat = brazo.tratamiento.trim();
  if (!trat) return null;
  if (medida === "OR") {
    const n = parseFloat(brazo.total ?? "");
    const e = parseFloat(brazo.eventos ?? "");
    return Number.isFinite(e) && Number.isFinite(n) && n > 0 && e >= 0 && e <= n ? { trat, eventos: e, total: n } : null;
  }
  const media = parseFloat(brazo.media ?? "");
  const de = parseFloat(brazo.de ?? "");
  const n = parseFloat(brazo.n ?? "");
  if (!Number.isFinite(media) || !Number.isFinite(de) || !Number.isFinite(n) || n <= 0 || de < 0) return null;
  const v = (de * de) / n;
  return Number.isFinite(v) && v > 0 ? { trat, media, varMedia: v } : null;
}

interface EstudioPreparado {
  etiqueta: string;
  /** Índices de tratamiento de los brazos válidos, en el orden capturado. */
  trat: number[];
  /** Efecto del contraste brazo(a) − brazo(b) (antisimétrica) y varianza propia del contraste (simétrica), entre los brazos del estudio. */
  Y: Mat;
  V: Mat;
}

interface Contraste {
  /** Posición del estudio en la lista de estudios preparados (identidad del estudio: dos estudios pueden compartir etiqueta). */
  estudio: number;
  i: number;
  j: number;
  /** Efecto del contraste: brazo(i) − brazo(j). */
  y: number;
  /** Peso de red (ya reducido si el estudio es multibrazo; incluye τ² si se pidió). */
  w: number;
  /** Varianza propia del contraste, sin τ² ni reducción multibrazo. */
  vPropia: number;
}

interface DatosRed {
  tratamientos: string[];
  estudios: EstudioPreparado[];
  descartados: { etiqueta: string; motivo: string }[];
}

function prepararDatos(estudiosRed: EstudioRed[], medida: MedidaRed): DatosRed {
  const nombres = new Set<string>();
  const crudos: { etiqueta: string; brazos: BrazoPreparado[]; Y: Mat; V: Mat }[] = [];
  const descartados: { etiqueta: string; motivo: string }[] = [];
  for (const e of estudiosRed.filter((x) => x.medida === medida)) {
    const etiqueta = e.etiqueta || `Estudio ${e.id}`;
    const brazos = e.brazos.map((b) => brazoValido(b, medida)).filter((b): b is BrazoPreparado => b !== null);
    if (brazos.length < 2) {
      descartados.push({ etiqueta, motivo: "menos de dos brazos con datos válidos" });
      continue;
    }
    if (new Set(brazos.map((b) => b.trat)).size !== brazos.length) {
      descartados.push({ etiqueta, motivo: "el mismo tratamiento aparece en más de un brazo" });
      continue;
    }
    const k = brazos.length;
    const Y = matriz(k, k), V = matriz(k, k);
    const incremento = medida === "OR" && brazos.some((b) => b.eventos === 0 || b.eventos === b.total) ? 0.5 : 0;
    let informativos = 0;
    for (let a = 0; a < k; a++) {
      for (let b = a + 1; b < k; b++) {
        if (medida === "OR") {
          const c = contrasteBinario(brazos[a].eventos!, brazos[a].total!, brazos[b].eventos!, brazos[b].total!, incremento);
          Y[a][b] = c.y; Y[b][a] = -c.y; V[a][b] = V[b][a] = c.v;
          if (!c.dobleCero) informativos++;
        } else {
          Y[a][b] = brazos[a].media! - brazos[b].media!; Y[b][a] = -Y[a][b];
          V[a][b] = V[b][a] = brazos[a].varMedia! + brazos[b].varMedia!;
          informativos++;
        }
      }
    }
    // Sin ningún contraste informativo el estudio no aporta nada al OR (Handbook 10.4.4.2). Con dos brazos, el único contraste es de doble
    // cero y se descarta. Con tres o más brazos, un contraste de doble cero se conserva con su corrección (equivale a `allstudies = TRUE` de
    // `pairwise`: `netmeta` no admite estudios multibrazo con contrastes descartados), porque los demás contrastes del estudio sí informan.
    if (informativos === 0) {
      descartados.push({ etiqueta, motivo: k === 2 ? "sin eventos (o con evento en todos) en los dos brazos: no informa el OR (Handbook 10.4.4.2)" : "sin contrastes informativos: sin eventos (o con evento en todos) en todos los brazos" });
      continue;
    }
    brazos.forEach((b) => nombres.add(b.trat));
    crudos.push({ etiqueta, brazos, Y, V });
  }
  const tratamientos = [...nombres].sort();
  const indice = new Map(tratamientos.map((t, i) => [t, i]));
  return {
    tratamientos,
    estudios: crudos.map((c) => ({ etiqueta: c.etiqueta, trat: c.brazos.map((b) => indice.get(b.trat)!), Y: c.Y, V: c.V })),
    descartados,
  };
}

// ---------------------------------------------------------------------------------------------
// Contrastes con los pesos de red (reducción de pesos de los estudios multibrazo)
// ---------------------------------------------------------------------------------------------

/**
 * Pesos de los k(k−1)/2 contrastes de un estudio a partir de las varianzas de sus contrastes R_ab (+ τ² en efectos aleatorios). Con el primer
 * brazo como referencia, la covarianza de los k − 1 contrastes es Σ_aa = R_0a y Σ_ab = (R_0a + R_0b − R_ab)/2 (identidad de polarización; con
 * varianzas aditivas, R_ab = vₐ + v_b, es Σ_ab = v₀ + δ_ab·vₐ); su inversa P es el Laplaciano del estudio con el nodo de referencia a tierra:
 * w_ab = −P_ab y w_0a = Σ_b P_ab. Con k = 2 es simplemente 1/R_01.
 */
function pesosDelEstudio(R: Mat, tau2: number): Mat | null {
  const k = R.length;
  const W = matriz(k, k);
  const r = (a: number, b: number) => (a === b ? 0 : R[a][b] + tau2);
  if (k === 2) {
    W[0][1] = W[1][0] = 1 / r(0, 1);
    return W;
  }
  const Sigma = matriz(k - 1, k - 1);
  for (let a = 1; a < k; a++) for (let b = 1; b < k; b++) Sigma[a - 1][b - 1] = a === b ? r(0, a) : (r(0, a) + r(0, b) - r(a, b)) / 2;
  const P = invertir(Sigma);
  if (!P) return null;
  for (let a = 0; a < k - 1; a++) {
    for (let b = 0; b < k - 1; b++) if (a !== b) W[a + 1][b + 1] = -P[a][b];
    const fila = P[a].reduce((s, x) => s + x, 0);
    W[0][a + 1] = W[a + 1][0] = fila;
  }
  return W;
}

function construirContrastes(datos: DatosRed, tau2: number): Contraste[] | null {
  const salida: Contraste[] = [];
  for (let s = 0; s < datos.estudios.length; s++) {
    const est = datos.estudios[s];
    const W = pesosDelEstudio(est.V, tau2);
    if (!W) return null;
    const k = est.trat.length;
    for (let a = 0; a < k; a++)
      for (let b = a + 1; b < k; b++)
        salida.push({ estudio: s, i: est.trat[a], j: est.trat[b], y: est.Y[a][b], w: W[a][b], vPropia: est.V[a][b] });
  }
  return salida;
}

// ---------------------------------------------------------------------------------------------
// Grafo: conectividad, Laplaciano y su pseudoinversa
// ---------------------------------------------------------------------------------------------

function esConexa(n: number, contrastes: Contraste[]): boolean {
  const padre = Array.from({ length: n }, (_, i) => i);
  const raiz = (x: number): number => (padre[x] === x ? x : (padre[x] = raiz(padre[x])));
  let componentes = n;
  for (const c of contrastes) {
    const a = raiz(c.i), b = raiz(c.j);
    if (a !== b) { padre[a] = b; componentes--; }
  }
  return componentes === 1;
}

/**
 * Pseudoinversa de Moore-Penrose del Laplaciano de un grafo conexo. Se invierte el Laplaciano reducido (sin el último nodo, definido
 * positivo) y se centra: L⁺ = C·G·C con C = I − 11ᵀ/n y G la inversa reducida rellena con ceros.
 */
function pseudoinversaLaplaciana(L: Mat): Mat | null {
  const n = L.length;
  const reducida = L.slice(0, n - 1).map((f) => f.slice(0, n - 1));
  // Con pesos negativos (estudios multibrazo con varianzas de contraste no aditivas) el Laplaciano reducido puede no ser definido positivo
  const inv = inversaSPD(reducida)?.inv ?? invertir(reducida);
  if (!inv) return null;
  const G = matriz(n, n);
  for (let i = 0; i < n - 1; i++) for (let j = 0; j < n - 1; j++) G[i][j] = inv[i][j];
  const mediaFila = G.map((f) => f.reduce((s, x) => s + x, 0) / n);
  const mediaTotal = mediaFila.reduce((s, x) => s + x, 0) / n;
  return G.map((f, i) => f.map((x, j) => x - mediaFila[i] - mediaFila[j] + mediaTotal));
}

interface Ajuste {
  n: number;
  contrastes: Contraste[];
  Lmas: Mat;
  theta: number[];
  q: number;
}

function ajustarRed(n: number, contrastes: Contraste[]): Ajuste | null {
  const L = matriz(n, n);
  const flujo = new Array<number>(n).fill(0);
  for (const c of contrastes) {
    L[c.i][c.i] += c.w; L[c.j][c.j] += c.w;
    L[c.i][c.j] -= c.w; L[c.j][c.i] -= c.w;
    flujo[c.i] += c.w * c.y; flujo[c.j] -= c.w * c.y;
  }
  const Lmas = pseudoinversaLaplaciana(L);
  if (!Lmas) return null;
  const theta = Lmas.map((f) => f.reduce((s, x, j) => s + x * flujo[j], 0));
  const q = contrastes.reduce((s, c) => s + c.w * (c.y - (theta[c.i] - theta[c.j])) ** 2, 0);
  return { n, contrastes, Lmas, theta, q };
}

/**
 * Traza de (I − H)·U·W, con H_ef = w_f·(L⁺_{i_e,i_f} − L⁺_{i_e,j_f} − L⁺_{j_e,i_f} + L⁺_{j_e,j_f}) y U_ef = ½·(B_e·B_f) entre
 * contrastes del mismo estudio (B_e·B_f = [i_e=i_f] − [i_e=j_f] − [j_e=i_f] + [j_e=j_f]; 1 en la diagonal). Es el denominador de τ².
 */
function trazaDeHeterogeneidad(aj: Ajuste): number {
  const { contrastes: C, Lmas } = aj;
  let traza = 0;
  for (let e = 0; e < C.length; e++) {
    for (let f = 0; f < C.length; f++) {
      if (C[e].estudio !== C[f].estudio) continue;
      const u = 0.5 * (+(C[e].i === C[f].i) - +(C[e].i === C[f].j) - +(C[e].j === C[f].i) + +(C[e].j === C[f].j));
      if (u === 0) continue;
      const h = C[f].w * (Lmas[C[e].i][C[f].i] - Lmas[C[e].i][C[f].j] - Lmas[C[e].j][C[f].i] + Lmas[C[e].j][C[f].j]);
      traza += ((e === f ? 1 : 0) - h) * u * C[e].w;
    }
  }
  return traza;
}

// ---------------------------------------------------------------------------------------------
// Descomposición de Q en heterogeneidad dentro de diseños e inconsistencia entre diseños
// ---------------------------------------------------------------------------------------------

export interface DesglosePorDiseno {
  tratamientos: string[];
  nEstudios: number;
  q: number;
  gl: number;
}

export interface DecomposicionQ {
  qDentro: number;
  glDentro: number;
  pDentro: number;
  qEntre: number;
  glEntre: number;
  pEntre: number;
  /** Solo diseños con ≥2 estudios (los únicos que aportan a Q-dentro). */
  disenos: DesglosePorDiseno[];
}

/** Resuelve A x = b (A simétrica definida positiva) por Cholesky. */
function resolverSPD(A: Mat, b: number[]): number[] | null {
  const L = cholesky(A);
  if (!L) return null;
  const n = b.length;
  const z = new Array<number>(n).fill(0);
  for (let i = 0; i < n; i++) z[i] = (b[i] - L[i].slice(0, i).reduce((s, x, k) => s + x * z[k], 0)) / L[i][i];
  const x = new Array<number>(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = z[i];
    for (let k = i + 1; k < n; k++) s -= L[k][i] * x[k];
    x[i] = s / L[i][i];
  }
  return x;
}

/**
 * Q dentro de cada diseño (Krahn, Binder y König 2013, ecs. 3, 4 y 10): los contrastes de cada estudio contra el primer tratamiento del
 * diseño (en orden alfabético) se combinan por mínimos cuadrados generalizados con la covarianza propia de cada estudio, y Q_d es la suma
 * de las distancias de Mahalanobis de cada estudio a ese promedio. Q entre diseños = Q_total − Q dentro (ec. 8).
 */
function descomponerQ(datos: DatosRed, qTotal: number, glTotal: number): DecomposicionQ {
  const grupos = new Map<string, { nombres: string[]; filas: { y: number[]; S: Mat }[] }>();
  for (const est of datos.estudios) {
    const orden = est.trat.map((tr, a) => ({ tr, a })).sort((x, z) => (datos.tratamientos[x.tr] < datos.tratamientos[z.tr] ? -1 : 1));
    const nombres = orden.map((o) => datos.tratamientos[o.tr]);
    const ref = orden[0].a;
    const resto = orden.slice(1).map((o) => o.a);
    const y = resto.map((a) => est.Y[a][ref]);
    const p = y.length;
    // Σ de los contrastes contra el primer tratamiento del diseño: Σ_aa = R_ref,a y Σ_ab = (R_ref,a + R_ref,b − R_ab)/2
    const S = matriz(p, p);
    for (let a = 0; a < p; a++) for (let b = 0; b < p; b++) S[a][b] = a === b ? est.V[ref][resto[a]] : (est.V[ref][resto[a]] + est.V[ref][resto[b]] - est.V[resto[a]][resto[b]]) / 2;
    const clave = nombres.join("|");
    const g = grupos.get(clave) ?? { nombres, filas: [] };
    g.filas.push({ y, S });
    grupos.set(clave, g);
  }

  let qDentro = 0, glDentro = 0;
  const disenos: DesglosePorDiseno[] = [];
  for (const g of grupos.values()) {
    const r = g.filas.length;
    if (r < 2) continue; // un solo estudio del diseño: no hay heterogeneidad que medir
    const p = g.filas[0].y.length;
    // Precisión acumulada y su producto con los contrastes: θ̂ = (Σ Sₛ⁻¹)⁻¹ Σ Sₛ⁻¹ yₛ
    const precisiones = g.filas.map((f) => inversaSPD(f.S));
    if (precisiones.some((x) => !x)) continue; // covarianza singular (datos degenerados): se omite este diseño
    const suma = matriz(p, p);
    const sumaPorY = new Array<number>(p).fill(0);
    g.filas.forEach((f, s) => {
      const P = precisiones[s]!.inv;
      for (let a = 0; a < p; a++) {
        for (let b = 0; b < p; b++) suma[a][b] += P[a][b];
        sumaPorY[a] += P[a].reduce((acc, x, b) => acc + x * f.y[b], 0);
      }
    });
    const media = resolverSPD(suma, sumaPorY);
    if (!media) continue;
    let qDiseno = 0;
    g.filas.forEach((f, s) => {
      const d = f.y.map((x, a) => x - media[a]);
      const P = precisiones[s]!.inv;
      for (let a = 0; a < p; a++) for (let b = 0; b < p; b++) qDiseno += d[a] * P[a][b] * d[b];
    });
    const glDiseno = (r - 1) * p;
    qDentro += qDiseno;
    glDentro += glDiseno;
    disenos.push({ tratamientos: g.nombres, nEstudios: r, q: qDiseno, gl: glDiseno });
  }

  const qEntre = Math.max(0, qTotal - qDentro);
  const glEntre = glTotal - glDentro;
  return {
    qDentro,
    glDentro,
    pDentro: glDentro > 0 ? chiCuadradoCola(qDentro, glDentro) : NaN,
    qEntre,
    glEntre,
    pEntre: glEntre > 0 ? chiCuadradoCola(qEntre, glEntre) : NaN,
    disenos,
  };
}

// ---------------------------------------------------------------------------------------------
// Resultado
// ---------------------------------------------------------------------------------------------

export interface EstimacionRed {
  tratamiento1: string;
  tratamiento2: string;
  /** Estimación en escala natural: OR para binario, diferencia de medias para continuo. */
  estimacion: number;
  icInferior: number;
  icSuperior: number;
  /** Diferencia y su SE en la escala de análisis (log para OR, lineal para DM) — tratamiento1−tratamiento2.
   *  Las usan `calcularPScores` y el node-splitting; no hace falta destransformarlas para comparar. */
  diff: number;
  se: number;
  /** true si al menos un estudio comparó estos dos tratamientos directamente. */
  esDirecta: boolean;
  nEstudiosDirectos: number;
}

/**
 * Node-splitting (inconsistencia local): separa, para una comparación con evidencia directa, cuánto del estimado de red viene de esa
 * evidencia directa y cuánto de la indirecta (inferida del resto de la red), y prueba si discrepan. Solo se calcula para pares con
 * evidencia directa Y con al menos una ruta indirecta (P < 1).
 */
export interface NodeSplitEstimacion {
  tratamiento1: string;
  tratamiento2: string;
  nEstudiosDirectos: number;
  /** Proporción de la información total de la red que viene de la evidencia directa (0–1). Solo la define el método de
   *  back-calculation; SIDDE no la usa (queda en NaN). */
  proporcionDirecta: number;
  estimacionDirecta: number;
  icInferiorDirecta: number;
  icSuperiorDirecta: number;
  /** Directo, en la escala de análisis (log/lineal) — para que `calcularNodeSplittingSIDDE` pueda reutilizarlo sin recalcular. */
  diffDirecto: number;
  seDirecto: number;
  estimacionIndirecta: number;
  icInferiorIndirecta: number;
  icSuperiorIndirecta: number;
  /** Prueba z de la diferencia directa vs. indirecta — p bajo sugiere inconsistencia en esta comparación. */
  z: number;
  pValor: number;
}

export interface ResultadoNMA {
  tratamientos: string[];
  medida: MedidaRed;
  escala: "log" | "lineal";
  k: number;
  m: number;
  n: number;
  modelo: "fijo" | "aleatorio";
  /** Todas las comparaciones i<j (n(n−1)/2), para la tabla de liga y el gráfico de red. */
  estimaciones: EstimacionRed[];
  /** Solo las comparaciones con evidencia directa Y ruta indirecta disponible. */
  nodeSplitting: NodeSplitEstimacion[];
  q: number;
  gl: number;
  pValor: number;
  i2: number;
  tau2: number;
  tau: number;
  /** Q_total descompuesto en heterogeneidad-dentro-de-diseño e inconsistencia-entre-diseños. `null` si no hay gl que repartir. */
  decomposicionQ: DecomposicionQ | null;
  /** Estudios que no entraron a la red, con el motivo (menos de dos brazos válidos, o el mismo tratamiento en dos brazos). */
  descartados: { etiqueta: string; motivo: string }[];
}

/** Función de distribución de la normal estándar, a partir de la cola de dos lados exacta. */
const normalAcumulada = (x: number): number => {
  const dosColas = normalDosColas(x);
  return x >= 0 ? 1 - 0.5 * dosColas : 0.5 * dosColas;
};

const aEscalaNatural = (x: number, escala: "log" | "lineal"): number => (escala === "log" ? Math.exp(x) : x);

/**
 * Calcula el meta-análisis en red para la medida indicada.
 * @returns `null` si hay < 3 tratamientos, < 2 estudios utilizables, o la red no está conexa.
 */
export function calcularNMA(
  estudiosRed: EstudioRed[],
  medida: MedidaRed,
  modelo: "fijo" | "aleatorio" = "aleatorio",
): ResultadoNMA | null {
  const datos = prepararDatos(estudiosRed, medida);
  const n = datos.tratamientos.length;
  if (n < 3) return null;

  const contrastesFijos = construirContrastes(datos, 0);
  if (!contrastesFijos || !esConexa(n, contrastesFijos)) return null;
  const fijo = ajustarRed(n, contrastesFijos);
  if (!fijo) return null;

  const k = datos.estudios.length;
  const glTotal = datos.estudios.reduce((s, e) => s + (e.trat.length - 1), 0) - (n - 1);
  const q = fijo.q;

  // τ² por momentos
  let tau2 = 0;
  if (glTotal > 0) {
    const traza = trazaDeHeterogeneidad(fijo);
    tau2 = traza > 0 ? Math.max(0, (q - glTotal) / traza) : 0;
  }
  const i2 = q > 0 ? Math.max(0, ((q - glTotal) / q) * 100) : 0;
  const pValor = glTotal > 0 ? chiCuadradoCola(q, glTotal) : NaN;

  let final = fijo;
  let contrastesFinales = contrastesFijos;
  if (modelo === "aleatorio" && tau2 > 0) {
    const cr = construirContrastes(datos, tau2);
    const aj = cr ? ajustarRed(n, cr) : null;
    if (!cr || !aj) return null;
    final = aj;
    contrastesFinales = cr;
  }

  const escala: "log" | "lineal" = medida === "OR" ? "log" : "lineal";
  const directos = new Map<number, number>(); // clave i·n + j (i < j) → nº de contrastes directos
  const clavePar = (a: number, b: number) => (a < b ? a * n + b : b * n + a);
  for (const c of contrastesFijos) directos.set(clavePar(c.i, c.j), (directos.get(clavePar(c.i, c.j)) ?? 0) + 1);

  const estimaciones: EstimacionRed[] = [];
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const diff = final.theta[i] - final.theta[j];
      const se = Math.sqrt(Math.max(0, final.Lmas[i][i] + final.Lmas[j][j] - 2 * final.Lmas[i][j]));
      const directo = directos.get(clavePar(i, j)) ?? 0;
      estimaciones.push({
        tratamiento1: datos.tratamientos[i],
        tratamiento2: datos.tratamientos[j],
        estimacion: aEscalaNatural(diff, escala),
        icInferior: aEscalaNatural(diff - Z95 * se, escala),
        icSuperior: aEscalaNatural(diff + Z95 * se, escala),
        diff,
        se,
        esDirecta: directo > 0,
        nEstudiosDirectos: directo,
      });
    }
  }

  // Node-splitting por back-calculation. La evidencia directa de un par es el promedio por inverso-varianza de SUS contrastes con la
  // varianza propia de cada uno (sin la reducción multibrazo, que solo existe para la red completa) más τ² en efectos aleatorios.
  // Convención comprobada contra la salida de `netsplit` (con estudios multibrazo la reducción de pesos daría otro directo).
  const nodeSplitting: NodeSplitEstimacion[] = [];
  const tau2Directo = modelo === "aleatorio" ? tau2 : 0;
  for (const est of estimaciones) {
    if (!est.esDirecta) continue;
    const i = datos.tratamientos.indexOf(est.tratamiento1);
    const j = datos.tratamientos.indexOf(est.tratamiento2);
    let pesoDirecto = 0, sumaPesada = 0;
    for (const c of contrastesFinales) {
      const mismo = c.i === i && c.j === j;
      const inverso = c.i === j && c.j === i;
      if (!mismo && !inverso) continue;
      const w = 1 / (c.vPropia + tau2Directo);
      pesoDirecto += w;
      sumaPesada += w * (mismo ? c.y : -c.y);
    }
    if (pesoDirecto === 0) continue;
    const directo = sumaPesada / pesoDirecto;
    const varDirecta = 1 / pesoDirecto;
    const varRed = est.se * est.se;
    const P = pesoDirecto * varRed; // proporción de evidencia directa = Var(red) / Var(directa)
    if (P >= 1 - 1e-9) continue; // toda la información es directa: no hay ruta indirecta que separar
    const indirecto = (est.diff - P * directo) / (1 - P);
    const varIndirecta = varRed / (1 - P);
    const seDirecto = Math.sqrt(varDirecta), seIndirecto = Math.sqrt(varIndirecta);
    const z = (directo - indirecto) / Math.sqrt(varDirecta + varIndirecta);
    nodeSplitting.push({
      tratamiento1: est.tratamiento1,
      tratamiento2: est.tratamiento2,
      nEstudiosDirectos: est.nEstudiosDirectos,
      proporcionDirecta: P,
      estimacionDirecta: aEscalaNatural(directo, escala),
      icInferiorDirecta: aEscalaNatural(directo - Z95 * seDirecto, escala),
      icSuperiorDirecta: aEscalaNatural(directo + Z95 * seDirecto, escala),
      diffDirecto: directo,
      seDirecto,
      estimacionIndirecta: aEscalaNatural(indirecto, escala),
      icInferiorIndirecta: aEscalaNatural(indirecto - Z95 * seIndirecto, escala),
      icSuperiorIndirecta: aEscalaNatural(indirecto + Z95 * seIndirecto, escala),
      z,
      pValor: chiCuadradoCola(z * z, 1),
    });
  }

  return {
    tratamientos: datos.tratamientos,
    medida,
    escala,
    k,
    m: contrastesFijos.length,
    n,
    modelo,
    estimaciones,
    nodeSplitting,
    q,
    gl: glTotal,
    pValor,
    i2,
    tau2,
    tau: Math.sqrt(tau2),
    decomposicionQ: glTotal > 0 ? descomponerQ(datos, q, glTotal) : null,
    descartados: datos.descartados,
  };
}

export interface PScoreTratamiento {
  tratamiento: string;
  pScore: number;
}

/**
 * Ranking por P-score (Rücker y Schwarzer 2015): para cada par, la probabilidad de que i sea mejor que j es Φ(−z) si un valor menor es
 * mejor (Φ(z) si es mayor), con z = diferencia/EE en la escala de análisis; el P-score de i es su promedio contra los demás tratamientos
 * (1 = seguramente el mejor; 0 = seguramente el peor). Es la versión frecuentista de SUCRA, sin remuestreo.
 *
 * @param menorEsMejor - true si un valor más bajo del desenlace es mejor (mortalidad, síntomas); false si uno más alto es mejor.
 */
export function calcularPScores(resultado: ResultadoNMA, menorEsMejor = true): PScoreTratamiento[] {
  const { tratamientos, estimaciones } = resultado;
  const n = tratamientos.length;
  return tratamientos.map((tratamiento, i) => {
    let suma = 0;
    for (let j = 0; j < n; j++) {
      if (j === i) continue;
      const par = estimaciones.find(
        (e) =>
          (e.tratamiento1 === tratamiento && e.tratamiento2 === tratamientos[j]) ||
          (e.tratamiento1 === tratamientos[j] && e.tratamiento2 === tratamiento),
      )!;
      const z = (par.tratamiento1 === tratamiento ? par.diff : -par.diff) / par.se;
      suma += normalAcumulada(menorEsMejor ? -z : z);
    }
    return { tratamiento, pScore: suma / (n - 1) };
  });
}

/**
 * Node-splitting por SIDDE (Efthimiou y cols. 2019): para cada comparación con estimado directo, se quitan los estudios completos que
 * comparan directamente el par (un estudio multibrazo que lo incluya también sale) y se reajusta la red con el resto; la estimación de red
 * de ese reajuste del par es el estimado indirecto. El estimado directo es el del back-calculation, que no depende del método.
 *
 * @returns Solo las comparaciones donde el reajuste sin evidencia directa sigue siendo una red válida (≥ 3 tratamientos, conexa, y los
 *   dos tratamientos del par siguen presentes).
 */
export function calcularNodeSplittingSIDDE(
  estudiosRed: EstudioRed[],
  medida: MedidaRed,
  modelo: "fijo" | "aleatorio" = "aleatorio",
): NodeSplitEstimacion[] {
  const completo = calcularNMA(estudiosRed, medida, modelo);
  if (!completo) return [];
  const resultado: NodeSplitEstimacion[] = [];
  for (const base of completo.nodeSplitting) {
    const sinDirecta = estudiosRed.filter((e) => {
      if (e.medida !== medida) return true;
      const tratados = new Set(e.brazos.map((b) => b.tratamiento.trim()));
      return !(tratados.has(base.tratamiento1) && tratados.has(base.tratamiento2));
    });
    const red = calcularNMA(sinDirecta, medida, modelo);
    if (!red || !red.tratamientos.includes(base.tratamiento1) || !red.tratamientos.includes(base.tratamiento2)) continue;
    const par = red.estimaciones.find(
      (e) =>
        (e.tratamiento1 === base.tratamiento1 && e.tratamiento2 === base.tratamiento2) ||
        (e.tratamiento1 === base.tratamiento2 && e.tratamiento2 === base.tratamiento1),
    );
    if (!par) continue;
    const indirecto = par.tratamiento1 === base.tratamiento1 ? par.diff : -par.diff;
    const seIndirecto = par.se;
    const z = (base.diffDirecto - indirecto) / Math.sqrt(base.seDirecto ** 2 + seIndirecto ** 2);
    resultado.push({
      ...base,
      proporcionDirecta: NaN,
      estimacionIndirecta: aEscalaNatural(indirecto, completo.escala),
      icInferiorIndirecta: aEscalaNatural(indirecto - Z95 * seIndirecto, completo.escala),
      icSuperiorIndirecta: aEscalaNatural(indirecto + Z95 * seIndirecto, completo.escala),
      z,
      pValor: chiCuadradoCola(z * z, 1),
    });
  }
  return resultado;
}
