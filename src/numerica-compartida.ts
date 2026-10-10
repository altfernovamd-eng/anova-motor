// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Piezas numéricas compartidas por los modelos que integran efectos aleatorios o necesitan optimizar sin derivadas
 * (modelo bivariado de precisión diagnóstica, modelos mixtos logísticos, modelo multinivel): cuadratura de Gauss-Hermite
 * (Golub y Welsch 1969), Nelder-Mead (1965), Hessiano numérico por diferencias centrales e inversa por Gauss-Jordan.
 * Lógica pura, sin dependencias de la app salvo `funciones-especiales`.
 */
import { lnGamma } from "./funciones-especiales.ts";

export const sigmoide = (x: number) => (x >= 0 ? 1 / (1 + Math.exp(-x)) : Math.exp(x) / (1 + Math.exp(x)));
export const logit = (p: number) => Math.log(p / (1 - p));
/** log σ(x) sin desbordes. */
export const logSigmoide = (x: number) => (x >= 0 ? -Math.log1p(Math.exp(-x)) : x - Math.log1p(Math.exp(x)));
export const limitar = (x: number, lo: number, hi: number) => (x < lo ? lo : x > hi ? hi : x);
export const lnChoose = (n: number, k: number) => lnGamma(n + 1) - lnGamma(k + 1) - lnGamma(n - k + 1);

/** Autovalores y autovectores de una matriz simétrica por el método de Jacobi cíclico. Devuelve los autovalores y los vectores en columnas. */
export function jacobiSimetrica(a: number[][]): { valores: number[]; vectores: number[][] } {
  const n = a.length;
  const A = a.map((f) => f.slice());
  const V: number[][] = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
  for (let barrido = 0; barrido < 100; barrido++) {
    let fuera = 0;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) fuera += A[i][j] * A[i][j];
    if (fuera < 1e-30) break;
    for (let p = 0; p < n - 1; p++) {
      for (let q = p + 1; q < n; q++) {
        if (Math.abs(A[p][q]) < 1e-300) continue;
        const theta = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1);
        const s = t * c;
        for (let k = 0; k < n; k++) {
          const akp = A[k][p];
          const akq = A[k][q];
          A[k][p] = c * akp - s * akq;
          A[k][q] = s * akp + c * akq;
        }
        for (let k = 0; k < n; k++) {
          const apk = A[p][k];
          const aqk = A[q][k];
          A[p][k] = c * apk - s * aqk;
          A[q][k] = s * apk + c * aqk;
        }
        for (let k = 0; k < n; k++) {
          const vkp = V[k][p];
          const vkq = V[k][q];
          V[k][p] = c * vkp - s * vkq;
          V[k][q] = s * vkp + c * vkq;
        }
      }
    }
  }
  return { valores: A.map((f, i) => f[i]), vectores: V };
}

export interface NodosGH {
  /** Nodos de la cuadratura para la normal estándar, z. */
  z: number[];
  /** Logaritmos de los pesos normalizados (suman 1). */
  logw: number[];
}
const cacheGH = new Map<number, NodosGH>();

/**
 * Nodos y pesos de Gauss-Hermite para la normal estándar con `q` puntos (Golub-Welsch: valores propios de la matriz de Jacobi de los
 * polinomios de Hermite). Integran exactamente los polinomios hasta grado 2q − 1 contra la densidad normal estándar.
 */
export function nodosGaussHermite(q: number): NodosGH {
  // Con muchos nodos y un integrando estrecho, la cuadratura ADAPTATIVA se vuelve inestable (los nodos de las colas multiplican un
  // cociente integrando/normal que crece): se comprobó que 21 a 61 nodos coinciden a 1e-14 y que desde ~100 el resultado se corrompe.
  if (!Number.isInteger(q) || q < 2 || q > 64) throw new RangeError(`Número de nodos de Gauss-Hermite no admitido (${q}): use entre 2 y 64.`);
  const guardado = cacheGH.get(q);
  if (guardado) return guardado;
  const J = Array.from({ length: q }, () => new Array<number>(q).fill(0));
  for (let i = 1; i < q; i++) {
    J[i][i - 1] = J[i - 1][i] = Math.sqrt(i); // para la normal estándar (polinomios de Hermite probabilistas): √i
  }
  const { valores, vectores } = jacobiSimetrica(J);
  const orden = valores.map((v, i) => [v, i] as const).sort((a, b) => a[0] - b[0]);
  const z = orden.map(([v]) => v);
  // El peso del nodo j es el cuadrado de la primera componente de su autovector (la medida de la normal estándar tiene masa 1).
  const w = orden.map(([, i]) => vectores[0][i] ** 2);
  const suma = w.reduce((s, x) => s + x, 0);
  const res = { z, logw: w.map((x) => Math.log(x / suma)) };
  cacheGH.set(q, res);
  return res;
}

/** Minimiza `f` por Nelder-Mead con reinicios; devuelve el mejor punto y su valor. */
export function nelderMead(
  f: (x: number[]) => number,
  x0: number[],
  opc: { paso?: number; tolX?: number; tolF?: number; maxEval?: number; reinicios?: number } = {},
): { x: number[]; fx: number; evaluaciones: number; convergio: boolean } {
  const n = x0.length;
  const paso = opc.paso ?? 0.5;
  const tolX = opc.tolX ?? 1e-10;
  const tolF = opc.tolF ?? 1e-12;
  const maxEval = opc.maxEval ?? 4000;
  const reinicios = opc.reinicios ?? 4;
  let mejor = x0.slice();
  let fMejor = f(mejor);
  let evals = 1;
  let convergio = false;
  for (let r = 0; r <= reinicios; r++) {
    const simplex: number[][] = [mejor.slice()];
    for (let i = 0; i < n; i++) {
      const v = mejor.slice();
      v[i] += paso * (r === 0 ? 1 : 0.2 / (r + 1)) * (v[i] === 0 ? 1 : 1);
      simplex.push(v);
    }
    let fs = simplex.map((v) => (evals++, f(v)));
    let it = 0;
    while (evals < maxEval) {
      const idx = fs.map((v, i) => [v, i] as const).sort((a, b) => a[0] - b[0]).map(([, i]) => i);
      simplex.splice(0, simplex.length, ...idx.map((i) => simplex[i]));
      fs = idx.map((i) => fs[i]);
      const tamano = Math.max(...simplex.slice(1).map((v) => Math.max(...v.map((x, j) => Math.abs(x - simplex[0][j])))));
      if (tamano < tolX && Math.abs(fs[n] - fs[0]) < tolF) {
        convergio = true;
        break;
      }
      const centro = simplex[0].map((_, j) => simplex.slice(0, n).reduce((s, v) => s + v[j], 0) / n);
      const reflejo = centro.map((c, j) => c + (c - simplex[n][j]));
      const fr = (evals++, f(reflejo));
      if (fr < fs[0]) {
        const expansion = centro.map((c, j) => c + 2 * (c - simplex[n][j]));
        const fe = (evals++, f(expansion));
        if (fe < fr) {
          simplex[n] = expansion;
          fs[n] = fe;
        } else {
          simplex[n] = reflejo;
          fs[n] = fr;
        }
      } else if (fr < fs[n - 1]) {
        simplex[n] = reflejo;
        fs[n] = fr;
      } else {
        const fuera = fr < fs[n];
        const contraccion = centro.map((c, j) => c + 0.5 * ((fuera ? reflejo[j] : simplex[n][j]) - c));
        const fc = (evals++, f(contraccion));
        if (fc < (fuera ? fr : fs[n])) {
          simplex[n] = contraccion;
          fs[n] = fc;
        } else {
          for (let i = 1; i <= n; i++) {
            simplex[i] = simplex[i].map((x, j) => simplex[0][j] + 0.5 * (x - simplex[0][j]));
            fs[i] = (evals++, f(simplex[i]));
          }
        }
      }
      it++;
      if (it > 100000) break;
    }
    const iMin = fs.indexOf(Math.min(...fs));
    const mejoro = fs[iMin] < fMejor - 1e-15;
    if (fs[iMin] <= fMejor) {
      mejor = simplex[iMin].slice();
      fMejor = fs[iMin];
    }
    if (!mejoro && r > 0 && convergio) break;
    if (evals >= maxEval) break;
  }
  return { x: mejor, fx: fMejor, evaluaciones: evals, convergio };
}

/** Hessiano numérico (diferencias centrales) de `f` en `x`, con paso relativo `h`. */
export function hessianoNumerico(f: (x: number[]) => number, x: number[], h = 1e-3): number[][] {
  const n = x.length;
  const H = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  const f0 = f(x);
  const pasos = x.map((v) => h * Math.max(1, Math.abs(v)));
  for (let i = 0; i < n; i++) {
    const xp = x.slice();
    const xm = x.slice();
    xp[i] += pasos[i];
    xm[i] -= pasos[i];
    H[i][i] = (f(xp) - 2 * f0 + f(xm)) / (pasos[i] * pasos[i]);
    for (let j = i + 1; j < n; j++) {
      const pp = x.slice();
      const pm = x.slice();
      const mp = x.slice();
      const mm = x.slice();
      pp[i] += pasos[i]; pp[j] += pasos[j];
      pm[i] += pasos[i]; pm[j] -= pasos[j];
      mp[i] -= pasos[i]; mp[j] += pasos[j];
      mm[i] -= pasos[i]; mm[j] -= pasos[j];
      H[i][j] = H[j][i] = (f(pp) - f(pm) - f(mp) + f(mm)) / (4 * pasos[i] * pasos[j]);
    }
  }
  return H;
}

/** Inversa de una matriz cuadrada por Gauss-Jordan con pivote parcial; `null` si es singular. */
export function invertir(a: number[][]): number[][] | null {
  const n = a.length;
  const M = a.map((f, i) => [...f, ...Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))]);
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    if (Math.abs(M[p][c]) < 1e-300) return null;
    [M[c], M[p]] = [M[p], M[c]];
    const piv = M[c][c];
    for (let j = 0; j < 2 * n; j++) M[c][j] /= piv;
    for (let r = 0; r < n; r++) {
      if (r === c) continue;
      const fct = M[r][c];
      if (fct === 0) continue;
      for (let j = 0; j < 2 * n; j++) M[r][j] -= fct * M[c][j];
    }
  }
  return M.map((f) => f.slice(n));
}


// ---------------------------------------------------------------------------------------------
// Álgebra lineal mínima para matrices pequeñas (denso, sin dependencias)
// ---------------------------------------------------------------------------------------------


export type Mat = number[][];

export const matriz = (f: number, c: number, v = 0): Mat => Array.from({ length: f }, () => new Array<number>(c).fill(v));
export const mult = (A: Mat, B: Mat): Mat => {
  const f = A.length, c = B[0].length, n = B.length;
  const R = matriz(f, c);
  for (let i = 0; i < f; i++) for (let k = 0; k < n; k++) { const a = A[i][k]; if (a !== 0) for (let j = 0; j < c; j++) R[i][j] += a * B[k][j]; }
  return R;
};
export const trans = (A: Mat): Mat => (A.length === 0 ? [] : A[0].map((_, j) => A.map((fila) => fila[j])));
export const multVec = (A: Mat, v: number[]): number[] => A.map((fila) => fila.reduce((s, a, j) => s + a * v[j], 0));

/** Factor de Cholesky inferior L (A = L L'); `null` si A no es definida positiva. */
export function cholesky(A: Mat): Mat | null {
  const n = A.length;
  const L = matriz(n, n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let s = A[i][j];
      for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      if (i === j) {
        if (!(s > 0)) return null;
        L[i][i] = Math.sqrt(s);
      } else L[i][j] = s / L[j][j];
    }
  }
  return L;
}

/** Inversa y log-determinante de una matriz simétrica definida positiva por Cholesky. */
export function inversaSPD(A: Mat): { inv: Mat; logDet: number } | null {
  const L = cholesky(A);
  if (!L) return null;
  const n = A.length;
  const Li = matriz(n, n);
  for (let j = 0; j < n; j++) {
    for (let i = j; i < n; i++) {
      let s = i === j ? 1 : 0;
      for (let k = j; k < i; k++) s -= L[i][k] * Li[k][j];
      Li[i][j] = s / L[i][i];
    }
  }
  const inv = mult(trans(Li), Li);
  let logDet = 0;
  for (let i = 0; i < n; i++) logDet += 2 * Math.log(L[i][i]);
  return { inv, logDet };
}

/** Raíz cuadrada inversa simétrica (con inversa generalizada: los valores propios casi nulos se tratan como 0). */
export function raizInversaSimetrica(B: Mat): Mat {
  const { valores, vectores } = jacobiSimetrica(B);
  const tope = Math.max(...valores.map(Math.abs), 1e-300);
  const n = B.length;
  const R = matriz(n, n);
  for (let k = 0; k < n; k++) {
    const lam = valores[k];
    if (!(lam > 1e-12 * tope)) continue;
    const f = 1 / Math.sqrt(lam);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) R[i][j] += f * vectores[i][k] * vectores[j][k];
  }
  return R;
}

