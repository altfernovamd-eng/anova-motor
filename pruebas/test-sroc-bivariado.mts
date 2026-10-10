// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Pruebas del modelo bivariado de precisión diagnóstica (sroc-bivariado.ts) y de su integración en calcularSROC.
 * La comparación numérica con R está en validacion/comparar-bivariado.mts; aquí van las propiedades que no dependen de R:
 * exactitud de la cuadratura, simetrías, casos límite, geometría de la curva y de las regiones, honestidad de los avisos.
 *   npm run test:sroc-bivariado
 */
import { existsSync, readFileSync } from "node:fs";
import {
  ajustarBivariado,
  ajustarBivariadoConCache,
  areaBajoCurva,
  curvaSROC,
  elipseROC,
  logVerosimilitudBinomial,
  nodosGaussHermite,
  regionesROC,
  sensibilidadEnCurva,
  tablasUtilizables,
  type ParametrosBivariado,
  type TablaDiagnostica,
} from "../src/sroc-bivariado.ts";
import { calcularSROC } from "../src/statistics.ts";
import type { Estudio } from "../src/tipos.ts";

let fallos = 0;
const ok = (c: boolean, m: string) => { console.log(`${c ? "✓" : "✗"} ${m}`); if (!c) fallos++; };
const cerca = (a: number, b: number, tol: number) => Math.abs(a - b) <= tol;
const expit = (x: number) => 1 / (1 + Math.exp(-x));
const logit = (p: number) => Math.log(p / (1 - p));

const T = (vp: number, fp: number, fn: number, vn: number): TablaDiagnostica => ({ vp, fp, fn, vn });
// Conjunto con heterogeneidad y correlación negativa (efecto de umbral), tamaños moderados; el óptimo es interior (ρ ≈ −0,91)
const datos: TablaDiagnostica[] = [
  T(45, 5, 5, 45), T(40, 8, 10, 42), T(90, 3, 12, 95), T(30, 4, 6, 60), T(70, 9, 8, 80),
  T(55, 14, 5, 70), T(25, 2, 11, 52), T(62, 8, 4, 40), T(33, 3, 9, 71), T(48, 11, 6, 55),
];
const P0: ParametrosBivariado = { muSens: 1.6, muEsp: 1.9, sdSens: 0.7, sdEsp: 0.5, rho: -0.4 };

// ---------- Cuadratura de Gauss-Hermite ----------
for (const q of [5, 9, 15, 21]) {
  const n = nodosGaussHermite(q);
  const w = n.logw.map(Math.exp);
  const m = (k: number) => n.z.reduce((s, z, i) => s + w[i] * z ** k, 0);
  ok(cerca(w.reduce((a, b) => a + b, 0), 1, 1e-12), `GH ${q}: los pesos suman 1`);
  ok(cerca(m(2), 1, 1e-10) && cerca(m(4), 3, 1e-9) && (q < 5 || cerca(m(6), 15, 1e-8)), `GH ${q}: momentos E[z²]=1, E[z⁴]=3${q >= 5 ? ", E[z⁶]=15" : ""}`);
  ok(cerca(m(1), 0, 1e-12) && cerca(m(3), 0, 1e-12), `GH ${q}: simetría (momentos impares nulos)`);
}

// ---------- Log-verosimilitud ----------
{
  // Con varianza entre estudios ≈ 0 la verosimilitud marginal es la binomial simple en σ(μ)
  const p: ParametrosBivariado = { muSens: 1.2, muEsp: 0.8, sdSens: 1e-3, sdEsp: 1e-3, rho: 0 };
  const lnC = (n: number, k: number) => { let s = 0; for (let i = 1; i <= k; i++) s += Math.log((n - k + i) / i); return s; };
  let directo = 0;
  for (const t of datos) directo += lnC(t.vp + t.fn, t.vp) + t.vp * Math.log(expit(p.muSens)) + t.fn * Math.log(1 - expit(p.muSens)) + lnC(t.vn + t.fp, t.vn) + t.vn * Math.log(expit(p.muEsp)) + t.fp * Math.log(1 - expit(p.muEsp));
  ok(cerca(logVerosimilitudBinomial(datos, p, 21), directo, 0.05), `σ→0: la verosimilitud marginal tiende a la binomial simple (${logVerosimilitudBinomial(datos, p, 21).toFixed(4)} vs ${directo.toFixed(4)})`);
  ok(cerca(logVerosimilitudBinomial(datos, P0, 21), logVerosimilitudBinomial(datos, P0, 41), 1e-7), "21 y 41 nodos dan la misma log-verosimilitud (cuadratura convergida)");
  const inv = datos.map((t) => T(t.vn, t.fn, t.fp, t.vp)); // intercambia los papeles: enfermos ↔ sanos
  const Pinv: ParametrosBivariado = { muSens: P0.muEsp, muEsp: P0.muSens, sdSens: P0.sdEsp, sdEsp: P0.sdSens, rho: P0.rho };
  ok(cerca(logVerosimilitudBinomial(datos, P0, 21), logVerosimilitudBinomial(inv, Pinv, 21), 1e-9), "simetría: intercambiar sensibilidad y especificidad deja la verosimilitud igual");
  const rev = [...datos].reverse();
  ok(cerca(logVerosimilitudBinomial(datos, P0, 15), logVerosimilitudBinomial(rev, P0, 15), 1e-10), "el orden de los estudios no importa");
  const conRho = (r: number) => logVerosimilitudBinomial(datos, { ...P0, rho: r }, 21);
  ok(Number.isFinite(conRho(0.999)) && Number.isFinite(conRho(-0.999)) && conRho(-0.4) !== conRho(0.4), "ρ muy cercano a ±1 no da NaN y el signo de ρ importa");
}

// ---------- Ajuste ----------
const bin = ajustarBivariado(datos, { metodo: "binomial" })!;
const nor = ajustarBivariado(datos, { metodo: "normal" })!;
ok(!!bin && bin.metodo === "binomial" && bin.estimacion === "ML" && bin.k === 10, "binomial: método ML con k=10");
ok(!!nor && nor.metodo === "normal" && nor.estimacion === "REML", "normal: método REML");
ok(bin.sensibilidad.icInferior < bin.sensibilidad.estimacion && bin.sensibilidad.estimacion < bin.sensibilidad.icSuperior && bin.especificidad.icInferior < bin.especificidad.estimacion && bin.especificidad.estimacion < bin.especificidad.icSuperior, "los IC contienen la estimación");
ok(cerca(bin.fpr.estimacion, 1 - bin.especificidad.estimacion, 1e-12) && cerca(bin.fpr.icInferior, 1 - bin.especificidad.icSuperior, 1e-12), "FPR = 1 − especificidad (IC invertido)");
ok(bin.vcovMu[0][0] > 0 && bin.vcovMu[1][1] > 0 && bin.vcovMu[0][0] * bin.vcovMu[1][1] > bin.vcovMu[0][1] ** 2, "la covarianza de las medias es definida positiva");
ok(Math.abs(bin.parametros.rho) < 1 && bin.parametros.sdSens > 0 && bin.parametros.sdEsp > 0, "parámetros dentro del rango");
// el óptimo no es mejorable: ninguna perturbación aumenta la verosimilitud
{
  const base = bin.logVerosimilitud;
  const pert = (d: Partial<ParametrosBivariado>): number => logVerosimilitudBinomial(datos, { ...bin.parametros, ...d }, 21);
  const q = bin.parametros;
  const mejores = [pert({ muSens: q.muSens + 0.02 }), pert({ muSens: q.muSens - 0.02 }), pert({ muEsp: q.muEsp + 0.02 }), pert({ muEsp: q.muEsp - 0.02 }),
    pert({ sdSens: q.sdSens * 1.05 }), pert({ sdSens: q.sdSens / 1.05 }), pert({ sdEsp: q.sdEsp * 1.05 }), pert({ sdEsp: q.sdEsp / 1.05 }),
    pert({ rho: Math.tanh(Math.atanh(q.rho) + 0.05) }), pert({ rho: Math.tanh(Math.atanh(q.rho) - 0.05) })];
  ok(mejores.every((v) => v <= base + 1e-8), "es un máximo local de la verosimilitud (10 perturbaciones, ninguna mejora)");
  ok(cerca(base, logVerosimilitudBinomial(datos, bin.parametros, 21), 1e-9), "logVerosimilitud informada = la recalculada en el óptimo");
}
// simetría del ajuste
{
  const inv = datos.map((t) => T(t.vn, t.fn, t.fp, t.vp));
  const b2 = ajustarBivariado(inv, { metodo: "binomial" })!;
  ok(cerca(b2.parametros.muSens, bin.parametros.muEsp, 2e-3) && cerca(b2.parametros.muEsp, bin.parametros.muSens, 2e-3), "ajuste simétrico: las medias se intercambian");
  ok(cerca(b2.parametros.sdSens, bin.parametros.sdEsp, 5e-3) && cerca(b2.parametros.rho, bin.parametros.rho, 1e-2), "ajuste simétrico: σ se intercambian y ρ se conserva");
  ok(cerca(b2.logVerosimilitud, bin.logVerosimilitud, 1e-6), "ajuste simétrico: misma log-verosimilitud");
  const rev = ajustarBivariado([...datos].reverse(), { metodo: "binomial" })!;
  ok(cerca(rev.parametros.muSens, bin.parametros.muSens, 1e-3) && cerca(rev.logVerosimilitud, bin.logVerosimilitud, 1e-6), "el orden de los estudios no cambia el ajuste");
}
// con estudios grandes, normal y binomial dan casi lo mismo
{
  const grandes = datos.map((t) => T(t.vp * 20, t.fp * 20, t.fn * 20, t.vn * 20));
  const b = ajustarBivariado(grandes, { metodo: "binomial" })!;
  const n = ajustarBivariado(grandes, { metodo: "normal" })!;
  ok(cerca(b.sensibilidad.estimacion, n.sensibilidad.estimacion, 5e-3) && cerca(b.especificidad.estimacion, n.especificidad.estimacion, 5e-3), `estudios grandes: binomial y normal coinciden (sens ${b.sensibilidad.estimacion.toFixed(3)} vs ${n.sensibilidad.estimacion.toFixed(3)})`);
  ok(cerca(b.parametros.sdSens, n.parametros.sdSens, 0.05) && cerca(b.parametros.sdEsp, n.parametros.sdEsp, 0.05), "estudios grandes: σ coinciden");
}
// pocos estudios y datos no utilizables
ok(ajustarBivariado(datos.slice(0, 3)) === null, "3 estudios → null (no inventa)");
ok(ajustarBivariado(datos.slice(0, 4)) !== null && ajustarBivariado(datos.slice(0, 4))!.advertencias.some((a) => /pocos|solo 4/i.test(a)), "4 estudios → ajusta y avisa de la poca precisión");
ok(tablasUtilizables([T(1, 2, 3, 4), T(0, 0, 0, 5), T(3, 0, 0, 0), T(1.5, 1, 1, 1), T(-1, 1, 1, 1), T(NaN, 1, 1, 1)]).length === 1, "tablas inutilizables (sin enfermos o sin sanos, no enteras, negativas, NaN) se descartan");
ok(ajustarBivariado([...datos.slice(0, 3), T(0, 0, 0, 20), T(0, 0, 0, 9)]) === null, "con 3 estudios utilizables tras descartar → null");
// celdas en cero
{
  const conCeros = [T(0, 3, 10, 50), T(12, 0, 2, 60), T(8, 4, 0, 40), T(20, 6, 5, 0 + 80), T(15, 1, 3, 70), T(9, 2, 4, 30)];
  const z = ajustarBivariado(conCeros, { metodo: "binomial" })!;
  const zn = ajustarBivariado(conCeros, { metodo: "normal" })!;
  ok(!!z && [z.parametros.muSens, z.parametros.muEsp, z.logVerosimilitud, z.vcovMu[0][0]].every(Number.isFinite), "celdas en cero: el modelo binomial ajusta sin NaN");
  ok(!!zn && [zn.parametros.muSens, zn.parametros.muEsp, zn.vcovMu[0][0]].every(Number.isFinite), "celdas en cero: la aproximación normal (con +0,5) también");
}
// sin heterogeneidad: varianza en el límite, declarada
{
  const h = [T(50, 5, 10, 95), T(50, 5, 10, 95), T(100, 10, 20, 190), T(25, 2, 5, 48), T(50, 5, 10, 95), T(75, 7, 15, 143)];
  const b = ajustarBivariado(h, { metodo: "binomial" })!;
  ok(b.limites.includes("sdSens") || b.limites.includes("sdEsp"), `sin heterogeneidad real: la varianza queda en el límite y se declara (${b.limites.join(", ") || "ninguno"})`);
  ok(b.advertencias.some((a) => /límite|nula|prácticamente/i.test(a)), "…con una advertencia legible");
  ok(Number.isFinite(b.sensibilidad.icInferior) && b.sensibilidad.icSuperior < 1 && b.sensibilidad.icInferior > 0, "…y aun así el IC de la sensibilidad es finito y propio");
}

// ---------- Curva, AUC y regiones ----------
{
  const p = bin.parametros;
  ok(cerca(sensibilidadEnCurva(p, expit(-p.muEsp)), expit(p.muSens), 1e-12), "la curva SROC pasa por el punto resumen");
  ok(cerca(sensibilidadEnCurva(P0, 0.05), 0.534392, 1e-6) && cerca(sensibilidadEnCurva(P0, 0.3), 0.955798, 1e-6), "curva SROC: valores calculados a mano con la pendiente σ_sens/σ_esp (P0, FPR 0,05 y 0,3)");
  const c = curvaSROC(p, 0.02, 0.6, 60);
  ok(c.length === 60 && cerca(c[0].fpr, 0.02, 1e-9) && cerca(c[59].fpr, 0.6, 1e-9), "curva: 60 puntos entre los FPR pedidos");
  ok(c.every((q, i) => i === 0 || (q.fpr > c[i - 1].fpr && q.tpr >= c[i - 1].tpr)), "la curva crece con el FPR (σ>0)");
  ok(cerca(areaBajoCurva({ muSens: 0, muEsp: 0, sdSens: 1, sdEsp: 1, rho: 0 }), 0.5, 1e-8), "AUC = 0,5 cuando la prueba no discrimina (TPR = FPR)");
  const aucSimple = (d: number) => { const n = 400000, a = -40, b = 40, h = (b - a) / n; let s = 0; for (let i = 0; i <= n; i++) { const u = a + i * h, f = expit(u); s += (i === 0 || i === n ? 1 : i % 2 ? 4 : 2) * expit(u + d) * f * (1 - f); } return (s * h) / 3; };
  ok(cerca(areaBajoCurva({ muSens: 2, muEsp: 0, sdSens: 1, sdEsp: 1, rho: 0 }), aucSimple(2), 1e-8), `AUC con pendiente 1 = integral independiente (Simpson): ${areaBajoCurva({ muSens: 2, muEsp: 0, sdSens: 1, sdEsp: 1, rho: 0 }).toFixed(6)}`);
  ok(areaBajoCurva({ ...P0, muSens: 3 }) > areaBajoCurva(P0) && areaBajoCurva(P0) > 0.5 && areaBajoCurva(P0) < 1, "AUC aumenta con la sensibilidad media y está en (0,5; 1)");
  ok(bin.auc > 0.5 && bin.auc < 1 && cerca(bin.auc, areaBajoCurva(bin.parametros), 1e-12), `AUC del resultado = la de la curva (${bin.auc.toFixed(4)})`);
  const [fMin, fMax] = bin.rangoFprObservado;
  ok(fMin > 0 && fMax < 1 && fMin < fMax && fMin >= Math.min(...datos.map((t) => t.fp / (t.fp + t.vn))) - 1e-12, `rango de FPR observado (${fMin.toFixed(3)}–${fMax.toFixed(3)}) contenido en los datos`);
}
{
  // elipse: todos los puntos cumplen d'V⁻¹d = χ²₂(0,95); crece con el nivel
  const V: [[number, number], [number, number]] = [[0.09, -0.03], [-0.03, 0.05]];
  const det = V[0][0] * V[1][1] - V[0][1] ** 2;
  const forma = (q: { fpr: number; tpr: number }, c: [number, number]) => {
    const d1 = logit(q.tpr) - c[0], d2 = -logit(q.fpr) - c[1];
    return (V[1][1] * d1 * d1 - 2 * V[0][1] * d1 * d2 + V[0][0] * d2 * d2) / det;
  };
  const e95 = elipseROC([1.5, 2], V, 0.95, 200);
  ok(e95.length === 200 && e95.every((q) => cerca(forma(q, [1.5, 2]), -2 * Math.log(0.05), 1e-9)), "elipse 95%: todos los puntos cumplen la forma cuadrática χ²₂(0,95) = 5,9915");
  const e80 = elipseROC([1.5, 2], V, 0.8, 200);
  ok(e80.every((q) => cerca(forma(q, [1.5, 2]), -2 * Math.log(0.2), 1e-9)), "elipse 80%: χ²₂(0,80)");
  const ancho = (e: { fpr: number }[]) => Math.max(...e.map((q) => q.fpr)) - Math.min(...e.map((q) => q.fpr));
  ok(ancho(e95) > ancho(e80), "la elipse de mayor nivel es más grande");
  const R = regionesROC(bin, 0.95);
  ok(R.confianza.length > 50 && R.prediccion.length > 50, "regiones: confianza y predicción calculadas");
  const anchoC = ancho(R.confianza), anchoP = ancho(R.prediccion);
  ok(anchoP > anchoC, `la región de predicción (${anchoP.toFixed(3)}) es mayor que la de confianza (${anchoC.toFixed(3)})`);
  // cada punto de la región de confianza cae DENTRO de la de predicción (misma forma, V más grande)
  const centro: [number, number] = [bin.parametros.muSens, bin.parametros.muEsp];
  const s11 = bin.parametros.sdSens ** 2, s22 = bin.parametros.sdEsp ** 2, s12 = bin.parametros.rho * bin.parametros.sdSens * bin.parametros.sdEsp;
  const Vp = [[bin.vcovMu[0][0] + s11, bin.vcovMu[0][1] + s12], [bin.vcovMu[1][0] + s12, bin.vcovMu[1][1] + s22]];
  const dp = Vp[0][0] * Vp[1][1] - Vp[0][1] ** 2;
  const dentro = R.confianza.every((q) => { const d1 = logit(q.tpr) - centro[0], d2 = -logit(q.fpr) - centro[1]; return (Vp[1][1] * d1 * d1 - 2 * Vp[0][1] * d1 * d2 + Vp[0][0] * d2 * d2) / dp < -2 * Math.log(0.05); });
  ok(dentro, "la región de confianza está contenida en la de predicción");
}

// ---------- Cache ----------
{
  const a = ajustarBivariadoConCache(datos, "binomial");
  const b = ajustarBivariadoConCache(datos, "binomial");
  const c = ajustarBivariadoConCache(datos, "normal");
  ok(a === b && a !== c && a?.metodo === "binomial" && c?.metodo === "normal", "caché: mismo objeto para las mismas tablas y método; otro método, otro ajuste");
  const d = ajustarBivariadoConCache(datos.map((t, i) => (i === 0 ? T(t.vp + 1, t.fp, t.fn, t.vn) : t)), "binomial");
  ok(d !== a, "caché: cambiar un dato invalida el resultado");
}

// ---------- Integración en calcularSROC ----------
const est = (id: number, t: TablaDiagnostica): Estudio =>
  ({ id, autor: `D${id}`, año: String(2010 + id), diseño: "Transversal", n: String(t.vp + t.fp + t.fn + t.vn), medida: "OR", valor: "", icInferior: "", icSuperior: "", datosDiagnosticos: { vp: String(t.vp), fp: String(t.fp), fn: String(t.fn), vn: String(t.vn) } }) as Estudio;
{
  const estudios = datos.map((t, i) => est(i + 1, t));
  const s = calcularSROC(estudios, "aleatorio", "DL")!;
  ok(!!s.bivariado && s.bivariado.metodo === "binomial" && !!s.bivariadoNormal && s.bivariadoNormal.metodo === "normal", "calcularSROC: devuelve el bivariado binomial y la comparación normal");
  ok(cerca(s.bivariado!.parametros.muSens, bin.parametros.muSens, 1e-6), "calcularSROC usa el mismo ajuste que el módulo");
  ok(s.regresion !== null && s.dor.gl === 9, "calcularSROC conserva Moses-Littenberg y el DOR (exploratorios)");
  const sin = calcularSROC(estudios, "aleatorio", "DL", { bivariado: false })!;
  ok(sin.bivariado === null && sin.bivariadoNormal === null && sin.puntos.length === 10, "{bivariado:false} no ajusta (para conteos rápidos)");
  const tres = calcularSROC(estudios.slice(0, 3), "aleatorio", "DL")!;
  ok(tres.bivariado === null && tres.puntos.length === 3, "3 estudios: sin bivariado");
  const sinDatos = estudios.map((e, i) => (i === 4 ? { ...e, datosDiagnosticos: { vp: "", fp: "", fn: "", vn: "" } } : e));
  ok(calcularSROC(sinDatos, "aleatorio", "DL")!.puntos.length === 9, "un estudio sin 2×2 completa no entra");
}

// ---------- Regresión con datos reales (mada::AuditC) ----------
const rutaRef = new URL("../validacion/referencias/referencia-bivariado.json", import.meta.url);
if (existsSync(rutaRef)) {
  const ref = JSON.parse(readFileSync(rutaRef, "utf8"));
  const au = ref.conjuntos.find((c: { id: string }) => c.id === "auditc");
  const tab = au.datos.map((d: { vp: number; fp: number; fn: number; vn: number }) => T(d.vp, d.fp, d.fn, d.vn));
  const b = ajustarBivariado(tab, { metodo: "binomial" })!;
  const rb = au.binomial;
  ok(cerca(b.logVerosimilitud, rb.loglik, 2e-5), `AuditC (k=${tab.length}): log-verosimilitud ${b.logVerosimilitud.toFixed(6)} = R ${rb.loglik.toFixed(6)}`);
  ok(cerca(b.parametros.muSens, rb.par[0], 2e-3) && cerca(b.parametros.muEsp, rb.par[1], 2e-3), "AuditC: μ iguales a la implementación independiente de R");
  ok(b.sensibilidad.estimacion > 0.9 && b.sensibilidad.estimacion < 0.95 && b.especificidad.estimacion > 0.75 && b.especificidad.estimacion < 0.8, "AuditC: sensibilidad ≈ 92 % y especificidad ≈ 78 %");
} else console.log("— (referencia de R no encontrada: se omite la regresión con AuditC)");

if (fallos) {
  console.error(`\n${fallos} prueba(s) fallaron.`);
  process.exit(1);
}
console.log("\nTodas las pruebas pasaron.");
