// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Validación diferencial del modelo bivariado de precisión diagnóstica (`sroc-bivariado.ts`) contra R.
 * La invoca `comparar.mts` con sus ayudantes. Referencia: `referencias/referencia-bivariado.json`
 * (generada por `generar-referencia-bivariado.R`: 30 conjuntos simulados + `mada::AuditC`).
 *
 * Qué se contrasta y contra qué (cada pieza contra la fuente que corresponde, con su tolerancia):
 *  1. Aproximación normal (REML): contra el óptimo EXACTO de la misma función objetivo (R propio, `optim` estricto) y, con tolerancia
 *     más holgada, contra `mada::reitsma` (mvmeta se detiene antes de converger del todo).
 *  2. Modelo binomial exacto (ML): contra una implementación INDEPENDIENTE en R (Gauss-Hermite adaptativa con `eigen`, `optim`, `numDeriv`):
 *     log-verosimilitud, parámetros y covarianza de las medias; contra `lme4::glmer` (Laplace) solo como comprobación de que el modelo
 *     está definido igual (tolerancia de aproximación).
 *  3. La cuadratura: la log-verosimilitud contra fuerza bruta (`integrate` anidado).
 *  4. Curva SROC, AUC y elipses, evaluadas con los parámetros de `mada` para aislarlas del ajuste: contra `sroc()`, `AUC()` y `ROCellipse()`.
 */
import { readFileSync } from "node:fs";
import {
  ajustarBivariado,
  areaBajoCurva,
  elipseROC,
  logVerosimilitudBinomial,
  sensibilidadEnCurva,
  type ParametrosBivariado,
  type TablaDiagnostica,
} from "../src/sroc-bivariado.ts";
import type { AyudasComparar } from "./comparar-avanzada.mts";

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type J = any;

const CHI2_95 = -2 * Math.log(1 - 0.95);
const logit = (p: number) => Math.log(p / (1 - p));
const expit = (x: number) => 1 / (1 + Math.exp(-x));

export function compararBivariado(h: AyudasComparar) {
  const { cmp, sinPar, arr } = h;
  const ref: J = JSON.parse(readFileSync(new URL("./referencias/referencia-bivariado.json", import.meta.url), "utf8"));

  for (const e of ref.conjuntos as J[]) {
    const tablas: TablaDiagnostica[] = (e.datos as J[]).map((d) => ({ vp: d.vp, fp: d.fp, fn: d.fn, vn: d.vn }));
    const id = `${e.id} (k=${e.k}, ${e.tipo})`;

    // ---------- 1) Aproximación normal (REML) ----------
    const nor = ajustarBivariado(tablas, { metodo: "normal" });
    if (!nor) sinPar.push(`bivariado normal ${id}: aNova no ajusta`);
    else {
      const ex = e.reml_exacto;
      const sdA = Math.exp(ex.par[0]), sdB = Math.exp(ex.par[1]), rho = Math.tanh(ex.par[2]);
      const enLimite = sdA < 2e-3 || sdB < 2e-3 || Math.abs(ex.par[2]) > 6.8;
      if (!enLimite) {
        cmp("SROC bivariado · normal REML vs óptimo exacto de R · σ", id, nor.parametros.sdSens, sdA, 1e-4, 1e-5);
        cmp("SROC bivariado · normal REML vs óptimo exacto de R · σ", id, nor.parametros.sdEsp, sdB, 1e-4, 1e-5);
        cmp("SROC bivariado · normal REML vs óptimo exacto de R · ρ", id, nor.parametros.rho, rho, 1e-3, 1e-4);
      } else {
        cmp("SROC bivariado · normal REML en el límite (σ o ρ) · mismo límite", id, nor.limites.length > 0 ? 1 : 0, 1, 1e-12, 1e-12);
      }
      if (e.mada) {
        const m = e.mada;
        const coef = arr(m.coef), vc = arr(m.vcov), psi = arr(m.Psi);
        // mada trabaja con (logit sens, logit FPR): la especificidad es el signo contrario y las covarianzas cruzadas cambian de signo.
        cmp("SROC bivariado · normal REML vs mada · μ sensibilidad", id, nor.parametros.muSens, coef[0], 2e-3, 2e-4);
        cmp("SROC bivariado · normal REML vs mada · μ especificidad", id, nor.parametros.muEsp, -coef[1], 2e-3, 2e-4);
        cmp("SROC bivariado · normal REML vs mada · var(μ)", id, nor.vcovMu[0][0], vc[0], 1e-2, 1e-4);
        cmp("SROC bivariado · normal REML vs mada · var(μ)", id, nor.vcovMu[1][1], vc[3], 1e-2, 1e-4);
        cmp("SROC bivariado · normal REML vs mada · cov(μ)", id, nor.vcovMu[0][1], -vc[1], 1e-2, 5e-4);
        if (!enLimite) {
          // mvmeta (mada) se detiene antes de converger: con una varianza casi nula (< 0,05) su σ puede quedar a ~1e-2 del óptimo exacto
          // (que aNova sí alcanza: ver la comparación con el óptimo exacto). Fuera de esa zona, la tolerancia es estrecha.
          const tolSd = Math.min(sdA, sdB) < 0.05 ? 2e-2 : 5e-3;
          cmp("SROC bivariado · normal REML vs mada · σ", id, nor.parametros.sdSens, Math.sqrt(psi[0]), 5e-2, tolSd);
          cmp("SROC bivariado · normal REML vs mada · σ", id, nor.parametros.sdEsp, Math.sqrt(psi[3]), 5e-2, tolSd);
        }
        const sens = arr(m.sens), fpr = arr(m.fpr);
        cmp("SROC bivariado · normal REML vs mada · sensibilidad e IC", id, nor.sensibilidad.estimacion, sens[0], 1e-3, 1e-4);
        cmp("SROC bivariado · normal REML vs mada · sensibilidad e IC", id, nor.sensibilidad.icInferior, sens[1], 5e-3, 5e-4);
        cmp("SROC bivariado · normal REML vs mada · sensibilidad e IC", id, nor.sensibilidad.icSuperior, sens[2], 5e-3, 5e-4);
        cmp("SROC bivariado · normal REML vs mada · FPR e IC", id, nor.fpr.estimacion, fpr[0], 1e-3, 1e-4);
        cmp("SROC bivariado · normal REML vs mada · FPR e IC", id, nor.fpr.icInferior, fpr[1], 5e-3, 5e-4);
        cmp("SROC bivariado · normal REML vs mada · FPR e IC", id, nor.fpr.icSuperior, fpr[2], 5e-3, 5e-4);
      }
    }

    // ---------- 2) Modelo binomial exacto (ML) ----------
    const bin = ajustarBivariado(tablas, { metodo: "binomial" });
    const rb = e.binomial;
    if (!bin || bin.metodo !== "binomial") sinPar.push(`bivariado binomial ${id}: aNova ${bin ? "recurre a la aproximación normal" : "no ajusta"}`);
    else {
      cmp("SROC bivariado · binomial exacto vs R independiente · logLik", id, bin.logVerosimilitud, rb.loglik, 1e-9, 2e-5);
      const par = arr(rb.par);
      cmp("SROC bivariado · binomial exacto vs R independiente · μ", id, bin.parametros.muSens, par[0], 1e-3, 1e-3);
      cmp("SROC bivariado · binomial exacto vs R independiente · μ", id, bin.parametros.muEsp, par[1], 1e-3, 1e-3);
      const cotaSd = rb.limite;
      if (!cotaSd) {
        cmp("SROC bivariado · binomial exacto vs R independiente · σ", id, bin.parametros.sdSens, Math.exp(par[2]), 5e-3, 2e-3);
        cmp("SROC bivariado · binomial exacto vs R independiente · σ", id, bin.parametros.sdEsp, Math.exp(par[3]), 5e-3, 2e-3);
        cmp("SROC bivariado · binomial exacto vs R independiente · ρ", id, bin.parametros.rho, Math.tanh(par[4]), 1e-2, 5e-3);
      } else {
        cmp("SROC bivariado · binomial exacto en el límite (σ o ρ) · mismo límite", id, bin.limites.length > 0 ? 1 : 0, 1, 1e-12, 1e-12);
      }
      const vc = rb.vcovMu ? arr(rb.vcovMu) : null;
      if (vc) {
        cmp("SROC bivariado · binomial exacto vs R independiente · var(μ)", id, bin.vcovMu[0][0], vc[0], 2e-2, 1e-4);
        cmp("SROC bivariado · binomial exacto vs R independiente · var(μ)", id, bin.vcovMu[1][1], vc[3], 2e-2, 1e-4);
        cmp("SROC bivariado · binomial exacto vs R independiente · cov(μ)", id, bin.vcovMu[0][1], vc[1], 5e-2, 1e-4);
      }
      // Sensibilidad/especificidad resumen e IC (Wald en logit) a partir de μ y vcov de R (identidad algebraica)
      cmp("SROC bivariado · binomial exacto · sensibilidad resumen", id, bin.sensibilidad.estimacion, expit(par[0]), 1e-3, 5e-4);
      cmp("SROC bivariado · binomial exacto · especificidad resumen", id, bin.especificidad.estimacion, expit(par[1]), 1e-3, 5e-4);
      // lme4::glmer (Laplace): comprobación del modelo, no de la aproximación numérica
      if (e.glmer && !cotaSd) {
        const g = e.glmer;
        const fx = arr(g.fixef);
        cmp("SROC bivariado · binomial exacto vs glmer (Laplace, solo comprobación) · μ", id, bin.parametros.muSens, fx[0], 0.1, 0.1);
        cmp("SROC bivariado · binomial exacto vs glmer (Laplace, solo comprobación) · μ", id, bin.parametros.muEsp, fx[1], 0.1, 0.1);
        const sd = arr(g.sd);
        cmp("SROC bivariado · binomial exacto vs glmer (Laplace, solo comprobación) · σ", id, bin.parametros.sdSens, sd[0], 0.2, 0.1);
        cmp("SROC bivariado · binomial exacto vs glmer (Laplace, solo comprobación) · σ", id, bin.parametros.sdEsp, sd[1], 0.2, 0.1);
      }
    }

    // ---------- 3) Log-verosimilitud contra la implementación independiente de R en el óptimo (misma cuadratura, otro código) ----------
    if (rb) {
      const p: ParametrosBivariado = { muSens: arr(rb.par)[0], muEsp: arr(rb.par)[1], sdSens: Math.exp(arr(rb.par)[2]), sdEsp: Math.exp(arr(rb.par)[3]), rho: Math.tanh(arr(rb.par)[4]) };
      cmp("SROC bivariado · log-verosimilitud en el óptimo de R (GH 21) vs R independiente", id, logVerosimilitudBinomial(tablas, p, 21), rb.loglik, 1e-11, 1e-7);
    }

    // ---------- 4) Curva, AUC y elipse con los parámetros de mada ----------
    if (e.mada) {
      const m = e.mada;
      const coef = arr(m.coef), psi = arr(m.Psi), vc = arr(m.vcov);
      const sA = Math.sqrt(psi[0]), sB = Math.sqrt(psi[3]);
      const rho = -psi[1] / (sA * sB);
      const p: ParametrosBivariado = { muSens: coef[0], muEsp: -coef[1], sdSens: sA, sdEsp: sB, rho };
      const fprs = arr(m.fprs), rg = arr(m.sroc_rg);
      fprs.forEach((f, i) => cmp("SROC bivariado · curva SROC vs mada::sroc (Rutter-Gatsonis)", `${e.id} FPR=${f}`, sensibilidadEnCurva(p, f), rg[i], 1e-7, 1e-8));
      cmp("SROC bivariado · AUC vs mada::AUC (rejilla fina)", id, areaBajoCurva(p), m.auc_fino, 5e-4, 5e-4);
      // Elipse de confianza: los puntos de mada deben cumplir d'V⁻¹d = χ²₂(0,95) con la V de mada, y los de aNova también;
      // además la extensión (mín./máx. de FPR y sensibilidad) de ambas curvas debe coincidir.
      const V = [[vc[0], -vc[1]], [-vc[2], vc[3]]]; // (logit sens, logit esp)
      const det = V[0][0] * V[1][1] - V[0][1] * V[1][0];
      const forma = (fpr: number, tpr: number) => {
        const d1 = logit(tpr) - p.muSens, d2 = -logit(fpr) - p.muEsp;
        return (V[1][1] * d1 * d1 - 2 * V[0][1] * d1 * d2 + V[0][0] * d2 * d2) / det;
      };
      const mel: number[][] = (m.elipse as J[]).map((r: J) => arr(r));
      // mada devuelve una matriz (n×2): en JSON es un arreglo de filas [fpr, sens]
      for (const [f, s] of mel) cmp("SROC bivariado · elipse de mada cumple la forma cuadrática (χ²₂ 95 %)", id, forma(f, s), CHI2_95, 1e-6, 1e-6);
      const anova = elipseROC([p.muSens, p.muEsp], V as [[number, number], [number, number]], 0.95, 400);
      for (const q of anova) cmp("SROC bivariado · elipse de aNova cumple la forma cuadrática (χ²₂ 95 %)", id, forma(q.fpr, q.tpr), CHI2_95, 1e-9, 1e-9);
      const ext = (xs: number[]) => [Math.min(...xs), Math.max(...xs)];
      const [fa0, fa1] = ext(anova.map((q) => q.fpr)), [ta0, ta1] = ext(anova.map((q) => q.tpr));
      const [fm0, fm1] = ext(mel.map((r) => r[0])), [tm0, tm1] = ext(mel.map((r) => r[1]));
      cmp("SROC bivariado · extensión de la elipse vs mada", `${id} FPR mín.`, fa0, fm0, 1e-2, 2e-3);
      cmp("SROC bivariado · extensión de la elipse vs mada", `${id} FPR máx.`, fa1, fm1, 1e-2, 2e-3);
      cmp("SROC bivariado · extensión de la elipse vs mada", `${id} sens. mín.`, ta0, tm0, 1e-2, 2e-3);
      cmp("SROC bivariado · extensión de la elipse vs mada", `${id} sens. máx.`, ta1, tm1, 1e-2, 2e-3);
    }
  }
  // ---------- Cuadratura contra FUERZA BRUTA (integrate anidado por tramos) en un conjunto pequeño y 5 puntos del espacio de parámetros ----------
  const tb = (ref.bruta.datos as J[]).map((d) => ({ vp: d.vp, fp: d.fp, fn: d.fn, vn: d.vn }));
  for (const [i, b] of (ref.bruta.puntos as J[]).entries()) {
    const par = arr(b.par);
    const p: ParametrosBivariado = { muSens: par[0], muEsp: par[1], sdSens: Math.exp(par[2]), sdEsp: Math.exp(par[3]), rho: Math.tanh(par[4]) };
    cmp("SROC bivariado · log-verosimilitud (GH 21) vs fuerza bruta (integrate)", `punto ${i + 1}`, logVerosimilitudBinomial(tb, p, 21), b.loglik_bruta, 1e-9, 5e-8);
    cmp("SROC bivariado · log-verosimilitud (GH 21) vs R GH 21 independiente", `punto ${i + 1}`, logVerosimilitudBinomial(tb, p, 21), b.loglik_gh21, 1e-11, 1e-8);
  }
  return { conjuntos: (ref.conjuntos as J[]).length, mada: ref.paquetes.mada as string, lme4: ref.paquetes.lme4 as string };
}
