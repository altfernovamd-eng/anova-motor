// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Validación diferencial de los modelos mixtos logísticos (`glmm.ts`) contra R.
 * La invoca `comparar.mts` con sus ayudantes. Referencia: `referencias/referencia-glmm.json`
 * (generada por `generar-referencia-glmm.R`: 40 conjuntos de proporciones y hasta 40 de tablas 2×2).
 *
 * Qué se contrasta:
 *  1. Contra una implementación INDEPENDIENTE en R (Gauss-Hermite adaptativa propia, `optim`, `numDeriv`): log-verosimilitud, μ, τ y el
 *     error estándar de μ (inversa del Hessiano completo).
 *  2. Contra `metafor::rma.glmm` como caja negra (PLO para proporciones; CM.EL para OR, solo donde converge): μ, error estándar y τ².
 *     El error estándar de metafor (vía lme4) está condicionado a τ, el de aNova no: se compara con tolerancia de aproximación.
 *  3. La verosimilitud marginal contra FUERZA BRUTA (`integrate` por tramos) en un conjunto pequeño de cada familia; las probabilidades
 *     de la hipergeométrica no central de la referencia salen de BiasedUrn::dFNCHypergeo.
 */
import { readFileSync } from "node:fs";
import { ajustarGLMMOdds, ajustarGLMMProporcion, logVerosimilitudOdds, logVerosimilitudProporcion } from "../src/glmm.ts";
import type { AyudasComparar } from "./comparar-avanzada.mts";

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type J = any;

export function compararGLMM(h: AyudasComparar) {
  const { cmp, sinPar, arr } = h;
  const ref: J = JSON.parse(readFileSync(new URL("./referencias/referencia-glmm.json", import.meta.url), "utf8"));

  for (const e of ref.conjuntos as J[]) {
    const esProp = e.id.startsWith("prop");
    const filas = e.datos as J[];
    const r = esProp
      ? ajustarGLMMProporcion(filas.map((d) => ({ x: d.x, n: d.n })))
      : ajustarGLMMOdds(filas.map((d) => ({ eventosExp: d.a, totalExp: d.n1, eventosCtrl: d.c, totalCtrl: d.n2 })));
    const id = `${e.id} (k=${filas.length})`;
    const fam = esProp ? "proporciones (binomial-normal)" : "OR (hipergeométrica-normal)";
    if (!r) { sinPar.push(`GLMM ${fam} ${id}: aNova no ajusta`); continue; }
    const R = e.indep;
    // Cerca de τ = 0 la verosimilitud es casi plana: el BFGS de R puede detenerse en τ ≈ 0,003–0,007 con una verosimilitud algo MENOR que la de
    // aNova, que llega al límite. Se acepta mientras la verosimilitud de aNova no sea peor y μ y el EE coincidan.
    const cercaDeCero = R.limite || R.tau < 0.02;
    cmp(`GLMM ${fam} · vs R independiente · log-verosimilitud`, id, r.logVerosimilitud, R.loglik, 1e-9, cercaDeCero ? 2e-4 : 2e-5);
    cmp(`GLMM ${fam} · vs R independiente · μ`, id, r.mu, R.mu, 1e-4, 5e-4);
    if (r.limites.includes("tauAlto")) {
      cmp(`GLMM ${fam} · τ en el tope (30) · mismo tope`, id, R.tau > 29 ? 1 : 0, 1, 1e-12, 1e-12);
    } else if (r.limites.includes("tau")) {
      cmp(`GLMM ${fam} · τ en el límite inferior · R también cerca de 0`, id, cercaDeCero ? 1 : 0, 1, 1e-12, 1e-12);
      cmp(`GLMM ${fam} · vs R independiente · EE de μ`, id, r.se, R.se, 5e-3, 2e-4);
    } else if (!cercaDeCero) {
      cmp(`GLMM ${fam} · vs R independiente · τ`, id, r.tau, R.tau, 2e-3, 2e-3);
      cmp(`GLMM ${fam} · vs R independiente · EE de μ`, id, r.se, R.se, 2e-3, 1e-4);
    } else {
      cmp(`GLMM ${fam} · vs R independiente · EE de μ`, id, r.se, R.se, 5e-3, 2e-4);
    }
    if (e.metafor) {
      const M = e.metafor;
      cmp(`GLMM ${fam} · vs metafor::rma.glmm (comprobación) · μ`, id, r.mu, arr(M.mu)[0], 5e-3, 5e-3);
      // metafor a veces no devuelve el error estándar (Hessiano no invertible): entonces no hay con qué comparar
      if (M.se !== null && M.se !== undefined) cmp(`GLMM ${fam} · vs metafor::rma.glmm (comprobación) · EE de μ`, id, r.se, M.se, 0.1, 5e-3);
      cmp(`GLMM ${fam} · vs metafor::rma.glmm (comprobación) · τ²`, id, r.tau2, M.tau2, 0.1, 0.02);
    }
  }
  for (const [fam, bloque] of Object.entries(ref.bruta as J)) {
    const b = bloque as J;
    const filas = b.datos as J[];
    const esProp = fam === "proporcion";
    const nombre = esProp ? "proporciones" : "OR";
    (b.puntos as J[]).forEach((pt, i) => {
      const par = arr(pt.par);
      const tau = Math.exp(par[1]);
      const ll = esProp
        ? logVerosimilitudProporcion(filas.map((d) => ({ x: d.x, n: d.n })), par[0], tau)
        : logVerosimilitudOdds(filas.map((d) => ({ eventosExp: d.a, totalExp: d.n1, eventosCtrl: d.c, totalCtrl: d.n2 })), par[0], tau);
      cmp(`GLMM ${nombre} · log-verosimilitud (GH 21) vs R GH 21 independiente`, `punto ${i + 1}`, ll, pt.gh, 1e-11, 1e-8);
      cmp(`GLMM ${nombre} · log-verosimilitud (GH 21) vs fuerza bruta (integrate)`, `punto ${i + 1}`, ll, pt.bruta, 1e-8, 5e-7);
    });
  }
  return { conjuntos: (ref.conjuntos as J[]).length, metafor: ref.paquetes.metafor as string, BiasedUrn: ref.paquetes.BiasedUrn as string };
}
