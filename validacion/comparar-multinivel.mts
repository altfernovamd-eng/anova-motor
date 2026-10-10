// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Validación diferencial del modelo multinivel con varianza robusta (`multinivel.ts`) contra R.
 * La invoca `comparar.mts` con sus ayudantes. Referencia: `referencias/referencia-multinivel.json`
 * (generada por `generar-referencia-multinivel.R`: 60 conjuntos simulados con metafor::rma.mv y clubSandwich).
 *
 * Qué se contrasta:
 *  1. El ajuste completo (REML + CR2 con Satterthwaite) contra `rma.mv` y `coef_test(vcov = "CR2")`: componentes de varianza,
 *     coeficientes, covarianza del modelo, log-verosimilitud restringida, error estándar robusto, grados de libertad, p e IC.
 *  2. La parte robusta de aNova evaluada con los componentes de varianza de metafor (para aislar el estimador sándwich del optimizador):
 *     debe coincidir con clubSandwich a precisión de máquina.
 */
import { readFileSync } from "node:fs";
import { ajustarMultinivel, inferenciaRobusta, type EfectoDependiente } from "../src/multinivel.ts";
import type { AyudasComparar } from "./comparar-avanzada.mts";

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type J = any;

export function compararMultinivel(h: AyudasComparar) {
  const { cmp, sinPar, arr } = h;
  const ref: J = JSON.parse(readFileSync(new URL("./referencias/referencia-multinivel.json", import.meta.url), "utf8"));

  for (const e of ref.conjuntos as J[]) {
    const mod = e.moderador !== "ninguno";
    const datos: EfectoDependiente[] = (e.datos as J[]).map((d) => ({ estudio: String(d.estudio), y: d.y, v: d.v, x: mod ? [d.x1] : undefined }));
    const opc = { rho: e.rho, nombresModeradores: mod ? ["x1"] : [] };
    const id = `${e.id} (${e.moderador}, ρ=${e.rho}, k=${datos.length})`;
    const s2 = arr(e.sigma2);
    const beta = arr(e.beta), vc = arr(e.vcov), se = arr(e.se_cr2), gl = arr(e.gl), p = arr(e.p), lo = arr(e.ic_inf), hi = arr(e.ic_sup);
    const P = beta.length;

    // 1) Ajuste completo
    const a = ajustarMultinivel(datos, opc);
    if (!a) sinPar.push(`multinivel ${id}: aNova no ajusta`);
    else {
      cmp("Multinivel · REML vs metafor::rma.mv · σ² entre estudios", id, a.sigma2Estudio, s2[0], 1e-4, 2e-6);
      cmp("Multinivel · REML vs metafor::rma.mv · σ² entre efectos", id, a.sigma2Efecto, s2[1], 1e-4, 2e-6);
      // metafor se detiene antes de converger del todo (σ² a ~4e-6 del óptimo): su log-verosimilitud puede quedar ~1e-6 por debajo de la de aNova
      cmp("Multinivel · REML vs metafor::rma.mv · log-verosimilitud", id, a.logVerosimilitudREML, e.ll, 1e-8, 5e-6);
      a.coeficientes.forEach((c, i) => {
        cmp("Multinivel · REML vs metafor::rma.mv · coeficientes", `${id} β${i}`, c.estimacion, beta[i], 1e-5, 1e-6);
        cmp("Multinivel · REML vs metafor::rma.mv · EE del modelo", `${id} β${i}`, c.seModelo, Math.sqrt(vc[i * P + i]), 1e-5, 1e-6);
        cmp("Multinivel · ajuste completo vs clubSandwich · EE robusto", `${id} β${i}`, c.seRobusto, se[i], 1e-4, 1e-6);
        cmp("Multinivel · ajuste completo vs clubSandwich · gl de Satterthwaite", `${id} β${i}`, c.gl, gl[i], 1e-3, 1e-4);
        cmp("Multinivel · ajuste completo vs clubSandwich · valor p", `${id} β${i}`, c.pValor, p[i], 1e-3, 1e-6);
        cmp("Multinivel · ajuste completo vs clubSandwich · IC", `${id} β${i}`, c.icInferior, lo[i], 1e-4, 1e-5);
        cmp("Multinivel · ajuste completo vs clubSandwich · IC", `${id} β${i}`, c.icSuperior, hi[i], 1e-4, 1e-5);
      });
    }

    // 2) Sándwich aislado: con los componentes de varianza de metafor
    const r = inferenciaRobusta(datos, s2[0], s2[1], opc);
    if (!r) sinPar.push(`multinivel (sándwich) ${id}: aNova no calcula`);
    else {
      r.beta.forEach((b, i) => {
        cmp("Multinivel · sándwich CR2 aislado vs clubSandwich · EE robusto", `${id} β${i}`, r.seRobusto[i], se[i], 1e-8, 1e-10);
        cmp("Multinivel · sándwich CR2 aislado vs clubSandwich · gl de Satterthwaite", `${id} β${i}`, r.gl[i], gl[i], 1e-7, 1e-8);
        cmp("Multinivel · sándwich CR2 aislado vs clubSandwich · coeficientes", `${id} β${i}`, b, beta[i], 1e-8, 1e-10);
      });
    }
  }
  return { conjuntos: (ref.conjuntos as J[]).length, metafor: ref.paquetes.metafor as string, clubSandwich: ref.paquetes.clubSandwich as string };
}
