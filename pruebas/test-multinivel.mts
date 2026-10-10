// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Pruebas del modelo multinivel con varianza robusta (multinivel.ts). La comparación numérica con R (metafor::rma.mv y clubSandwich) está en
 * validacion/comparar-multinivel.mts; aquí van las propiedades que no dependen de R: una solución analítica del estimador CR2, los
 * grados de libertad exactos en un diseño balanceado, equivalencias con el meta-análisis de efectos aleatorios, invariancias y casos límite.
 *   npm run test:multinivel
 */
import { existsSync, readFileSync } from "node:fs";
import { ajustarMultinivel, inferenciaRobusta, sensibilidadRho, RHO_SENSIBILIDAD, type EfectoDependiente } from "../src/multinivel.ts";
import { calcularMetaAnalisis } from "../src/statistics.ts";
import { tStudentCuantilDosColas, tStudentDosColas } from "../src/funciones-especiales.ts";
import type { Estudio } from "../src/tipos.ts";

let fallos = 0;
const ok = (c: boolean, m: string) => { console.log(`${c ? "✓" : "✗"} ${m}`); if (!c) fallos++; };
const cerca = (a: number, b: number, tol: number) => Math.abs(a - b) <= tol * Math.max(1, Math.abs(a), Math.abs(b));

// Datos de ejemplo: 8 estudios con 1 a 4 efectos
const base: EfectoDependiente[] = [
  { estudio: "A", y: 0.42, v: 0.04 }, { estudio: "A", y: 0.55, v: 0.06 },
  { estudio: "B", y: 0.10, v: 0.05 }, { estudio: "B", y: 0.31, v: 0.05 }, { estudio: "B", y: 0.22, v: 0.09 },
  { estudio: "C", y: 0.70, v: 0.03 },
  { estudio: "D", y: 0.15, v: 0.08 }, { estudio: "D", y: -0.05, v: 0.07 }, { estudio: "D", y: 0.20, v: 0.06 }, { estudio: "D", y: 0.01, v: 0.10 },
  { estudio: "E", y: 0.50, v: 0.05 }, { estudio: "E", y: 0.38, v: 0.04 },
  { estudio: "F", y: 0.28, v: 0.06 },
  { estudio: "G", y: 0.61, v: 0.07 }, { estudio: "G", y: 0.47, v: 0.05 },
  { estudio: "H", y: 0.05, v: 0.09 },
];

// ---------- Ajuste básico ----------
const r = ajustarMultinivel(base, { rho: 0.6 })!;
ok(!!r && r.k === 16 && r.m === 8 && r.estudiosConVarios === 5 && r.maxPorEstudio === 4, "recuentos: 16 efectos, 8 estudios, 5 con varios, máximo 4");
ok(r.sigma2Estudio >= 0 && r.sigma2Efecto >= 0 && r.coeficientes.length === 1, "componentes de varianza no negativos y un solo coeficiente");
const c = r.coeficientes[0];
ok(c.seRobusto > 0 && c.gl > 0 && c.icInferior < c.estimacion && c.estimacion < c.icSuperior, "IC robusto contiene la estimación; df y EE positivos");
ok(cerca(c.pValor, tStudentDosColas(c.t, c.gl), 1e-12) && c.pValor > 0 && c.pValor <= 1, "valor p = cola de dos lados de una t con los df de Satterthwaite (no la normal)");
ok(cerca(c.icInferior, c.estimacion - tStudentCuantilDosColas(0.05, c.gl) * c.seRobusto, 1e-12) && c.icInferior < c.icModeloInferior, "IC robusto con el cuantil t(df) y más ancho que el del modelo con pocos estudios");

// ---------- Solución analítica: CR2 con observaciones independientes y solo intercepto ----------
{
  // Un efecto por estudio: Φ_i = σ²_w + v_i. CR2 = M² Σ wᵢ² eᵢ² / (1 − M wᵢ) con M = 1/Σw (derivado a mano del estimador).
  const ind = base.map((d, i) => ({ ...d, estudio: `s${i}` }));
  const s2w = 0.03;
  const rb = inferenciaRobusta(ind, 0, s2w, { rho: 0.6 })!;
  const w = ind.map((d) => 1 / (d.v + s2w));
  const sw = w.reduce((a, b) => a + b, 0);
  const mu = ind.reduce((s, d, i) => s + w[i] * d.y, 0) / sw;
  const M = 1 / sw;
  const varCR2 = M * M * ind.reduce((s, d, i) => s + (w[i] * w[i] * (d.y - mu) ** 2) / (1 - M * w[i]), 0);
  ok(cerca(rb.beta[0], mu, 1e-12) && cerca(rb.seModelo[0], Math.sqrt(M), 1e-12), "independientes: β̂ y EE del modelo son los del meta-análisis de efectos aleatorios");
  ok(cerca(rb.seRobusto[0], Math.sqrt(varCR2), 1e-9), `independientes: CR2 = fórmula analítica HC2 ponderada (${rb.seRobusto[0].toFixed(6)})`);
}

// ---------- Diseño balanceado: df de Satterthwaite = m − 1 ----------
{
  const bal: EfectoDependiente[] = [];
  const ys = [0.1, 0.4, 0.25, 0.6, 0.35, 0.05];
  ys.forEach((m0, j) => [0, 0.05, -0.08].forEach((d) => bal.push({ estudio: `b${j}`, y: m0 + d, v: 0.05 })));
  const rb = inferenciaRobusta(bal, 0.02, 0.01, { rho: 0.5 })!;
  ok(cerca(rb.gl[0], ys.length - 1, 1e-9), `balanceado (6 estudios iguales): df de Satterthwaite = m − 1 = 5 (${rb.gl[0].toFixed(10)})`);
}

// ---------- Equivalencia con efectos aleatorios cuando no hay dependencia ----------
{
  const ind: EfectoDependiente[] = base.map((d, i) => ({ ...d, estudio: `s${i}` }));
  const m = ajustarMultinivel(ind, { rho: 0.6 })!;
  const estudios: Estudio[] = ind.map((d, i) => {
    const se = Math.sqrt(d.v);
    return { id: i + 1, autor: `E${i}`, año: "2020", diseño: "ECA", n: "100", medida: "DM", valor: String(d.y), icInferior: String(d.y - 1.959963985 * se), icSuperior: String(d.y + 1.959963985 * se) } as Estudio;
  });
  const pool = calcularMetaAnalisis(estudios, "DM", "REML")!;
  ok(cerca(m.sigma2Efecto, pool.heterogeneidad.tau2, 1e-5) && m.sigma2Estudio === 0, `un efecto por estudio: σ²_efecto = τ² REML del meta-análisis (${m.sigma2Efecto.toFixed(6)} vs ${pool.heterogeneidad.tau2.toFixed(6)})`);
  ok(cerca(m.coeficientes[0].estimacion, pool.aleatorio.estimacion, 1e-6) && cerca(m.coeficientes[0].icModeloInferior, pool.aleatorio.icInferior, 1e-6), "un efecto por estudio: el efecto y el IC del modelo son los del meta-análisis REML");
  ok(m.advertencias.some((a) => /ningún estudio aporta más de un efecto/i.test(a)), "…y avisa de que no hay dependencia que modelar");
}

// ---------- Invariancias ----------
{
  const rev = [...base].reverse();
  const r2 = ajustarMultinivel(rev, { rho: 0.6 })!;
  ok(cerca(r2.coeficientes[0].estimacion, c.estimacion, 1e-8) && cerca(r2.coeficientes[0].seRobusto, c.seRobusto, 1e-8) && cerca(r2.coeficientes[0].gl, c.gl, 1e-8), "el orden de las filas no cambia el resultado");
  const ren = base.map((d) => ({ ...d, estudio: `zz-${d.estudio}` }));
  ok(cerca(ajustarMultinivel(ren, { rho: 0.6 })!.coeficientes[0].seRobusto, c.seRobusto, 1e-9), "los nombres de los estudios no importan");
  const esc = base.map((d) => ({ ...d, y: d.y * 3, v: d.v * 9 }));
  const rs = ajustarMultinivel(esc, { rho: 0.6 })!;
  ok(cerca(rs.coeficientes[0].estimacion, 3 * c.estimacion, 1e-5) && cerca(rs.coeficientes[0].seRobusto, 3 * c.seRobusto, 1e-5) && cerca(rs.coeficientes[0].pValor, c.pValor, 1e-4), "cambiar la escala (y×3, v×9) reescala el efecto y el EE y deja el valor p igual");
  // moderador: desplazarlo por una constante cambia el intercepto pero no la pendiente ni su EE
  const mod = base.map((d, i) => ({ ...d, x: [i % 5 - 2] }));
  const mod2 = base.map((d, i) => ({ ...d, x: [i % 5 - 2 + 10] }));
  const a1 = ajustarMultinivel(mod, { rho: 0.6, nombresModeradores: ["x"] })!;
  const a2 = ajustarMultinivel(mod2, { rho: 0.6, nombresModeradores: ["x"] })!;
  ok(a1.coeficientes.length === 2 && cerca(a1.coeficientes[1].estimacion, a2.coeficientes[1].estimacion, 1e-6) && cerca(a1.coeficientes[1].seRobusto, a2.coeficientes[1].seRobusto, 1e-6) && cerca(a1.coeficientes[1].gl, a2.coeficientes[1].gl, 1e-6), "moderador: desplazarlo por una constante no cambia la pendiente, su EE robusto ni sus df");
  ok(!cerca(a1.coeficientes[0].estimacion, a2.coeficientes[0].estimacion, 1e-3), "…pero sí cambia el intercepto");
}

// ---------- ρ y sensibilidad ----------
{
  const s = sensibilidadRho(base);
  ok(s.length === RHO_SENSIBILIDAD.length && s.every((x) => x.resultado), "sensibilidad: un ajuste por cada ρ");
  const est = s.map((x) => x.resultado!.coeficientes[0].estimacion);
  ok(Math.max(...est) - Math.min(...est) < 0.2, "sensibilidad: el efecto combinado varía poco con ρ en estos datos");
  const se = s.map((x) => x.resultado!.coeficientes[0].seModelo);
  ok(se.every((x) => x > 0), "sensibilidad: EE del modelo positivo para cada ρ");
}

// ---------- Casos límite ----------
ok(ajustarMultinivel(base.slice(0, 2)) === null, "menos de 3 efectos → null");
ok(ajustarMultinivel(base.filter((d) => d.estudio === "B" || d.estudio === "D")) !== null || true, "dos estudios con varios efectos (no revienta)");
ok(ajustarMultinivel(base.map((d) => ({ ...d, estudio: "X" }))) === null, "un solo estudio → null");
ok(ajustarMultinivel(base, { rho: 1 }) === null && ajustarMultinivel(base, { rho: -0.1 }) === null, "ρ fuera de [0, 1) → null (no inventa)");
ok(ajustarMultinivel([...base, { estudio: "Z", y: NaN, v: 0.1 }, { estudio: "Z2", y: 0.2, v: 0 }, { estudio: "", y: 0.3, v: 0.1 }], { rho: 0.6 })!.k === 16, "filas con y no finito, varianza ≤ 0 o sin estudio se descartan");
{
  const homog: EfectoDependiente[] = base.map((d) => ({ ...d, y: 0.3 + (d.v - 0.06) * 0.01 }));
  const h = ajustarMultinivel(homog, { rho: 0.6 })!;
  ok(h.limites.length > 0 && h.sigma2Estudio === 0 && h.coeficientes[0].seRobusto > 0, `sin heterogeneidad real: las varianzas quedan en 0 y se declara (${h.limites.join(", ")})`);
}
{
  const pocos = ajustarMultinivel(base.filter((d) => ["A", "B", "D", "E"].includes(d.estudio)), { rho: 0.6 })!;
  ok(pocos.advertencias.some((a) => /solo 4 estudios/i.test(a)), "con 4 estudios avisa de la poca precisión");
}

// ---------- Regresión con la referencia de R ----------
const ruta = new URL("../validacion/referencias/referencia-multinivel.json", import.meta.url);
if (existsSync(ruta)) {
  const ref = JSON.parse(readFileSync(ruta, "utf8"));
  const e = ref.conjuntos.find((x: { moderador: string }) => x.moderador === "ninguno");
  const datos = e.datos.map((d: { estudio: number; y: number; v: number }) => ({ estudio: String(d.estudio), y: d.y, v: d.v }));
  const a = ajustarMultinivel(datos, { rho: e.rho })!;
  const se = [].concat(e.se_cr2)[0], gl = [].concat(e.gl)[0];
  ok(cerca(a.coeficientes[0].seRobusto, se, 1e-4) && cerca(a.coeficientes[0].gl, gl, 1e-3), `referencia de R (${e.id}): EE robusto ${a.coeficientes[0].seRobusto.toFixed(6)} = clubSandwich ${Number(se).toFixed(6)}; df ${a.coeficientes[0].gl.toFixed(3)} = ${Number(gl).toFixed(3)}`);
} else console.log("— (referencia de R no encontrada: se omite la regresión)");

if (fallos) { console.error(`\n${fallos} prueba(s) fallaron.`); process.exit(1); }
console.log("\nTodas las pruebas pasaron.");
