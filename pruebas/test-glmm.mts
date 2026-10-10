// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Pruebas de los modelos mixtos logísticos (glmm.ts), del modelo multinivel integrado en el Paso 7 y de su llegada al manuscrito. La
 * comparación numérica con R está en validacion/comparar-glmm.mts y comparar-multinivel.mts; aquí van las propiedades que no
 * dependen de R: soluciones calculadas a mano, simetrías, casos límite, honestidad de lo que se declara y la integración con los estudios.
 *   npm run test:glmm
 */
import { existsSync, readFileSync } from "node:fs";
import {
  ajustarGLMMOdds,
  ajustarGLMMProporcion,
  logVerosimilitudOdds,
  logVerosimilitudProporcion,
  type DatosTabla,
} from "../src/glmm.ts";
import { metaAnalisisMH } from "../src/statistics.ts";
import { normalDosColas } from "../src/funciones-especiales.ts";
import type { Estudio } from "../src/tipos.ts";

let fallos = 0;
const ok = (c: boolean, m: string) => { console.log(`${c ? "✓" : "✗"} ${m}`); if (!c) fallos++; };
const cerca = (a: number, b: number, tol: number) => Math.abs(a - b) <= tol * Math.max(1, Math.abs(a), Math.abs(b));
const expit = (x: number) => 1 / (1 + Math.exp(-x));

const lnC = (n: number, k: number) => { let s = 0; for (let i = 1; i <= k; i++) s += Math.log((n - k + i) / i); return s; };
const T = (a: number, n1: number, c: number, n2: number): DatosTabla => ({ eventosExp: a, totalExp: n1, eventosCtrl: c, totalCtrl: n2 });

// ---------- Verosimilitud: soluciones a mano ----------
{
  const x = [3, 9, 0, 14, 5], n = [20, 31, 12, 40, 18];
  const datos = x.map((v, i) => ({ x: v, n: n[i] }));
  const mu = -0.8;
  let directo = 0;
  for (let i = 0; i < x.length; i++) directo += lnC(n[i], x[i]) + x[i] * Math.log(expit(mu)) + (n[i] - x[i]) * Math.log(1 - expit(mu));
  ok(cerca(logVerosimilitudProporcion(datos, mu, 1e-6), directo, 1e-9), `τ→0: la verosimilitud marginal es la binomial simple (${directo.toFixed(6)})`);
  ok(cerca(logVerosimilitudProporcion(datos, mu, 0.9, 21), logVerosimilitudProporcion(datos, mu, 0.9, 41), 1e-9), "21 y 41 nodos dan la misma verosimilitud (cuadratura convergida)");
  let lanza = false; try { logVerosimilitudProporcion(datos, mu, 0.9, 101); } catch { lanza = true; }
  ok(lanza, "pedir más de 64 nodos falla con un error claro en vez de dar un número corrupto");
  // integrando estrecho (muchos pacientes, τ grande): solo converge si la cuadratura se centra bien en el máximo (cuadratura adaptativa)
  const duro = [T(100, 2000, 150, 2000), T(40, 900, 70, 950)];
  ok(cerca(logVerosimilitudOdds(duro, 0.3, 2, 21), logVerosimilitudOdds(duro, 0.3, 2, 61), 1e-10), "integrando estrecho (n = 2000, τ = 2): 21 y 61 nodos coinciden gracias al centrado adaptativo");
  const duroP = [{ x: 400, n: 2000 }, { x: 350, n: 1800 }];
  ok(cerca(logVerosimilitudProporcion(duroP, -1.3, 2, 21), logVerosimilitudProporcion(duroP, -1.3, 2, 61), 1e-10), "…y en proporciones");
  // un estudio, hipergeométrica: P(a=2 | n1=5, n2=5, s=3) con log-OR θ
  const th = 0.7;
  const u = [0, 1, 2, 3];
  const z = u.reduce((s, v) => s + Math.exp(lnC(5, v) + lnC(5, 3 - v) + th * v), 0);
  const p = Math.exp(lnC(5, 2) + lnC(5, 1) + th * 2) / z;
  ok(cerca(logVerosimilitudOdds([T(2, 5, 1, 5)], th, 1e-6), Math.log(p), 1e-9), `hipergeométrica no central calculada a mano: P = ${p.toFixed(6)}`);
  // con θ = 0 es la hipergeométrica central
  const pc = (lnC(5, 2) + lnC(5, 1)) - lnC(10, 3);
  ok(cerca(logVerosimilitudOdds([T(2, 5, 1, 5)], 0, 1e-6), pc, 1e-9), "con OR = 1 es la hipergeométrica central C(5,2)C(5,1)/C(10,3)");
}

// ---------- Ajuste: proporciones ----------
const xs = [11, 31, 41, 72, 15, 95, 68, 26, 11, 15], ns = [68, 87, 123, 183, 56, 182, 190, 139, 133, 31];
const prop = ajustarGLMMProporcion(xs.map((x, i) => ({ x, n: ns[i] })))!;
ok(!!prop && prop.modelo === "binomial-normal" && prop.k === 10 && prop.kInformativos === 10, "proporciones: k = 10 estudios");
ok(prop.tau2 > 0 && prop.limites.length === 0 && prop.se > 0 && prop.pValor === null, "τ² > 0, error estándar positivo y sin valor p (no hay nulo natural)");
ok(cerca(prop.estimacion, expit(prop.mu), 1e-12) && cerca(prop.icInferior, expit(prop.muIc[0]), 1e-12) && cerca((prop.muIc[0] + prop.muIc[1]) / 2, prop.mu, 1e-12), "la proporción es expit(μ) y el IC es simétrico en logit");
{
  // es un máximo: ninguna perturbación mejora la verosimilitud
  const dat = xs.map((x, i) => ({ x, n: ns[i] }));
  const base = logVerosimilitudProporcion(dat, prop.mu, prop.tau);
  const v = [logVerosimilitudProporcion(dat, prop.mu + 0.01, prop.tau), logVerosimilitudProporcion(dat, prop.mu - 0.01, prop.tau), logVerosimilitudProporcion(dat, prop.mu, prop.tau * 1.02), logVerosimilitudProporcion(dat, prop.mu, prop.tau / 1.02)];
  ok(v.every((l) => l <= base + 1e-8) && cerca(base, prop.logVerosimilitud, 1e-9), "es un máximo local de la verosimilitud y la informada es la recalculada");
  const rev = ajustarGLMMProporcion(dat.slice().reverse())!;
  ok(cerca(rev.mu, prop.mu, 1e-5) && cerca(rev.tau, prop.tau, 1e-3), "el orden de los estudios no cambia el ajuste");
}
{
  // sin heterogeneidad: τ en el límite y μ = logit(proporción global)
  const n = [40, 60, 80, 100], x = [10, 15, 20, 25];
  const h = ajustarGLMMProporcion(x.map((v, i) => ({ x: v, n: n[i] })))!;
  ok(h.limites.includes("tau") && cerca(h.mu, Math.log(0.25 / 0.75), 1e-5), `proporciones idénticas: τ en el límite y μ = logit(0,25) = ${h.mu.toFixed(6)}`);
  ok(h.advertencias.some((a) => /límite inferior/i.test(a)) && Number.isFinite(h.se) && h.se > 0, "…y se declara, con error estándar finito");
}
{
  const conCero = ajustarGLMMProporcion([{ x: 0, n: 30 }, { x: 4, n: 50 }, { x: 9, n: 60 }, { x: 2, n: 25 }])!;
  ok(conCero.kInformativos === 4 && Number.isFinite(conCero.mu) && conCero.estimacion > 0 && conCero.estimacion < 0.2, "una proporción con 0 eventos entra sin corrección y el modelo ajusta");
  const malos = ajustarGLMMProporcion([{ x: 1.5, n: 30 }, { x: 5, n: 4 }, { x: 4, n: 50 }, { x: 9, n: 60 }, { x: 2, n: 25 }])!;
  ok(malos.excluidos.length === 2 && malos.k === 3 && malos.excluidos.some((e) => /no enteros/.test(e.motivo)) && malos.excluidos.some((e) => /imposibles/.test(e.motivo)), "recuentos no enteros o imposibles se excluyen con su motivo");
}
ok(ajustarGLMMProporcion([{ x: 1, n: 10 }]) === null, "un solo estudio → null");

// ---------- Ajuste: odds ratio ----------
const tablas = [T(5, 120, 12, 118), T(3, 80, 9, 82), T(0, 50, 4, 52), T(9, 300, 20, 310), T(2, 60, 2, 61), T(7, 150, 15, 149), T(1, 40, 6, 43), T(12, 400, 21, 395)];
const orr = ajustarGLMMOdds(tablas)!;
ok(!!orr && orr.modelo === "hipergeometrica-normal" && orr.k === 8 && orr.kInformativos === 8 && orr.sinInformacion.length === 0, "OR: 8 estudios, todos con información (el de 0 eventos entra sin corrección)");
ok(orr.estimacion < 1 && orr.pValor !== null && cerca(orr.pValor, normalDosColas(orr.mu / orr.se), 1e-12) && cerca(orr.estimacion, Math.exp(orr.mu), 1e-12), "OR < 1 en estos datos; p de Wald consistente con μ/EE");
{
  const inv = ajustarGLMMOdds(tablas.map((t) => T(t.eventosCtrl, t.totalCtrl, t.eventosExp, t.totalExp)))!;
  ok(cerca(inv.mu, -orr.mu, 1e-4) && cerca(inv.tau, orr.tau, 1e-3) && cerca(inv.logVerosimilitud, orr.logVerosimilitud, 1e-8) && cerca(inv.se, orr.se, 1e-3), "intercambiar experimental y control invierte el signo de μ y deja τ, EE y verosimilitud iguales");
  const rev = ajustarGLMMOdds(tablas.slice().reverse())!;
  ok(cerca(rev.mu, orr.mu, 1e-5), "el orden de los estudios no cambia el OR");
  const conDobleCero = ajustarGLMMOdds([...tablas, T(0, 30, 0, 32), T(30, 30, 32, 32)])!;
  ok(conDobleCero.sinInformacion.length === 2 && conDobleCero.k === 10 && conDobleCero.kInformativos === 8 && cerca(conDobleCero.mu, orr.mu, 1e-9) && cerca(conDobleCero.logVerosimilitud, orr.logVerosimilitud, 1e-9), "sin eventos en ningún brazo (o con evento en todos): se declaran y no cambian el ajuste");
  const frac = ajustarGLMMOdds([...tablas, T(3, 40, 5 / 2, 51 / 2)])!;
  ok(frac.excluidos.length === 1 && /no enteros/.test(frac.excluidos[0].motivo), "grupo control dividido (recuentos no enteros) → excluido y declarado, no redondeado");
  ok(ajustarGLMMOdds([T(0, 30, 0, 32), T(0, 20, 0, 25), T(3, 40, 5, 41)]) === null, "menos de 2 estudios con información → null");
}
{
  // muestras grandes y poco heterogéneo: el OR del modelo mixto ≈ Mantel-Haenszel
  const grandes = [T(60, 500, 80, 500), T(45, 400, 62, 410), T(110, 900, 140, 880), T(30, 300, 41, 310), T(75, 600, 99, 590), T(52, 450, 70, 445)];
  const g = ajustarGLMMOdds(grandes)!;
  const est: Estudio[] = grandes.map((t, i) => ({ id: i + 1, autor: `G${i}`, año: "2020", diseño: "ECA", n: String(t.totalExp + t.totalCtrl), medida: "OR", valor: "1", icInferior: "0.5", icSuperior: "2", crudos: { tipo: "2x2", eventosTrat: String(t.eventosExp), totalTrat: String(t.totalExp), eventosCtrl: String(t.eventosCtrl), totalCtrl: String(t.totalCtrl) } }) as Estudio);
  const mh = metaAnalisisMH(est, "OR")!;
  ok(Math.abs(Math.log(g.estimacion) - Math.log(mh.estimacion)) < 0.03, `muestras grandes: OR del modelo mixto ${g.estimacion.toFixed(3)} ≈ Mantel-Haenszel ${mh.estimacion.toFixed(3)}`);
}

// ---------- Valores de la implementación independiente de R (copiados de referencia-glmm.json, para no depender del archivo) ----------
{
  // El EE sale de la inversa del Hessiano COMPLETO (μ, ln τ), no condicionado a τ: con τ grande la diferencia llega al 15 %
  const p30 = [[21, 23], [18, 20], [22, 24], [29, 29], [27, 30], [9, 9], [14, 20], [8, 8], [27, 27], [36, 40], [12, 12], [19, 20]].map(([x, n]) => ({ x, n }));
  const r30 = ajustarGLMMProporcion(p30)!;
  ok(cerca(r30.logVerosimilitud, -19.7065439608742, 2e-5) && cerca(r30.mu, 2.75430057465626, 5e-4) && cerca(r30.tau, 0.767960293314666, 2e-3), `R independiente (prop30): μ ${r30.mu.toFixed(5)}, τ ${r30.tau.toFixed(4)}, logLik ${r30.logVerosimilitud.toFixed(5)}`);
  ok(cerca(r30.se, 0.405387787283757, 2e-3), `R independiente (prop30): EE de μ ${r30.se.toFixed(5)} = 0.40539 (el condicionado a τ daría ≈ 15 % menos)`);
  const o35 = ajustarGLMMOdds([[1, 267, 5, 217], [1, 253, 6, 387], [1, 256, 10, 388], [10, 373, 28, 464], [3, 221, 7, 302], [4, 382, 3, 235], [1, 162, 15, 194], [0, 64, 3, 76], [0, 119, 10, 171], [2, 153, 2, 118], [1, 158, 13, 200], [4, 264, 8, 252]].map(([a, n1, c, n2]) => T(a, n1, c, n2)))!;
  ok(cerca(o35.logVerosimilitud, -19.1998458598464, 2e-5) && cerca(o35.mu, -1.2965399391896, 5e-4) && cerca(o35.se, 0.274661592748168, 2e-3), `R independiente (or35, con ceros): ln OR ${o35.mu.toFixed(5)}, EE ${o35.se.toFixed(5)}`);
}

// ---------- Regresión con la referencia de R ----------
const ruta = new URL("../validacion/referencias/referencia-glmm.json", import.meta.url);
if (existsSync(ruta) && JSON.parse(readFileSync(ruta, "utf8")).conjuntos.length > 0) {
  const ref = JSON.parse(readFileSync(ruta, "utf8"));
  const p = ref.conjuntos.find((c: { id: string }) => c.id.startsWith("prop"));
  const o = ref.conjuntos.find((c: { id: string; k_informativos?: number }) => c.id.startsWith("or") && (c.k_informativos ?? 0) >= 5);
  const rp = ajustarGLMMProporcion(p.datos.map((d: { x: number; n: number }) => ({ x: d.x, n: d.n })))!;
  ok(cerca(rp.logVerosimilitud, p.indep.loglik, 2e-5) && cerca(rp.mu, p.indep.mu, 5e-4), `referencia de R (${p.id}): μ ${rp.mu.toFixed(5)} = ${p.indep.mu.toFixed(5)}`);
  // el error estándar de la inversa del Hessiano completo (no el condicionado a τ): en conjuntos con τ grande difieren
  for (const c of ref.conjuntos as { id: string; datos: { x: number; n: number }[]; indep: { tau: number; se: number; limite: boolean } }[]) {
    if (!c.id.startsWith("prop") || c.indep.limite || c.indep.tau < 1) continue;
    const r1 = ajustarGLMMProporcion(c.datos.map((d) => ({ x: d.x, n: d.n })))!;
    ok(cerca(r1.se, c.indep.se, 3e-3), `referencia de R (${c.id}, τ = ${c.indep.tau.toFixed(2)}): EE de μ ${r1.se.toFixed(5)} = ${c.indep.se.toFixed(5)}`);
    break;
  }
  const ro = ajustarGLMMOdds(o.datos.map((d: { a: number; n1: number; c: number; n2: number }) => T(d.a, d.n1, d.c, d.n2)))!;
  ok(cerca(ro.logVerosimilitud, o.indep.loglik, 2e-5) && cerca(ro.mu, o.indep.mu, 5e-4), `referencia de R (${o.id}): ln OR ${ro.mu.toFixed(5)} = ${o.indep.mu.toFixed(5)}`);
} else console.log("— (referencia de R no encontrada: se omite la regresión)");

if (fallos) { console.error(`\n${fallos} prueba(s) fallaron.`); process.exit(1); }
console.log("\nTodas las pruebas pasaron.");
