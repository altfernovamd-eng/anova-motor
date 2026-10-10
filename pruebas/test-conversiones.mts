// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Conversiones de datos del Cochrane Handbook (cap. 6 y 10.6) contra R.
 *
 *   npm run test:conversiones
 *
 * Referencia: validacion/referencias/referencia-conversiones.json (generar-referencia-conversiones.R). Las identidades
 * exactas se prueban contra datos simulados con la verdad conocida (R calcula el IC, la t y el p con t.test y aNova debe
 * recuperar la DE verdadera); la combinación de grupos y la DE del cambio, contra el cálculo directo sobre los datos completos;
 * las aproximaciones de Wan et al. (2014), contra `meta::metamean` como caja negra.
 */
import { readFileSync } from "node:fs";
import {
  combinarGrupos,
  correlacionDesdeSdCambio,
  dmeDesdeOR,
  eeDesdeICAbsoluto,
  eeDesdePAbsoluto,
  eeDiferenciaDesdeIC,
  eeDiferenciaDesdeP,
  eeDiferenciaDesdeT,
  eeLogDesdeICRazon,
  eeLogDesdePRazon,
  icDesdeEEAbsoluto,
  icDesdeEELogRazon,
  mediaSdWan,
  orDesdeDME,
  sdCambioDesdeCorrelacion,
  sdDesdeEE,
  sdDesdeEEDiferencia,
  sdDesdeICMedia,
  sdDesdeIQRNormal,
} from "../src/conversiones-datos.ts";

const ref = JSON.parse(readFileSync(new URL("../validacion/referencias/referencia-conversiones.json", import.meta.url), "utf8"));

let fallos = 0;
let comprobaciones = 0;
const ok = (cond: boolean, msg: string) => {
  comprobaciones++;
  if (!cond) {
    fallos++;
    console.log("  ✗", msg);
  }
};
const rel = (a: number, b: number) => (b === 0 ? Math.abs(a) : Math.abs(a - b) / Math.abs(b));

const peor: Record<string, number> = {};
function mide(clave: string, valor: number, esperado: number, tol: number, contexto: string) {
  const e = rel(valor, esperado);
  peor[clave] = Math.max(peor[clave] ?? 0, e);
  ok(e <= tol, `${clave}: error relativo ${e.toExponential(2)} (aNova ${valor}, R ${esperado}) en ${contexto}`);
}

console.log("— Contra datos simulados con la verdad conocida —");
for (const c of ref.icMedia) {
  const ctx = `n=${c.n}, nivel=${c.nivel}`;
  const r = sdDesdeICMedia(c.inf, c.sup, c.n, { nivel: c.nivel, distribucion: "t" })!;
  mide("DE desde el IC de una media (t)", r.sd, c.sd, 1e-9, ctx);
  mide("EE desde el IC de una media (t)", r.se, c.ee, 1e-9, ctx);
  const rz = sdDesdeICMedia(c.infZ, c.supZ, c.n, { nivel: c.nivel, distribucion: "normal" })!;
  mide("DE desde el IC de una media (normal)", rz.sd, c.sd, 1e-9, ctx);
  mide("DE desde el EE", sdDesdeEE(c.ee, c.n)!, c.sd, 1e-12, ctx);
}
for (const d of ref.dif) {
  const ctx = `n=${d.n1}/${d.n2}, nivel=${d.nivel}`;
  const ee1 = eeDiferenciaDesdeIC(d.inf, d.sup, d.n1, d.n2, { nivel: d.nivel })!;
  mide("EE de la diferencia desde su IC", ee1, d.ee, 1e-9, ctx);
  mide("EE de la diferencia desde la t", eeDiferenciaDesdeT(d.md, d.t)!, d.ee, 1e-12, ctx);
  if (d.p > 1e-12) mide("EE de la diferencia desde el p", eeDiferenciaDesdeP(d.md, d.p, d.n1, d.n2)!, d.ee, 1e-6, ctx);
  mide("DE intragrupo desde el EE de la diferencia", sdDesdeEEDiferencia(d.ee, d.n1, d.n2)!, d.sp, 1e-9, ctx);
  mide("DE intragrupo desde el IC de la diferencia", sdDesdeEEDiferencia(ee1, d.n1, d.n2)!, d.sp, 1e-9, ctx);
}
for (const c of ref.comb) {
  const r = combinarGrupos({ n: c.n1, media: c.m1, sd: c.sd1 }, { n: c.n2, media: c.m2, sd: c.sd2 })!;
  const ctx = `n=${c.n1}+${c.n2}`;
  ok(r.n === c.n, `combinar grupos: N ${r.n} ≠ ${c.n} en ${ctx}`);
  mide("combinar grupos: media", r.media, c.m, 1e-12, ctx);
  mide("combinar grupos: DE", r.sd, c.sd, 1e-10, ctx);
}
for (const c of ref.cambio) {
  const ctx = `sdBase=${c.sdBase.toFixed(2)}`;
  mide("DE del cambio desde la correlación", sdCambioDesdeCorrelacion(c.sdBase, c.sdFinal, c.cor)!, c.sdCambio, 1e-9, ctx);
  mide("correlación desde las tres DE", correlacionDesdeSdCambio(c.sdBase, c.sdFinal, c.sdCambio)!, c.cor, 1e-9, ctx);
}

console.log("— Wan et al. (2014) contra meta::metamean (caja negra) —");
const tolWan = (n: number) => (n >= 400 ? 1e-9 : n > 50 ? 3e-3 : null); // n ≤ 50: meta usa las tablas del artículo (ver nota)
for (const w of ref.wan) {
  const ctx = `n=${w.n}`;
  const c1 = mediaSdWan({ n: w.n, mediana: w.m, min: w.a, max: w.b })!;
  const c3 = mediaSdWan({ n: w.n, mediana: w.m, q1: w.q1, q3: w.q3 })!;
  const c2 = mediaSdWan({ n: w.n, mediana: w.m, q1: w.q1, q3: w.q3, min: w.a, max: w.b })!;
  ok(c1.escenario === "min-mediana-max" && c3.escenario === "cuartiles" && c2.escenario === "cinco-numeros", `Wan: escenarios correctos en ${ctx}`);
  // La media no depende de tablas: coincide con meta para cualquier n.
  mide("Wan · media (mín, mediana, máx)", c1.media, w.media1, 1e-12, ctx);
  mide("Wan · media (cuartiles)", c3.media, w.media2, 1e-12, ctx);
  mide("Wan · media (cinco números)", c2.media, w.media3, 1e-12, ctx);
  const tol = tolWan(w.n);
  if (tol !== null) {
    mide("Wan · DE (mín, mediana, máx), n > 50", c1.sd, w.sd1, tol, ctx);
    mide("Wan · DE (cuartiles), n > 50", c3.sd, w.sd2, tol, ctx);
    mide("Wan · DE (cinco números), n > 50", c2.sd, w.sd3, tol, ctx);
  }
  // Con n pequeño la aproximación de Blom y la tabla del artículo difieren, pero ambas son aproximaciones del mismo valor.
  if (w.n <= 50 && w.n >= 13) {
    ok(rel(c1.sd, w.sd1) < 0.02 && rel(c3.sd, w.sd2) < 0.02, `Wan: con n=${w.n} la diferencia con meta es menor del 2 % (las tablas del artículo frente a Blom)`);
  }
}

console.log("— Fórmulas simples contra R —");
for (const o of ref.orDme) {
  const r = dmeDesdeOR(o.or, o.ee);
  if (o.or === 1) {
    ok(!!r && Math.abs(r.dme) < 1e-15, "OR = 1 → DME = 0");
    continue;
  }
  mide("OR → DME", r!.dme, o.dme, 1e-12, `OR=${o.or}`);
  mide("OR → EE de la DME", r!.ee, o.eeDme, 1e-12, `OR=${o.or}`);
  const v = orDesdeDME(o.dme, o.eeDme)!;
  mide("DME → OR (ida y vuelta)", v.or, o.or, 1e-12, `OR=${o.or}`);
  mide("DME → EE del log OR (ida y vuelta)", v.eeLogOR, o.ee, 1e-12, `OR=${o.or}`);
}
for (const e of ref.eeIc) {
  mide("EE desde el IC de una medida absoluta", eeDesdeICAbsoluto(e.inf, e.sup, e.nivel)!, e.eeAbs, 1e-10, `nivel=${e.nivel}`);
  mide("EE del log desde el IC de una razón", eeLogDesdeICRazon(e.inf, e.sup, e.nivel)!, e.eeLog, 1e-10, `nivel=${e.nivel}`);
}
for (const e of ref.eeP) {
  mide("EE desde el p de una medida absoluta", eeDesdePAbsoluto(e.est, e.p)!, e.eeAbs, 1e-9, `p=${e.p}`);
  mide("EE del log desde el p de una razón", eeLogDesdePRazon(e.razon, e.p)!, e.eeLog, 1e-9, `p=${e.p}`);
}

console.log("— IC desde la estimación y el P —");
{
  // Ida y vuelta: del IC al EE y de ahí otra vez al IC; y de la estimación y el P al IC que deja ese mismo P en el límite.
  const abs = icDesdeEEAbsoluto(0.9, eeDesdeICAbsoluto(0.2, 1.6, 0.95)!, 0.95)!;
  ok(Math.abs(abs.inf - 0.2) < 1e-12 && Math.abs(abs.sup - 1.6) < 1e-12, "medida absoluta: IC → EE → IC reproduce el IC");
  const raz = icDesdeEELogRazon(Math.sqrt(0.62 * 1.91), eeLogDesdeICRazon(0.62, 1.91)!)!;
  ok(Math.abs(raz.inf - 0.62) < 1e-12 && Math.abs(raz.sup - 1.91) < 1e-12, "razón: IC → EE del log → IC reproduce el IC");
  const ee = eeDesdePAbsoluto(0.9, 0.05)!; // con p = 0,05 el límite inferior del IC 95 % cae en 0
  ok(Math.abs(icDesdeEEAbsoluto(0.9, ee)!.inf) < 1e-9, "con P = 0,05 el IC 95 % de una medida absoluta toca el 0");
  const eeR = eeLogDesdePRazon(1.8, 0.05)!;
  ok(Math.abs(icDesdeEELogRazon(1.8, eeR)!.inf - 1) < 1e-9, "con P = 0,05 el IC 95 % de una razón toca el 1");
  ok(icDesdeEEAbsoluto(1, 0) === null && icDesdeEELogRazon(-1, 0.3) === null, "EE ≤ 0 o razón ≤ 0 → null");
}

console.log("— Casos límite y números a mano —");
ok(Math.abs(sdDesdeIQRNormal(10, 24.5)! - 10.740740740740742) < 1e-12, "IQR / 1,35");
ok(sdDesdeIQRNormal(5, 5) === null && sdDesdeIQRNormal(7, 3) === null, "IQR vacío o invertido → null");
ok(Math.abs(sdDesdeEE(2, 25)! - 10) < 1e-15, "SD = 2 · √25 = 10");
ok(sdDesdeEE(-1, 10) === null && sdDesdeEE(1, 0) === null && sdDesdeEE(Number.NaN, 10) === null, "EE negativo, N = 0 o NaN → null");
ok(sdDesdeICMedia(5, 3, 20) === null && sdDesdeICMedia(3, 5, 1) === null && sdDesdeICMedia(3, 5, 20, { nivel: 1.2 }) === null, "IC invertido, N = 1 o nivel imposible → null");
ok(eeDiferenciaDesdeP(1, 0, 5, 5) === null && eeDiferenciaDesdeP(1, 1, 5, 5) === null && eeDiferenciaDesdeP(0, 0.04, 5, 5) === null, "p fuera de (0,1) o diferencia 0 → null (no se puede recuperar el EE)");
ok(eeDiferenciaDesdeT(2, 0) === null, "t = 0 → null");
ok(combinarGrupos({ n: 10, media: 5, sd: 1 }, { n: 0, media: 5, sd: 1 }) === null, "combinar con un grupo vacío → null");
{
  const igual = combinarGrupos({ n: 10, media: 4, sd: 2 }, { n: 10, media: 4, sd: 2 })!;
  ok(igual.n === 20 && igual.media === 4 && Math.abs(igual.sd - Math.sqrt((9 * 4 + 9 * 4) / 19)) < 1e-12, "dos grupos iguales: la media no cambia y la DE es la agrupada");
  const lejanos = combinarGrupos({ n: 10, media: 0, sd: 1 }, { n: 10, media: 10, sd: 1 })!;
  ok(lejanos.sd > 5, "dos grupos con medias muy distintas: la DE combinada incluye la diferencia entre medias");
}
ok(sdCambioDesdeCorrelacion(10, 12, 1.5) === null && sdCambioDesdeCorrelacion(10, 12, -1.1) === null, "correlación fuera de [−1, 1] → null");
ok(Math.abs(sdCambioDesdeCorrelacion(10, 10, 1)!) < 1e-12, "correlación 1 con DE iguales → DE del cambio 0");
ok(Math.abs(sdCambioDesdeCorrelacion(6, 8, 0)! - 10) < 1e-12, "correlación 0 → √(6² + 8²) = 10");
ok(correlacionDesdeSdCambio(5, 5, 50) === null, "tres DE incoherentes (correlación fuera de rango) → null");
ok(mediaSdWan({ n: 4, mediana: 5, q1: 3, q3: 8 }) === null, "n < 5 → null");
ok(mediaSdWan({ n: 30, mediana: 5, q1: 8, q3: 3 }) === null, "cuartiles invertidos → null");
ok(mediaSdWan({ n: 30, mediana: 5, q1: 3 }) === null, "un solo cuartil → null");
ok(mediaSdWan({ n: 30, mediana: 5 }) === null, "solo la mediana → null (no se inventa media ni DE)");
ok(mediaSdWan({ n: 30, mediana: 5, min: 1, q1: 3, q3: 8 }) === null, "combinación incompleta (mín y cuartiles sin máx) → null");
ok(mediaSdWan({ n: 10, mediana: 5, q1: 3, q3: 8 })!.advertencias.length > 0, "n < 15 → advierte");
ok(mediaSdWan({ n: 40, mediana: 5, min: 1, max: 11 })!.advertencias.some((a) => /rango/i.test(a)), "solo el rango → advierte que es inestable");
ok(dmeDesdeOR(-1, 0.3) === null && dmeDesdeOR(2, 0) === null && orDesdeDME(0.3, -1) === null, "OR ≤ 0 o EE ≤ 0 → null");

console.log("\nPeor error relativo por conversión:");
for (const [k, v] of Object.entries(peor).sort()) console.log(`  ${k.padEnd(52)} ${v.toExponential(2)}`);
console.log(`\n${comprobaciones - fallos}/${comprobaciones} comprobaciones correctas`);
if (fallos > 0) {
  console.log(`✗ ${fallos} fallos`);
  process.exit(1);
}
console.log("✓ Conversiones de datos de acuerdo con R");
