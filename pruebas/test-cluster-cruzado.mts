// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
// Comprobación de las fórmulas de ensayos clúster-aleatorizados y cruzados contra el Cochrane Handbook v6, capítulo 23.
// No hay un paquete de R que implemente este capítulo (son fórmulas de ajuste, no un método con software de referencia),
// así que la referencia son: (1) las fórmulas del capítulo, leídas literalmente de sus imágenes; (2) los ejemplos
// resueltos del propio capítulo (23.1.4.1, 23.1.5 y 23.2.7.4); y (3) las equivalencias que el capítulo afirma.
// Ejecutar: npm run test:cluster-cruzado   ·   Detalle y fuente: VALIDACION_ESTADISTICA.md §3.13
import {
  calcularMetaAnalisis,
  efectoDeDisenoCluster,
  efectoDesde2x2,
  efectoDesdeCruzado,
  efectoDesdeMedias,
} from "../src/statistics.ts";
import type { Estudio } from "../src/tipos.ts";

let fallos = 0;
const ok = (c: boolean, m: string) => { console.log(`${c ? "✓" : "✗"} ${m}`); if (!c) fallos++; };
const cerca = (a: number, b: number, tol: number) => Math.abs(a - b) <= tol;

// ------------------------------------------------------------------------------------------------
// 23.1 Clúster — DEFF = 1 + (M − 1) × ICC, M = tamaño medio del clúster (23.1.4)
// ------------------------------------------------------------------------------------------------
// Ejemplo 23.1.4.1: 10 aulas con 295 niños (tratamiento) y 11 con 330 (control); ICC = 0,02.
// M = (295 + 330) / (10 + 11) = 29,76 (el Handbook lo redondea a 29,8 y obtiene 1,576).
const deff = efectoDeDisenoCluster(295 + 330, 10 + 11, 0.02)!;
ok(cerca(deff, 1 + (625 / 21 - 1) * 0.02, 1e-12), "DEFF = 1 + (M − 1)·ICC con M = participantes totales / clústeres totales");
ok(cerca(deff, 1.576, 1e-3), `DEFF del ejemplo 23.1.4.1 = ${deff.toFixed(4)} (Handbook: 1,576, con M redondeado a 29,8)`);
ok(cerca(1 + (29.8 - 1) * 0.02, 1.576, 1e-12), "El 1,576 del Handbook sale de redondear M a 29,8 (comprobación de la aritmética del libro)");

// Validaciones: nunca se ajusta con datos inválidos
ok(efectoDeDisenoCluster(100, 10, -0.01) === null, "ICC negativo → null");
ok(efectoDeDisenoCluster(100, 10, 1.2) === null, "ICC mayor que 1 → null (es una correlación; antes se aceptaba)");
ok(efectoDeDisenoCluster(100, 10, 1) !== null && efectoDeDisenoCluster(100, 10, 0) === 1, "ICC = 1 válido; ICC = 0 → DEFF = 1 (sin ajuste)");
ok(efectoDeDisenoCluster(5, 10, 0.05) === null, "Menos de un participante por clúster (M < 1) → null");
ok(efectoDeDisenoCluster(100, 0, 0.05) === null && efectoDeDisenoCluster(NaN, 10, 0.05) === null, "Sin nº de clústeres o con N no numérico → null");

// 23.1.4 (tamaño efectivo) y 23.1.5 (inflar el EE) son «equivalentes»: se comprueba para cada medida con las celdas SIN redondear.
const dEj = 1.576; // el del Handbook, para reproducir su ejemplo
const t = efectoDesde2x2(63, 295, 84, 330, "OR")!;
ok(cerca(t.valor, 0.795, 5e-4) && cerca(t.icInferior, 0.548, 5e-4) && cerca(t.icSuperior, 1.154, 5e-4), `OR del ejemplo 23.1.5 = ${t.valor.toFixed(3)} (${t.icInferior.toFixed(3)} a ${t.icSuperior.toFixed(3)}); Handbook: 0,795 (0,548 a 1,154)`);
ok(cerca(Math.log(t.valor), -0.23, 5e-3) && cerca(t.se, 0.19, 5e-3), `ln(OR) = ${Math.log(t.valor).toFixed(3)} y EE = ${t.se.toFixed(3)}; Handbook: −0,23 y 0,19`);
ok(cerca(t.se * Math.sqrt(dEj), 0.24, 5e-3), `EE inflado = ${(t.se * Math.sqrt(dEj)).toFixed(3)}; Handbook: 0,24`);
// Datos binarios: eventos y participantes divididos por el mismo DEFF (los eventos del ejemplo: 40,0/187,2 y 53,3/209,4)
ok(cerca(63 / dEj, 40.0, 0.05) && cerca(295 / dEj, 187.2, 0.05) && cerca(84 / dEj, 53.3, 0.05) && cerca(330 / dEj, 209.4, 0.05), "Tamaños efectivos del ejemplo: 40,0/187,2 y 53,3/209,4");
for (const medida of ["OR", "RR", "RAR"] as const) {
  const original = efectoDesde2x2(63, 295, 84, 330, medida)!;
  const efectivo = efectoDesde2x2(63 / dEj, 295 / dEj, 84 / dEj, 330 / dEj, medida)!;
  ok(cerca(efectivo.se, original.se * Math.sqrt(dEj), 1e-12) && cerca(efectivo.valor, original.valor, 1e-12),
    `Equivalencia 23.1.4 ≡ 23.1.5 en ${medida}: EE con celdas/DEFF = EE × √DEFF (mismo efecto)`);
}
// Datos continuos: solo se reduce N; medias y DE no cambian
{
  const original = efectoDesdeMedias(10, 4, 120, 8, 4.5, 130)!;
  const efectivo = efectoDesdeMedias(10, 4, 120 / dEj, 8, 4.5, 130 / dEj)!;
  ok(cerca(efectivo.se, original.se * Math.sqrt(dEj), 1e-12), "Equivalencia 23.1.4 ≡ 23.1.5 en diferencia de medias: EE con N/DEFF = EE × √DEFF");
}

// Lo que el motor hace de verdad con un estudio clúster: infla el EE por √DEFF y deja intacto el efecto
{
  const est = (id: number, v: number, lo: number, hi: number, extra: Partial<Estudio> = {}): Estudio =>
    ({ id, autor: `E${id}`, año: "2020", diseño: "ECA", n: "", medida: "OR", valor: String(v), icInferior: String(lo), icSuperior: String(hi), ...extra }) as Estudio;
  const cluster = est(1, t.valor, t.icInferior, t.icSuperior, { n: "625", varianteEnsayo: { tipo: "cluster", numClusters: "21", icc: "0.02" } });
  const normal = est(2, 0.9, 0.6, 1.35);
  const conAjuste = calcularMetaAnalisis([cluster, normal], "OR", "DL")!;
  const sinAjuste = calcularMetaAnalisis([{ ...cluster, varianteEnsayo: undefined }, normal], "OR", "DL")!;
  const a = conAjuste.estudios.find((e) => e.id === 1)!, b = sinAjuste.estudios.find((e) => e.id === 1)!;
  ok(cerca(a.se, b.se * Math.sqrt(deff), 1e-12), "El motor multiplica el EE del estudio clúster por √DEFF (con DEFF exacto 1,5752)");
  ok(cerca(a.y, b.y, 1e-12), "El efecto (ln OR) del estudio clúster no cambia: solo cambia su precisión");
  ok(a.se > b.se && conAjuste.fijo.icSuperior - conAjuste.fijo.icInferior > sinAjuste.fijo.icSuperior - sinAjuste.fijo.icInferior,
    "Ajustar por clúster ensancha el IC combinado: sin ajuste el estudio recibiría más peso del debido (unit-of-analysis error)");
  const incompleto = calcularMetaAnalisis([{ ...cluster, varianteEnsayo: { tipo: "cluster", numClusters: "", icc: "0.02" } }, normal, est(3, 1.1, 0.7, 1.7)], "OR", "DL")!;
  ok(!incompleto.estudios.some((e) => e.id === 1), "Clúster con datos incompletos → el estudio se EXCLUYE (no se usa sin ajustar)");
  const iccInvalido = calcularMetaAnalisis([{ ...cluster, varianteEnsayo: { tipo: "cluster", numClusters: "21", icc: "1.5" } }, normal, est(3, 1.1, 0.7, 1.7)], "OR", "DL")!;
  ok(!iccInvalido.estudios.some((e) => e.id === 1), "Clúster con ICC > 1 → el estudio se EXCLUYE");
}

// ------------------------------------------------------------------------------------------------
// 23.2 Cruzados — MD = ME − MC ; SE(MD) = SDdiff / √N ; SDdiff = √(SDE² + SDC² − 2·Corr·SDE·SDC)  (23.2.7.1)
// ------------------------------------------------------------------------------------------------
// Ejemplo 23.2.7.4: N = 10; E: media 7,0 y DE 2,38; C: media 6,5 y DE 2,21.
{
  const conSd = efectoDesdeCruzado(7.0, 6.5, 10, 2, NaN, NaN, NaN)!;
  ok(cerca(conSd.valor, 0.5, 1e-12), "MD = ME − MC = 0,5");
  ok(cerca(conSd.se, 0.632, 5e-4), `SDdiff imputada = 2 → EE = 2/√10 = ${conSd.se.toFixed(4)}; Handbook: 0,632`);

  const conCorr = efectoDesdeCruzado(7.0, 6.5, 10, NaN, 2.38, 2.21, 0.68)!;
  const sdDiff = Math.sqrt(2.38 ** 2 + 2.21 ** 2 - 2 * 0.68 * 2.38 * 2.21);
  ok(cerca(sdDiff, 1.8426, 5e-5), `SDdiff con Corr = 0,68 = ${sdDiff.toFixed(4)}; Handbook: 1,8426 (el libro imprime «1,846» una línea antes: errata, el cálculo posterior usa 1,8426)`);
  ok(cerca(conCorr.se, 0.583, 5e-4), `EE con Corr = 0,68 = ${conCorr.se.toFixed(4)}; Handbook: 0,583`);
  ok(cerca(conCorr.se, sdDiff / Math.sqrt(10), 1e-12), "EE = SDdiff/√N");
  ok(cerca(conCorr.icInferior, 0.5 - 1.959963985 * conCorr.se, 1e-9), "IC95 = MD ± 1,96·EE");

  // Corr = (SDE² + SDC² − SDdiff²) / (2·SDE·SDC) (23.2.7.3): la inversa de la fórmula anterior
  const corr = (2.38 ** 2 + 2.21 ** 2 - sdDiff ** 2) / (2 * 2.38 * 2.21);
  ok(cerca(corr, 0.68, 1e-12), "La correlación calculada desde las tres DE (23.2.7.3) devuelve 0,68: las fórmulas son inversas");

  // La DE de las diferencias reportada tiene prioridad sobre la imputada
  const ambas = efectoDesdeCruzado(7.0, 6.5, 10, 2, 2.38, 2.21, 0.68)!;
  ok(cerca(ambas.se, 2 / Math.sqrt(10), 1e-12), "Si el estudio reporta SDdiff, se usa esa y no la imputada");

  // Sin información suficiente no se inventa nada
  ok(efectoDesdeCruzado(7, 6.5, 10, NaN, 2.38, 2.21, NaN) === null, "Sin SDdiff y sin correlación → null (no se inventa una correlación)");
  ok(efectoDesdeCruzado(7, 6.5, 10, NaN, 2.38, 2.21, 1.2) === null, "Correlación fuera de −1 a 1 → null");
  ok(efectoDesdeCruzado(7, 6.5, 1, 2, NaN, NaN, NaN) === null, "N < 2 → null");
  ok(efectoDesdeCruzado(7, 6.5, 10, NaN, 2, 2, 1) === null, "Corr = 1 con DE iguales (SDdiff = 0) → null en vez de un EE de cero");

  // Por qué no usar el análisis de grupos independientes (23.2.6): con Corr > 0 sobreestima el EE (IC demasiado ancho, estudio infraponderado)
  const independiente = efectoDesdeMedias(7.0, 2.38, 10, 6.5, 2.21, 10)!;
  ok(independiente.se > conCorr.se, `Con Corr = 0,68 el análisis de grupos independientes (EE ${independiente.se.toFixed(3)}) es más ancho que el pareado (${conCorr.se.toFixed(3)}), como dice 23.2.6`);
  const corr0 = efectoDesdeCruzado(7.0, 6.5, 10, NaN, 2.38, 2.21, 0)!;
  ok(cerca(corr0.se, independiente.se, 1e-12), "Con Corr = 0 el pareado coincide con el de grupos independientes (sin ventaja del diseño cruzado, 23.2.7.1)");
}

// Cruzado con desenlace distinto de DM: no se calcula (Handbook 23.2.5: binarios requieren un estadístico)
{
  // `calcularEfectoCrudos` solo calcula cruzado para DM; ver efecto-crudos.ts
  const { calcularEfectoCrudos } = await import("../src/efecto-crudos.ts");
  const crudos = { tipo: "cruzado" as const, mediaTrat: "7", mediaCtrl: "6.5", n: "10", sdDiff: "2", sdTrat: "", sdCtrl: "", correlacion: "" };
  ok(calcularEfectoCrudos("DM", crudos) !== null, "Cruzado + DM: se calcula");
  ok(calcularEfectoCrudos("OR", crudos) === null && calcularEfectoCrudos("DME", crudos) === null, "Cruzado + OR/DME: no se calcula (no se adivina)");
}

console.log(fallos === 0 ? "\nTodo bien." : `\n${fallos} fallo(s).`);
process.exit(fallos === 0 ? 0 : 1);
