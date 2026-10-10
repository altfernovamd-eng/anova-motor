// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Fernando Ramírez Nova
/**
 * Prueba de regresión: valida la concordancia inter-observador de aNova
 * (`src/lib/meta-analysis/kappa.ts`) contra el paquete `irr` de R — el estándar de
 * referencia para fiabilidad entre evaluadores — reproduciendo su salida publicada.
 *
 * Dataset: `diagnoses` (30 sujetos, 6 evaluadores, 5 categorías) del paquete `irr`.
 * Referencia:
 *   - Cohen  : kappa2(diagnoses[,2:3])  → κ = 0.631, z = 7.56   (2 evaluadores)
 *   - Fleiss : kappam.fleiss(diagnoses) → κ = 0.43,  z = 17.7   (6 evaluadores)
 *
 * Correrlo:  npm run benchmark:kappa   (o)   node pruebas/benchmark-kappa.mts
 * Sale con código 1 si alguna cifra deja de coincidir.
 */
import { cohenKappa, fleissKappa } from "../src/kappa.ts";

// diagnoses codificado por el dígito inicial de cada categoría (1..5); una cadena = 6 calificaciones del sujeto.
const DIAG = [
  "444444", "222555", "233335", "555555", "222444", "113333", "333355", "113334", "114444", "555555",
  "144444", "124444", "222333", "144444", "224445", "333335", "111455", "111112", "224444", "133555",
  "555555", "244444", "224555", "114444", "144445", "222224", "111155", "224444", "133333", "555555",
];

let fallos = 0;
let total = 0;
const check = (nombre: string, anova: number, meta: number, dp: number) => {
  const a = anova.toFixed(dp);
  const m = meta.toFixed(dp);
  const ok = a === m;
  if (!ok) fallos++;
  total++;
  console.log(`${nombre.padEnd(26)} ${a.padEnd(12)} ${m.padEnd(12)} ${ok ? "✓" : "✗ DIFIERE"}`);
};

// --- Cohen (κ2), columnas 2 y 3 (índices 1 y 2) ---
const cohen = cohenKappa(DIAG.map((s) => [s[1], s[2]] as [string, string]))!;
console.log(`\naNova  vs  R 'irr' — κ de Cohen (kappa2, 2 evaluadores)  —  diagnoses[,2:3]\n`);
console.log(`${"Cantidad".padEnd(26)} ${"aNova".padEnd(12)} ${"irr".padEnd(12)} Estado`);
console.log("-".repeat(62));
check("kappa", cohen.kappa, 0.631, 3);
check("z", cohen.z, 7.56, 2);
console.log("-".repeat(62));
console.log(`(p-valor aNova = ${cohen.pValor.toExponential(2)}; irr ≈ 4.04e-14 — ambos ≈ 0)`);

// --- Fleiss (κ), los 6 evaluadores ---
const fleiss = fleissKappa(DIAG.map((s) => s.split("")))!;
console.log(`\naNova  vs  R 'irr' — κ de Fleiss (kappam.fleiss, 6 evaluadores)  —  diagnoses\n`);
console.log(`${"Cantidad".padEnd(26)} ${"aNova".padEnd(12)} ${"irr".padEnd(12)} Estado`);
console.log("-".repeat(62));
check("kappa", fleiss.kappa, 0.43, 2);
check("z", fleiss.z, 17.7, 1);
console.log("-".repeat(62));

if (fallos === 0) {
  console.log(`\n✓ PASA: las ${total} cantidades coinciden con R (irr) a los decimales publicados.\n`);
  process.exit(0);
} else {
  console.error(`\n✗ FALLA: ${fallos} de ${total} cantidad(es) ya no coinciden con irr.\n`);
  process.exit(1);
}
