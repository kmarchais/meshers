import { readFile, writeFile } from "node:fs/promises";
import assert from "node:assert/strict";
const variants = {};
for (const name of ["baseline", "analytic", "active", "analytic-active", "baseline-repeat"]) {
  const log = await readFile(`dist/experiment-${name}.log`, "utf8");
  variants[name] = [...log.matchAll(/EXPERIMENT (\{[^\r\n]+\})/g)].map(m => JSON.parse(m[1]));
}
for (const before of variants["baseline-repeat"]) {
  const active = variants.active.find(r => r.case === before.case);
  if (before.error) { assert(active.error); continue; }
  assert.equal(active.fingerprint, before.fingerprint);
  for (const name of ["analytic", "analytic-active"]) {
    const after = variants[name].find(r => r.case === before.case);
    assert(!after.error && after.metrics.volume.target_met);
    assert.equal(after.tets, before.tets);
    assert.equal(after.metrics.bad_edges, 0);
    assert.deepEqual(after.metrics.periodic_matches, before.metrics.periodic_matches);
  }
  assert.equal(variants.analytic.find(r=>r.case===before.case).fingerprint,
    variants["analytic-active"].find(r=>r.case===before.case).fingerprint);
}
const gil = JSON.parse(await readFile("dist/gil-results.json", "utf8"));
await writeFile("optimizer-experiments-results.json", JSON.stringify({
  baseline_commit:"d4f276f",
  method:"Native release builds; three runs per case, medians. Variants run sequentially. Initial baseline had higher load and stopped on the coarse graded Split-P error; the complete repeated baseline is the reference. Timing differences of a few percent are inconclusive. GIL benchmark uses the rebuilt default Python extension, not optimizer experiments. Compiled and callback throughput probes use different resolutions and pass counts; only compare concurrency within each probe.",
  variants, gil,
}, null, 2));
console.log("Active meshes match baseline; combined meshes match analytic; successful analytic cases retain quality targets, counts, closure and periodic matching.");
