// Run against a saved baseline or current WASM; includes adapter diagnostics.
import { readFile, writeFile } from "node:fs/promises";
import { createMeshers } from "./meshers.js";
import assert from "node:assert/strict";
const engine = await createMeshers(await readFile(process.argv[2]));
const records = [];
for (const [name, shape, preset, repeat, resolution, grade] of [
  ["gyroid-fea", 0, 2, 1, 16, 0],
  ["split-p-fea", 1, 2, 1, 16, 0],
  ["graded-gyroid-fea", 0, 2, 2, 12, 0.3],
  ["gyroid-fast", 0, 0, 1, 16, 0],
  ["split-p-fast", 1, 0, 1, 16, 0],
  ["graded-gyroid-fast", 0, 0, 3, 20, 0.3],
  ["graded-split-p-fast", 1, 0, 3, 20, 0.3],
]) {
  const times = [];
  let mesh;
  for (let run = 0; run < 4; run++) {
    const start = performance.now();
    mesh = engine.generate({
      shape,
      preset,
      repeat,
      resolution,
      grade,
      thickness: 0.6,
    });
    assert.equal(mesh.metrics.bad_edges, 0);
    assert.deepEqual(mesh.metrics.periodic_matches, [
      grade ? null : true,
      true,
      true,
    ]);
    if (preset === 2) assert(mesh.metrics.volume.target_met);
    if (run) times.push((performance.now() - start) / 1000);
  }
  const record = {
    name,
    seconds: times.sort((a, b) => a - b)[1],
    times,
    triangles: mesh.triangles.length,
    tetrahedra: mesh.tetrahedra.length,
    metrics: mesh.metrics,
  };
  records.push(record);
  console.log(JSON.stringify(record));
}
await writeFile(
  process.argv[3],
  JSON.stringify({ runtime: process.version, records }, null, 2),
);
