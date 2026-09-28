import assert from "node:assert/strict";
import { readFile, writeFile } from "node:fs/promises";
import { createMeshers } from "./dist/package/index.js";
const engine = await createMeshers(
  await readFile("./dist/package/meshers.wasm"),
);
const records = [];
for (const shape of [0, 1, 2])
  for (const preset of [0, 1, 2]) {
    const start = performance.now();
    const mesh = engine.generate({
      shape,
      preset,
      resolution: [12, 20, 16][preset],
      repeat: 1,
      thickness: 0.6,
      grade: 0,
    });
    assert(mesh.triangles.length > 0);
    assert.equal(mesh.metrics.bad_edges, 0);
    assert.deepEqual(mesh.metrics.periodic_matches, [true, true, true]);
    assert.equal(mesh.tetrahedra.length > 0, preset === 2);
    if (preset === 2) assert(mesh.metrics.volume.target_met);
    const result = {
      shape,
      preset,
      seconds: (performance.now() - start) / 1000,
      triangles: mesh.triangles.length,
      tetrahedra: mesh.tetrahedra.length,
      metrics: mesh.metrics,
    };
    records.push(result);
    console.log(
      shape,
      preset,
      result.seconds.toFixed(3),
      result.triangles,
      result.tetrahedra,
    );
  }
for (const preset of [0, 1, 2]) {
  const mesh = engine.generate({
    shape: 0,
    preset,
    repeat: 2,
    resolution: 12,
    thickness: 0.6,
    grade: 0.3,
  });
  assert.equal(mesh.metrics.bad_edges, 0);
  assert.deepEqual(mesh.metrics.periodic_matches, [null, true, true]);
  records.push({
    case: "graded-2-cube",
    preset,
    triangles: mesh.triangles.length,
    tetrahedra: mesh.tetrahedra.length,
    metrics: mesh.metrics,
  });
}
assert.throws(() => engine.generate({ thickness: NaN }), /Choose/);
assert.throws(
  () => engine.generate({ preset: 2, repeat: 3, resolution: 32 }),
  /budget/,
);
const again = engine.generate();
assert(again.triangles.length > 0);
await writeFile(
  "smoke-results.json",
  JSON.stringify({ runtime: process.version, records }, null, 2),
);
console.log(
  "All WASM smoke checks passed, including graded domains and recovery after errors.",
);
