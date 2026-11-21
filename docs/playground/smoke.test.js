import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { test } from "node:test";

const bytes = await readFile(new URL("../site/assets/playground/meshers.wasm", import.meta.url));
const { instance: { exports: core } } = await WebAssembly.instantiate(bytes);

function generate(field, cells, width) {
  const pointer = core.generate(field, cells, width);
  const buffer = new Uint8Array(core.memory.buffer, pointer, core.response_length());
  return JSON.parse(new TextDecoder().decode(buffer));
}

test("the published module meshes both fields and returns wall and cap triangles", () => {
  for (const field of [0, 1]) {
    const coarse = generate(field, 6, 0.3);
    const fine = generate(field, 16, 0.8);
    assert.ok(fine.tets > coarse.tets && coarse.tets > 0);
    for (const mesh of [coarse, fine]) {
      assert.equal(mesh.faces.length, mesh.tags.length);
      assert.equal(new Set(mesh.tags).size, 7);
      assert.ok(mesh.faces.every(face => face.length === 3 &&
        face.every(node => Number.isInteger(node) && node >= 0 && node < mesh.points.length)));
    }
  }
  assert.ok(generate(0, 1000, 0.5).error);
});

test("the browser mesher cannot grow its linear memory beyond 32 MiB", () => {
  const pages = core.memory.buffer.byteLength / 65536;
  assert.throws(() => core.memory.grow(513 - pages), RangeError);
});

test("optimization improves the coarse gyroid shown in the explorer", () => {
  const mesh = generate(0, 8, 0.5);
  assert.ok(mesh.minQuality > 0.1); // Previously 0.041 with optimization disabled.
  const angles = mesh.faces.flatMap(face => face.map((node, i) => {
    const edge = end => mesh.points[end].map((x, axis) => x - mesh.points[node][axis]);
    const a = edge(face[(i + 1) % 3]);
    const b = edge(face[(i + 2) % 3]);
    const dot = a.reduce((sum, x, axis) => sum + x * b[axis], 0);
    return Math.acos(dot / (Math.hypot(...a) * Math.hypot(...b))) * 180 / Math.PI;
  }));
  assert.ok(Math.abs(mesh.minAngle - Math.min(...angles)) < 1e-6);
});
