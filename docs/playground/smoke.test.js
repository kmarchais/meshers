import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { test } from "node:test";
import { qualityFigure } from "../site/assets/guide-figures.mjs";
import { cutFigure, foldFigure, sagFigure } from "../site/assets/course-figures.mjs";
import { triangleSection } from "../site/assets/extraction-figures.mjs";
import { boundaryMesh } from "../site/assets/boundary-mesh.mjs";

const bytes = await readFile(new URL("../site/assets/playground/meshers.wasm", import.meta.url));

test("course controls reproduce the worked examples and the fold threshold", () => {
  assert.equal(cutFigure(50).value, "0.50 · 12.5% volume");
  assert.equal(sagFigure(60).value, "h = 1.20 · sag = 0.200");
  assert.equal(sagFigure(30).value, "h = 0.60 · sag = 0.046");
  assert.equal(foldFigure(24).value, "min dx/dt = 0.04 · increasing");
  assert.equal(foldFigure(25).value, "min dx/dt = 0.00 · singular endpoint");
  assert.equal(foldFigure(26).value, "min dx/dt = -0.04 · folded");
});
const { instance: { exports: core } } = await WebAssembly.instantiate(bytes);

test("surface lesson roots lie on the cutting plane and form the expected contour", () => {
  for (const [level, expected] of [[0.5, 3], [1.5, 4], [3.5, 0]]) {
    const roots = triangleSection(level);
    assert.equal(roots.length, expected);
    for (const [x, y, z] of roots) assert.ok(Math.abs(x + y + z - level) < 1e-12);
  }
  assert.deepEqual(triangleSection(1.5), [
    [0.75, 0.75, 0], [0.5, 0.5, 0.5], [1, 0.25, 0.25], [1, 0.5, 0],
  ]);
});

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

test("the starter has the same boundary, labels and metrics as the live default mesh", async () => {
  const url = new URL("../site/assets/playground/gyroid.json", import.meta.url);
  const starter = JSON.parse(await readFile(url, "utf8"));
  const { coords, indices, ...metrics } = boundaryMesh(generate(0, 10, 0.5));
  const { coords: initial, indices: faces, ...reported } = starter;
  assert.deepEqual(reported, metrics);
  assert.deepEqual(faces, Array.from(indices));
  assert.equal(initial.length, coords.length);
  assert.ok(coords.every((x, i) => Math.abs(x - initial[i]) <= 5.1e-8));
});

test("the quality illustration agrees with the volume and six edges of a Tet4", () => {
  for (const percent of [1, 5, 50, 100]) {
    const height = Math.sqrt(2 / 3) * percent / 100;
    const volume = Math.sqrt(3) / 4 * height / 3;
    const sumEdgeSquared = 3 + 3 * (1 / 3 + height ** 2);
    const mmg = 72 * Math.sqrt(3) * volume / sumEdgeSquared ** 1.5;
    assert.ok(qualityFigure(percent).graphic.includes("Adjusted · " + mmg.toFixed(3)));
  }
});
