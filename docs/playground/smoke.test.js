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
