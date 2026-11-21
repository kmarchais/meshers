// Generate the initial view with the same binary and boundary packing as the worker.
import { readFile, writeFile } from "node:fs/promises";
import { boundaryMesh } from "../site/assets/boundary-mesh.mjs";

const assets = new URL("../site/assets/playground/", import.meta.url);
const bytes = await readFile(new URL("meshers.wasm", assets));
const { instance: { exports: core } } = await WebAssembly.instantiate(bytes);
const address = core.generate(0, 10, 0.5);
const raw = new Uint8Array(core.memory.buffer, address, core.response_length());
const mesh = JSON.parse(new TextDecoder().decode(raw));
if (mesh.error) throw new Error(mesh.error);
const data = boundaryMesh(mesh);
await writeFile(new URL("gyroid.json", assets), JSON.stringify({
  ...data,
  coords: Array.from(data.coords, x => Number(x.toFixed(7))),
  indices: Array.from(data.indices),
}));
