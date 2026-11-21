import { boundaryMesh } from "./boundary-mesh.mjs";

const ready = fetch("playground/meshers.wasm")
  .then(response => {
    if (!response.ok) throw new Error("Could not load the mesher.");
    return response.arrayBuffer();
  })
  .then(bytes => WebAssembly.instantiate(bytes))
  .then(({ instance }) => instance.exports);

self.onmessage = async ({ data: { field, cells, width } }) => {
  try {
    const core = await ready;
    self.postMessage({ ready: true });
    const address = core.generate(field, cells, width);
    const bytes = new Uint8Array(core.memory.buffer, address, core.response_length());
    const mesh = JSON.parse(new TextDecoder().decode(bytes));
    if (mesh.error) throw new Error(mesh.error);
    const data = boundaryMesh(mesh);
    self.postMessage(data, [data.coords.buffer, data.indices.buffer]);
  } catch (error) {
    self.postMessage({ error: error.message });
  }
};
