import { createMeshers } from "./meshers.js";
const engine = createMeshers(
  Uint8Array.from(atob(WASM_BASE64), (c) => c.charCodeAt(0)),
);
self.onmessage = async ({ data }) => {
  try {
    const api = await engine;
    const start = performance.now();
    const mesh = api.generate(data);
    self.postMessage({
      ok: true,
      mesh,
      seconds: (performance.now() - start) / 1000,
      config: data,
    });
  } catch (error) {
    self.postMessage({ ok: false, error: error.message || String(error) });
  }
};
