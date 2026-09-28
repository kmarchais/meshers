import assert from "node:assert/strict";
import { readFile, writeFile } from "node:fs/promises";
import { Worker } from "node:worker_threads";
import { Script } from "node:vm";
import { stlBytes, vtuText } from "./exports.js";
const html = await readFile("dist/index.html", "utf8");
const script = html.slice(
  html.indexOf("<script>") + 8,
  html.lastIndexOf("</script>"),
);
new Script(script); // Verify the emitted inline bundle parses, not just its source.
assert(!html.includes("/* APP_BUNDLE */"));
const bundle = await readFile("dist/worker.bundle.js", "utf8");
function start() {
  return new Worker(
    `const {parentPort}=require('node:worker_threads');globalThis.self={postMessage:m=>parentPort.postMessage(m)};${bundle};parentPort.on('message',data=>self.onmessage({data}));`,
    { eval: true },
  );
}
function request(worker, data) {
  return new Promise((resolve, reject) => {
    worker.once("message", resolve);
    worker.once("error", reject);
    worker.postMessage(data);
  });
}
const worker = start();
const result = await request(worker, {
  shape: 0,
  preset: 2,
  repeat: 1,
  resolution: 16,
  thickness: 0.6,
  grade: 0,
});
assert(result.ok);
assert(result.mesh.metrics.volume.target_met);
const stl = stlBytes(result.mesh, 10);
assert.equal(
  new DataView(stl).getUint32(80, true),
  result.mesh.triangles.length,
);
assert.equal(stl.byteLength, 84 + 50 * result.mesh.triangles.length);
await writeFile("dist/check.stl", new Uint8Array(stl));
await writeFile("dist/check.vtu", vtuText(result.mesh, 10));
const failure = await request(worker, {
  shape: 0,
  preset: 2,
  repeat: 3,
  resolution: 32,
  thickness: 0.6,
  grade: 0,
});
assert.equal(failure.ok, false);
worker.postMessage({
  shape: 1,
  preset: 2,
  repeat: 2,
  resolution: 20,
  thickness: 0.6,
  grade: 0.3,
});
await worker.terminate();
const replacement = start();
const preview = await request(replacement, {
  shape: 0,
  preset: 0,
  repeat: 1,
  resolution: 12,
  thickness: 0.6,
  grade: 0,
});
assert(preview.ok);
await replacement.terminate();
console.log(
  "Inline bundle parses. Embedded worker generates volumes, returns errors, terminates and restarts. Actual UI export functions produced STL and VTU fixtures.",
);
