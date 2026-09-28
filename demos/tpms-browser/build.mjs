import { build } from "esbuild";
import { readFile, writeFile, mkdir, copyFile } from "node:fs/promises";
const wasm = await readFile(
  "../../target/wasm32-unknown-unknown/release/meshers_wasm.wasm",
);
const worker = await build({
  entryPoints: ["worker.js"],
  bundle: true,
  write: false,
  format: "iife",
  minify: true,
  logOverride: { "empty-import-meta": "silent" },
  define: { WASM_BASE64: JSON.stringify(wasm.toString("base64")) },
});
const app = await build({
  entryPoints: ["app.js"],
  bundle: true,
  write: false,
  format: "iife",
  minify: true,
  define: { WORKER_SOURCE: JSON.stringify(worker.outputFiles[0].text) },
});
await mkdir("dist/package", { recursive: true });
await writeFile("dist/worker.bundle.js", worker.outputFiles[0].text);
const template = await readFile("index.html", "utf8");
const thirdParty = await readFile("node_modules/three/LICENSE", "utf8");
await writeFile(
  "dist/index.html",
  template.replace("/* APP_BUNDLE */", () =>
    app.outputFiles[0].text.replaceAll("</script", "<\\/script"),
  ) + `\n<!-- Three.js license\n${thirdParty}\n-->`,
);
await writeFile("dist/THIRD_PARTY_NOTICES.txt", thirdParty);
await writeFile("dist/package/meshers.wasm", wasm);
await copyFile("meshers.js", "dist/package/index.js");
await writeFile(
  "dist/package/package.json",
  JSON.stringify(
    {
      name: "@meshers/tpms-wasm-experiment",
      version: "0.1.0-experiment",
      type: "module",
      exports: "./index.js",
      files: ["index.js", "meshers.wasm"],
      license: "MIT",
    },
    null,
    2,
  ),
);
await copyFile("../../LICENSE", "dist/package/LICENSE");
await copyFile("README.md", "dist/package/README.md");
console.log(
  `Built standalone dist/index.html and WASM package (${(wasm.length / 1048576).toFixed(2)} MiB wasm).`,
);
