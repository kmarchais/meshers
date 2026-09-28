# TPMS studio

Experimental local browser demo of meshers. The standalone HTML bundles the Rust WASM engine, Three.js, controls and worker code. It needs no server or network after download. Generation runs in a dedicated worker; Stop terminates that worker.

Quick view raymarches the implicit TPMS field in a GPU fragment shader, without generating a mesh or loading the WASM engine. Geometry, repeats, offset, grading and cut update live. The frame-rate display includes browser frame pacing and is not a mesh-generation benchmark. Finite ray steps and pixel tolerance can miss thin features. Mesh checks and exports are disabled in this view. Printing and FEA generate real meshes; the comparison button runs those two modes. The library's numeric preset 0 still generates a linear surface for callers that need it.

Build from the meshers repository:

```powershell
rustup target add wasm32-unknown-unknown
cargo build -p meshers-wasm --release --target wasm32-unknown-unknown
cd demos/tpms-browser
npm ci
npm run build
npm test
```

Open `index.html` and follow **Open TPMS studio**, or open `dist/index.html` directly. The reusable ES module package is in `dist/package`. `studio.template.html` is a build input and does not run on its own.

```js
import { createMeshers } from './package/index.js';
const meshers = await createMeshers();
const mesh = meshers.generate({shape: 0, preset: 0, repeat: 1,
  resolution: 12, thickness: 0.6, grade: 0});
```

Shapes: 0 gyroid, 1 split-P, 2 Schwarz P. Presets: 0 linear surface, 1 exact-intersection surface, 2 optimized tetrahedral volume. The API is synchronous; call it inside a worker. Import from a server when loading separate WASM files. Output coordinates are in unit-cell lengths. The standalone UI scales exports to millimetres using its cell-size control.

The presets default to 12, 20 and 16 samples per cell. For repeat `r`, the grid has `resolution*r-1` divisions per axis. Grading varies the implicit sheet offset from `thickness*(1-grade)` to `thickness*(1+grade)` along X. It disables X periodicity and keeps Y/Z periodicity. Offset is not physical wall thickness or relative density.

Printing uses exact roots without triangle polishing. FEA uses the real volume mesher, four optimization passes before and after adaptive refinement, and sampled geometry tolerance 0.01 unit cells. The demo reports candidates below minimum MMG quality 0.1 as target misses. A preset is not a printability or solver certification. Surface deviation is a first-order distance estimate at face centroids, not a maximum-error bound. Surface closure checks oriented edges, not self-intersection or vertex manifoldness. The volume core performs additional topology and element validity checks.

Browser limits: surfaces at most 64 grid divisions per axis, volumes at most 40 and 600,000 tetrahedra. The worker uses serial Rust execution. Timings include generation, validation and JSON serialization/parsing, but exclude WASM loading, worker transfer and rendering. The package is experimental and not published to npm.

Quick view uses four subpixel rays per pixel for antialiasing, including partial coverage at implicit silhouettes. This costs more GPU work than a single ray; frame rate depends on viewport size, geometry and hardware.

Three.js is MIT licensed. Its license is distributed with the built artifacts.

See [PERFORMANCE.md](PERFORMANCE.md) for optimizer and fast-extractor measurements, quality tradeoffs and reproduction commands.
