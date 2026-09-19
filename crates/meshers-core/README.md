# Meshers CPU core

CPU float64 tetrahedral meshing of implicit solids and constrained intersections.
General bounded fields are available through `implicit::generate`; see
[Rust API guide](https://github.com/kmarchais/meshers/blob/main/docs/site/api/rust.md).

```rust
use meshers_core::{generate_gyroid, GyroidOptions};

let mesh = generate_gyroid(GyroidOptions {
    resolution: 24,
    ..GyroidOptions::default()
})?;
let diagnostics = mesh.metrics(0.5);
assert!(diagnostics.periodic_nodes_and_triangles.iter().all(|&v| v));
# Ok::<(), meshers_core::MeshingError>(())
```

- `Mesh` owns float64 points, tetrahedron indices and oriented boundary triangles.
- `GyroidOptions` configures the current unit-box band generator. `threshold`
  is the raw field half-width. `threads=0` selects the legacy serial traversal;
  positive thread counts use colored parallel optimization. The sequences can
  differ, so compare like-for-like settings.
- `MeshingError` distinguishes invalid options from generation failure. Worker
  pool creation failures return errors rather than panicking.
- `metrics(threshold)` includes geometry-specific gyroid diagnostics; this is
  not yet a general imported-mesh quality API.
- On the specialized `GyroidOptions` API, `periodic=false` currently disables the final requirement only. The prototype
  still groups translated nodes. Use `implicit::Options` for true per-axis or nonperiodic geometry. The
  specialized entry point keeps its historical behavior for reproducibility.
- The experimental `accelerator` module is a narrow adapter boundary for the
  external optimizers. Its factory and engine contain no device types. Device
  acceptance is a trusted numerical extension and must validate its own moves.

The crate depends only on Rayon, Serde and serde_json and their Rust dependencies.
It forbids unsafe code. It requires no Python, OpenCascade, HDF5, CUDA or GPU
driver. No device backend is included in this crate.

From the repository root:

```sh
cargo test --release
cargo clippy --all-targets -- -D warnings
cargo run --release --example gyroid
```

The earlier 2D CAD crate remains an independent workspace at `crates/meshers`.
Use `cargo test --manifest-path crates/meshers/Cargo.toml` for that crate. Keeping
it separate also keeps its git/OpenCascade dependency out of CPU resolution.
The legacy CAD CI job remains separate from the new CPU platform matrix.
