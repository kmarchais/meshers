# Meshers CPU core

Use `implicit::generate` for a bounded scalar field, with `implicit::Options` to set per-axis periodicity. The crate returns float64 points, tetrahedra and oriented boundary triangles. [Rust API and example](https://kmarchais.github.io/meshers/api/rust/)

The specialized `generate_gyroid` API retains its original unit-box behavior. Its `threshold` is the field half-width. `periodic=false` disables the final pairing requirement but still groups translated nodes. Use `implicit::Options` for nonperiodic geometry.

## Run from the repository root

```sh
cargo test --release
cargo clippy --all-targets -- -D warnings
cargo run --release --example gyroid
```

The core forbids unsafe code. It runs on the CPU without Python or a GPU driver. The experimental `accelerator` module defines an interface for external optimizers.
