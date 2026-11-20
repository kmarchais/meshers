# Rust API

[Browse the Rust crate reference](../rust/meshers_core/index.html).

Start with `implicit::generate`, `implicit::Options`, `implicit::Region` and
`implicit::ScalarField`. A closure implements a field with numerical gradients;
implement `ScalarField` to supply gradients or batched evaluation.

```rust
use meshers_core::implicit::{generate, Options};

fn main() -> Result<(), meshers_core::MeshingError> {
    let field = |p: [f64; 3]| {
        (p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2)
            + (p[2] - 0.5).powi(2) - 0.3_f64.powi(2)
    };
    let output = generate(&field, Options::default())?;
    println!("{} tetrahedra", output.mesh.tets.len());
    Ok(())
}
```

The CPU core has no Python, VTK or CAD runtime dependency. The Python extension
owns field compilation and language bindings; export belongs to adapters.
