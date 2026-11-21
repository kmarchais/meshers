# Rust API

<ol class="mesh-flow">
  <li><strong>ScalarField</strong>Evaluate the implicit function</li>
  <li><strong>Options</strong>Set the region and sampling grid</li>
  <li><strong>Output</strong>Inspect the mesh and diagnostics</li>
</ol>

[Browse the Rust crate reference](../rust/meshers_core/index.html).

Pass a closure to `implicit::generate` for numerical gradients. Implement `implicit::ScalarField` to supply gradients or batched evaluation.

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

The CPU core runs independently of Python, VTK and CAD runtimes.
