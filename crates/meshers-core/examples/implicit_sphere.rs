//! Minimal custom implicit field; the documentation links to this example.
use meshers_core::implicit::{Options, generate};

fn main() -> Result<(), meshers_core::MeshingError> {
    let field = |p: [f64; 3]| {
        (p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2) + (p[2] - 0.5).powi(2) - 0.3_f64.powi(2)
    };
    let output = generate(&field, Options::default())?;
    println!("{} tetrahedra", output.mesh.tets.len());
    Ok(())
}
