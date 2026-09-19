use meshers_core::{GyroidOptions, MeshingError, generate_gyroid};
fn main() -> Result<(), MeshingError> {
    let mesh = generate_gyroid(GyroidOptions {
        resolution: 24,
        ..GyroidOptions::default()
    })?;
    println!(
        "{} points, {} tetrahedra",
        mesh.points.len(),
        mesh.tets.len()
    );
    Ok(())
}
