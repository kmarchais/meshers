use meshers_core::{
    Point,
    implicit::{Options, Region, generate},
};
use std::f64::consts::TAU;
fn gyroid(p: Point) -> f64 {
    let [x, y, z] = p.map(|v| v * TAU);
    x.sin() * y.cos() + y.sin() * z.cos() + z.sin() * x.cos()
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let name = std::env::args().nth(1).unwrap_or_else(|| "graded".into());
    let field = |p: Point| match name.as_str() {
        "sphere" => {
            ((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2) + (p[2] - 0.5).powi(2)).sqrt() - 0.3
        }
        "band" => gyroid(p),
        _ => {
            let w = p[0];
            let schwarz = p.iter().map(|v| (TAU * v).cos()).sum::<f64>();
            ((1. - w) * gyroid(p) + w * schwarz).abs() - (0.5 + 0.2 * p[2])
        }
    };
    if !["sphere", "band", "graded"].contains(&name.as_str()) {
        return Err("choose sphere, band or graded".into());
    }
    let options = Options {
        cells: [28; 3],
        geometry_tolerance: 0.04,
        region: if name == "band" {
            Region::Band {
                lower: -0.5,
                upper: 0.5,
            }
        } else {
            Region::Negative
        },
        periodic: [name == "band"; 3],
        ..Options::default()
    };
    let start = std::time::Instant::now();
    let result = generate(&field, options)?;
    println!(
        "{}",
        serde_json::json!({"case":name,"seconds":start.elapsed().as_secs_f64(),
        "points":result.mesh.points.len(),"tetrahedra":result.mesh.tets.len(),
        "volume":result.diagnostics.volume,"minimum_mmg_quality":result.diagnostics.minimum_mmg_quality,
        "sampled_surface_error":result.diagnostics.maximum_sampled_surface_error,
        "periodic_node_pairs":result.periodic_pairs.each_ref().map(|pairs|pairs.len())})
    );
    Ok(())
}
