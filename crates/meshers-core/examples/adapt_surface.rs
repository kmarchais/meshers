//! Adaptive Gyroid sheet surface with per-phase timing.
//! cargo run --release --features experimental-quads --example adapt_surface -- 1e-3 32 8
use meshers_core::{Point, quads::Band, remesh, triangles};

fn main() {
    let mut args = std::env::args().skip(1);
    let tolerance: f64 = args.next().map_or(1e-3, |v| v.parse().unwrap());
    let cells: usize = args.next().map_or(32, |v| v.parse().unwrap());
    let iterations: usize = args.next().map_or(8, |v| v.parse().unwrap());
    let field = |p: Point| {
        let t = std::f64::consts::TAU;
        ((t * p[0]).sin() * (t * p[1]).cos()
            + (t * p[1]).sin() * (t * p[2]).cos()
            + (t * p[2]).sin() * (t * p[0]).cos())
            / 0.55
    };
    let band = Band {
        field: &field,
        bounds: [[-0.5; 3], [0.5; 3]],
        levels: [-1., 1.],
    };
    let start = std::time::Instant::now();
    let mut mesh = triangles::extract(&band, cells).unwrap();
    let extraction = start.elapsed().as_secs_f64();
    let extracted = mesh.faces.len();
    let options = remesh::Options {
        tolerance,
        min_length: tolerance / 2.,
        max_length: 0.25,
        gradation: 0.3,
        iterations,
        balanced: std::env::args().any(|a| a == "balanced"),
        fast: std::env::args().any(|a| a == "fast"),
    };
    let stats = remesh::adapt(&mut mesh, &band, &options).unwrap();
    let total = start.elapsed().as_secs_f64();
    println!(
        "extracted {extracted} -> {} triangles, max error {:.2e}, total {total:.3}s",
        mesh.faces.len(),
        stats.maximum_error
    );
    println!("  {:15} {extraction:.3}s", "extraction");
    let names = [
        "sizing",
        "split",
        "collapse",
        "flip",
        "relax",
        "enforce+validate",
    ];
    for (name, s) in names.iter().zip(stats.seconds) {
        println!("  {name:15} {s:.3}s");
    }
}
