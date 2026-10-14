//! Periodic Gyroid sheet adaptation with per-phase timing.
//! cargo run --release --example adapt_gyroid -- 1e-3
use meshers_core::implicit::{Options, Region, generate};
use meshers_core::{Point, tet_adapt};

fn main() {
    let tolerance: f64 = std::env::args().nth(1).map_or(1e-3, |v| v.parse().unwrap());
    let field = |p: Point| {
        let t = std::f64::consts::TAU;
        ((t * p[0]).sin() * (t * p[1]).cos()
            + (t * p[1]).sin() * (t * p[2]).cos()
            + (t * p[2]).sin() * (t * p[0]).cos())
            / 0.55
    };
    let base = Options {
        bounds: [[-0.5; 3], [0.5; 3]],
        cells: [24; 3],
        region: Region::Band {
            lower: -1.,
            upper: 1.,
        },
        periodic: [true; 3],
        geometry_tolerance: tolerance,
        optimize_passes: 1,
        ..Default::default()
    };
    let start = std::time::Instant::now();
    let out = generate(
        &field,
        Options {
            adapt: Some(tet_adapt::Options {
                tolerance,
                min_length: tolerance,
                max_length: 0.1,
                ..Default::default()
            }),
            ..base
        },
    )
    .unwrap();
    let total = start.elapsed().as_secs_f64();
    let stats = out.diagnostics.adapt.unwrap();
    println!(
        "tets {} total {total:.2}s q_min {:.3}",
        out.mesh.tets.len(),
        out.diagnostics.minimum_mmg_quality
    );
    let names = [
        "sizing",
        "split",
        "collapse",
        "edge removal",
        "face swap",
        "relocate",
        "enforce+unfold",
    ];
    for (name, s) in names.iter().zip(stats.seconds) {
        println!("  {name:15} {s:.3}s");
    }
}
