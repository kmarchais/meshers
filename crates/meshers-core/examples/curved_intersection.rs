//! Experimental native intersection runner; dimensions match the Python specimen.
use meshers_core::{
    Point,
    implicit::{self, intersection},
};
use serde::Serialize;
use std::{fs::File, io::BufWriter, time::Instant};
#[derive(Serialize)]
struct Export<'a> {
    points: &'a [Point],
    parameters: &'a [Point],
    tetrahedra: &'a [[usize; 4]],
    surface: &'a [[usize; 3]],
    constraints: &'a [u8],
    seconds: f64,
    minimum_quality: f64,
    elements_below_quality_01: usize,
    quality: &'a intersection::QualityDiagnostics,
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    let mode = args.get(1).map(String::as_str).unwrap_or("cylinder");
    if !["cylinder", "torus"].contains(&mode) {
        return Err("expected cylinder or torus".into());
    }
    let path = args.get(2).ok_or("expected output JSON path")?;
    let snap = args
        .get(3)
        .map(|v| v.parse::<f64>())
        .transpose()?
        .unwrap_or(0.15);
    let passes = args
        .get(4)
        .map(|v| v.parse::<usize>())
        .transpose()?
        .unwrap_or(0);
    let minimum_quality = args
        .get(5)
        .map(|v| v.parse::<f64>())
        .transpose()?
        .unwrap_or(0.);
    let gyroid = |p: Point| {
        let [a, b, c] =
            [p[0] / 0.75, p[1] / 0.75, p[2]].map(|v| std::f64::consts::TAU * (v + 0.125));
        a.sin() * b.cos() + b.sin() * c.cos() + c.sin() * a.cos()
    };
    let upper = |p| gyroid(p) - 0.55;
    let lower = |p| -gyroid(p) - 0.55;
    let cylinder = |p: Point| p[0].hypot(p[1]) - 0.75;
    let radius = 9. / std::f64::consts::PI;
    let map = |p: Point| {
        if mode == "torus" {
            let a = p[2] / radius;
            [(radius + p[0]) * a.cos(), (radius + p[0]) * a.sin(), -p[1]]
        } else {
            p
        }
    };
    let start = Instant::now();
    let result = intersection::generate(
        &[&upper, &lower, &cylinder],
        map,
        implicit::Options {
            bounds: [[-0.8, -0.8, 0.], [0.8, 0.8, 3.]],
            cells: [48, 48, 80],
            periodic: [false, false, true],
            snap,
            optimize_passes: passes,
            minimum_quality,
            geometry_tolerance: 0.,
            max_tetrahedra: 4_000_000,
            ..Default::default()
        },
    )?;
    let mut minimum_quality: f64 = 1.;
    let mut below = 0;
    for t in &result.mesh.tets {
        let q = meshers_core::quality(t.map(|i| result.mesh.points[i])).powf(1.5);
        minimum_quality = minimum_quality.min(q);
        below += usize::from(q < 0.1);
    }
    let seconds = start.elapsed().as_secs_f64();
    serde_json::to_writer(
        BufWriter::new(File::create(path)?),
        &Export {
            points: &result.mesh.points,
            parameters: &result.parameters,
            tetrahedra: &result.mesh.tets,
            surface: &result.mesh.surface,
            constraints: &result.constraints,
            seconds,
            minimum_quality,
            elements_below_quality_01: below,
            quality: &result.quality,
        },
    )?;
    println!(
        "{mode}: {} tetrahedra, min quality {minimum_quality:e}, {below} below 0.1, {seconds:.2}s",
        result.mesh.tets.len()
    );
    Ok(())
}
