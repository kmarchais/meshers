//! Inspect the exact f64 input parsed by the PMY example, without meshing.
use meshers_core::{Point, Tet};
use serde::Deserialize;
use serde_json::json;
use std::{
    fs::File,
    io::{BufReader, BufWriter},
};

// Keep the typed deserialization and reader path identical to pmy_background.
#[derive(Deserialize)]
struct Input {
    points: Vec<Point>,
    tetrahedra: Vec<Tet>,
    levels: [f64; 2],
    probes: Vec<Point>,
    builtin_background: bool,
    #[serde(default)]
    trace_snapping: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 3 {
        return Err("expected input JSON and output JSON paths".into());
    }
    let input: Input = serde_json::from_reader(BufReader::new(File::open(&args[1])?))?;
    let point_bits: Vec<_> = input.points.iter().map(|p| p.map(f64::to_bits)).collect();
    let probe_bits: Vec<_> = input.probes.iter().map(|p| p.map(f64::to_bits)).collect();
    let serialized = serde_json::to_vec(&input.points)?;
    let parsed_again: Vec<Point> = serde_json::from_slice(&serialized)?;
    let repeated_bits: Vec<_> = parsed_again.iter().map(|p| p.map(f64::to_bits)).collect();
    let output = json!({
        "scope": "Exact typed serde_json reader parse; no mesh generation or changes to the input",
        "points": input.points,
        "point_bits": point_bits,
        "tetrahedra": input.tetrahedra,
        "level_bits": input.levels.map(f64::to_bits),
        "probe_bits": probe_bits,
        "builtin_background": input.builtin_background,
        "trace_snapping": input.trace_snapping,
        "serialize_reparse_point_bits": repeated_bits,
    });
    serde_json::to_writer(BufWriter::new(File::create(&args[2])?), &output)?;
    Ok(())
}
