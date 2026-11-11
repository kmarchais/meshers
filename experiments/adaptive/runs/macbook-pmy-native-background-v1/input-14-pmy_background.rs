//! Bounded PMY experiment runner for a frozen, independently preflighted background.

#[cfg(feature = "experimental-quads")]
mod experiment {
    use meshers_core::{
        Point, Tet,
        implicit::{self, Options, Region, ScalarField},
    };
    use serde::Deserialize;
    use serde_json::json;
    use std::{
        fs::File,
        io::{BufReader, BufWriter},
        time::Instant,
    };

    #[derive(Deserialize)]
    struct Input {
        points: Vec<Point>,
        tetrahedra: Vec<Tet>,
        levels: [f64; 2],
        probes: Vec<Point>,
        builtin_background: bool,
    }

    struct Pmy;
    impl ScalarField for Pmy {
        fn value(&self, p: Point) -> f64 {
            let [x, y, z] = p.map(|v| v * std::f64::consts::TAU);
            2. * x.cos() * y.cos() * z.cos()
                + (2. * x).sin() * y.sin()
                + x.sin() * (2. * z).sin()
                + (2. * y).sin() * z.sin()
        }
        fn gradient(&self, p: Point) -> Option<Point> {
            let [x, y, z] = p.map(|v| v * std::f64::consts::TAU);
            Some(
                [
                    -2. * x.sin() * y.cos() * z.cos()
                        + 2. * (2. * x).cos() * y.sin()
                        + x.cos() * (2. * z).sin(),
                    -2. * x.cos() * y.sin() * z.cos()
                        + (2. * x).sin() * y.cos()
                        + 2. * (2. * y).cos() * z.sin(),
                    -2. * x.cos() * y.cos() * z.sin()
                        + 2. * x.sin() * (2. * z).cos()
                        + (2. * y).sin() * z.cos(),
                ]
                .map(|v| v * std::f64::consts::TAU),
            )
        }
    }

    pub fn run() -> Result<(), Box<dyn std::error::Error>> {
        let args: Vec<_> = std::env::args().collect();
        if args.len() != 3 {
            return Err("expected input JSON and output JSON paths".into());
        }
        let input: Input = serde_json::from_reader(BufReader::new(File::open(&args[1])?))?;
        let options = Options {
            bounds: [[-0.5; 3], [0.5; 3]],
            cells: [16; 3],
            region: Region::Band {
                lower: input.levels[0],
                upper: input.levels[1],
            },
            periodic: [true; 3],
            geometry_tolerance: 1.,
            minimum_quality: 0.1,
            max_tetrahedra: 250_000,
            optimize_passes: 4,
            snap: 0.2,
            threads: 1,
            order: 2,
            curved_fit: true,
            ..Default::default()
        };
        let start = Instant::now();
        let output = if input.builtin_background {
            implicit::generate(&Pmy, options)
        } else {
            implicit::generate_with_background(&Pmy, options, &input.points, &input.tetrahedra)
        };
        let seconds = start.elapsed().as_secs_f64();
        let mut record = json!({
            "native_generation_s": seconds,
            "builtin_background": input.builtin_background,
            "probe_values": Pmy.values(&input.probes),
            "probe_gradients": Pmy.gradients(&input.probes),
            "options": {"cells": 16, "order": 2, "adapt": false,
                "snap": 0.2, "optimize_passes": 4, "curved_fit": true,
                "geometry_tolerance": 1., "minimum_quality": 0.1,
                "max_tetrahedra": 250_000, "threads": 1,
                "periodic": [true, true, true], "levels": input.levels},
            "admitted": false,
            "scope": "Native experiment; independent unchanged admission gates are applied afterward."
        });
        match output {
            Ok(output) => {
                let quadratic = output.quadratic.ok_or("expected quadratic output")?;
                record["status"] = json!("returned");
                record["points"] = json!(quadratic.elements.points);
                record["tetrahedra"] = json!(quadratic.elements.tets);
                record["surface"] = json!(quadratic.elements.surface);
                record["boundary_tags"] = json!(output.boundary_tags);
                record["periodic_pairs"] = json!(quadratic.periodic_pairs);
                record["diagnostics"] = json!({
                    "volume_linear": output.diagnostics.volume,
                    "minimum_mmg_quality": output.diagnostics.minimum_mmg_quality,
                    "elements_below_quality_01": output.diagnostics.elements_below_quality_01,
                    "reverted_snap_vertices": output.diagnostics.reverted_snap_vertices,
                    "maximum_sampled_surface_error": output.diagnostics.maximum_sampled_surface_error,
                    "minimum_jacobian_ratio": quadratic.elements.minimum_jacobian_ratio,
                    "curved_edges": quadratic.elements.curved_edges,
                    "straightened_edges": quadratic.elements.straightened_edges,
                });
            }
            Err(error) => {
                record["status"] = json!("native-rejected");
                record["error"] = json!(error.to_string());
            }
        }
        serde_json::to_writer(BufWriter::new(File::create(&args[2])?), &record)?;
        println!("{}: {:.3} s", record["status"], seconds);
        Ok(())
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    #[cfg(feature = "experimental-quads")]
    return experiment::run();
    #[cfg(not(feature = "experimental-quads"))]
    Err("pmy_background requires experimental-quads".into())
}
