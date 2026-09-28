//! Small numeric WASM ABI. Output JSON remains owned by Rust until the next call.
use meshers_core::{
    Point, cross, dot,
    implicit::{self, ScalarField},
    norm, sub,
    surface_band::Band,
    triangles,
};
use serde_json::{Value, json};
use std::{
    cell::RefCell,
    collections::{BTreeSet, HashMap},
    f64::consts::TAU,
};

thread_local! { static OUTPUT: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) }; }

struct Tpms {
    shape: u32,
    thickness: f64,
    grade: f64,
    repeat: f64,
}
impl Tpms {
    fn raw(&self, p: Point) -> f64 {
        let [x, y, z] = p.map(|v| v * TAU);
        match self.shape {
            0 => x.sin() * y.cos() + y.sin() * z.cos() + z.sin() * x.cos(),
            1 => {
                1.1 * ((2. * x).sin() * y.cos() * z.sin()
                    + (2. * y).sin() * z.cos() * x.sin()
                    + (2. * z).sin() * x.cos() * y.sin())
                    - 0.2
                        * ((2. * x).cos() * (2. * y).cos()
                            + (2. * y).cos() * (2. * z).cos()
                            + (2. * z).cos() * (2. * x).cos())
                    - 0.4 * ((2. * x).cos() + (2. * y).cos() + (2. * z).cos())
            }
            _ => x.cos() + y.cos() + z.cos(),
        }
    }
}
impl ScalarField for Tpms {
    fn value(&self, p: Point) -> f64 {
        self.raw(p) / (0.5 * self.thickness * (1. + 2. * self.grade * p[0] / self.repeat))
    }
    fn gradient(&self, p: Point) -> Option<Point> {
        let q = p.map(|v| v * TAU);
        let s = q.map(f64::sin);
        let c = q.map(f64::cos);
        let g = match self.shape {
            0 => std::array::from_fn(|a| c[a] * c[(a + 1) % 3] - s[(a + 2) % 3] * s[a]),
            1 => {
                let s2 = q.map(|v| (2. * v).sin());
                let c2 = q.map(|v| (2. * v).cos());
                std::array::from_fn(|a| {
                    let b = (a + 1) % 3;
                    let d = (a + 2) % 3;
                    1.1 * (2. * c2[a] * c[b] * s[d] + s2[b] * c[d] * c[a] - s2[d] * s[a] * s[b])
                        + 0.4 * s2[a] * (c2[b] + c2[d])
                        + 0.8 * s2[a]
                })
            }
            _ => s.map(|v| -v),
        };
        let width = 0.5 * self.thickness * (1. + 2. * self.grade * p[0] / self.repeat);
        let mut result = g.map(|v| v * TAU / width);
        // Quotient rule for the X-dependent sheet width.
        if self.grade != 0. {
            result[0] -= self.raw(p) * self.thickness * self.grade / self.repeat / width.powi(2);
        }
        Some(result)
    }
}

fn percentile(v: &mut [f64], fraction: f64) -> f64 {
    if v.is_empty() {
        return 0.;
    }
    v.sort_by(f64::total_cmp);
    v[((v.len() - 1) as f64 * fraction).round() as usize]
}

fn generate(
    shape: u32,
    preset: u32,
    repeat: u32,
    resolution: u32,
    thickness: f64,
    grade: f64,
) -> Result<Value, String> {
    if shape > 2
        || preset > 2
        || !(1..=3).contains(&repeat)
        || !(8..=32).contains(&resolution)
        || !thickness.is_finite()
        || !(0.2..=1.2).contains(&thickness)
        || !grade.is_finite()
        || grade.abs() > 0.7
    {
        return Err(
            "Choose 1–3 cells, 8–32 samples per cell, thickness 0.2–1.2, and grading up to 70%."
                .into(),
        );
    }
    let cells = (repeat * resolution - 1) as usize;
    if cells > 64 || (preset == 2 && cells > 40) {
        return Err("Browser budget exceeded. Reduce sampling or the number of cells. Surfaces allow 64 grid divisions; volumes allow 40.".into());
    }
    let r = repeat as f64;
    let field = Tpms {
        shape,
        thickness,
        grade,
        repeat: r,
    };
    let bounds = [[-r / 2.; 3], [r / 2.; 3]];
    let periodic = [grade == 0., true, true];
    let (points, faces, tets, labels, volume_metrics) = if preset == 2 {
        let output = implicit::generate(
            &field,
            implicit::Options {
                bounds,
                cells: [cells; 3],
                region: implicit::Region::Band {
                    lower: -1.,
                    upper: 1.,
                },
                periodic,
                geometry_tolerance: 0.01,
                minimum_quality: 0.,
                max_tetrahedra: 600_000,
                optimize_passes: 4,
                threads: 0,
                ..Default::default()
            },
        )
        .map_err(|e| e.to_string())?;
        let d = output.diagnostics;
        let metrics = json!({"minimum_quality":d.minimum_mmg_quality,"below_01":d.elements_below_quality_01,
            "sampled_error":d.maximum_sampled_surface_error,"volume":d.volume,
            "target_met":d.minimum_mmg_quality>=0.1 && d.maximum_sampled_surface_error<=0.01});
        // Convert volume boundary tags to surface labels: implicit=0, box=2..7.
        let labels = output
            .boundary_tags
            .iter()
            .map(|&v| if v == 0 { 0 } else { v + 1 })
            .collect::<Vec<_>>();
        (
            output.mesh.points,
            output.mesh.surface,
            output.mesh.tets,
            labels,
            metrics,
        )
    } else {
        let band = Band {
            field: &field,
            bounds,
            levels: [-1., 1.],
        };
        let mesh = triangles::extract_with_edge_refinement(&band, cells, preset == 1)
            .map_err(|e| e.to_string())?;
        (
            mesh.points,
            mesh.faces,
            Vec::new(),
            mesh.labels,
            Value::Null,
        )
    };
    let mut angles = Vec::with_capacity(faces.len());
    let mut areas = Vec::with_capacity(faces.len());
    let mut distances = Vec::new();
    let mut edge_counts = HashMap::<(usize, usize), (usize, i32)>::new();
    for (face, &label) in faces.iter().zip(&labels) {
        let [a, b, c] = face.map(|i| points[i]);
        let area = norm(cross(sub(b, a), sub(c, a))) * 0.5;
        if !area.is_finite() || area <= 0. {
            return Err("The mesh contains a degenerate face. Try finer sampling.".into());
        }
        areas.push(area);
        let p = [a, b, c];
        let angle = (0..3)
            .map(|i| {
                let u = sub(p[(i + 1) % 3], p[i]);
                let v = sub(p[(i + 2) % 3], p[i]);
                (dot(u, v) / (norm(u) * norm(v)))
                    .clamp(-1., 1.)
                    .acos()
                    .to_degrees()
            })
            .fold(180., f64::min);
        angles.push(angle);
        if label < 2 {
            let center = std::array::from_fn(|i| (a[i] + b[i] + c[i]) / 3.);
            distances.push(
                (field.value(center).abs() - 1.).abs()
                    / norm(field.gradient(center).unwrap()).max(1e-12),
            );
        }
        for i in 0..3 {
            let a = face[i];
            let b = face[(i + 1) % 3];
            let entry = edge_counts.entry((a.min(b), a.max(b))).or_default();
            entry.0 += 1;
            entry.1 += if a < b { 1 } else { -1 };
        }
    }
    let bad_edges = edge_counts
        .values()
        .filter(|&&(n, w)| n != 2 || w != 0)
        .count();
    let cap_key = |face: &[usize; 3], axis: usize| {
        let mut key = face.map(|i| {
            let mut p = points[i];
            p[axis] = 0.;
            p.map(|x| (x * 1e8).round() as i64)
        });
        key.sort_unstable();
        key
    };
    let periodic_matches: [Option<bool>; 3] = std::array::from_fn(|axis| {
        if !periodic[axis] {
            return None;
        }
        let sides = [0, 1].map(|side| {
            faces
                .iter()
                .zip(&labels)
                .filter(|&(_, l)| *l == (2 + 2 * axis + side) as u8)
                .map(|(f, _)| cap_key(f, axis))
                .collect::<BTreeSet<_>>()
        });
        Some(!sides[0].is_empty() && sides[0] == sides[1])
    });
    let mean = areas.iter().sum::<f64>() / areas.len() as f64;
    let cv =
        (areas.iter().map(|a| (a - mean).powi(2)).sum::<f64>() / areas.len() as f64).sqrt() / mean;
    let mut sorted_angles = angles.clone();
    let p01 = percentile(&mut sorted_angles, 0.01);
    let hist: Vec<usize> = (0..12)
        .map(|bin| {
            angles
                .iter()
                .filter(|&&a| (a / 5.).floor().min(11.) as usize == bin)
                .count()
        })
        .collect();
    Ok(
        json!({"points":points,"triangles":faces,"tetrahedra":tets,"angles":angles,"labels":labels,
        "metrics":{"angle_p01":p01,"angle_min":sorted_angles[0],"area_cv":cv,
            "area_p01":percentile(&mut areas,0.01),"area_median":percentile(&mut areas,0.5),"area_p99":percentile(&mut areas,0.99),
            "centroid_distance_p95":percentile(&mut distances,0.95),"bad_edges":bad_edges,"periodic_matches":periodic_matches,
            "histogram":hist,"cells":cells,"volume":volume_metrics}}),
    )
}

#[unsafe(no_mangle)]
pub extern "C" fn meshers_generate(
    shape: u32,
    preset: u32,
    repeat: u32,
    resolution: u32,
    thickness: f64,
    grade: f64,
) -> u32 {
    let result = match generate(shape, preset, repeat, resolution, thickness, grade) {
        Ok(v) => json!({"ok":true,"mesh":v}),
        Err(e) => json!({"ok":false,"error":e}),
    };
    let output = serde_json::to_vec(&result).unwrap();
    let length = output.len() as u32;
    OUTPUT.with(|v| *v.borrow_mut() = output);
    length
}
#[unsafe(no_mangle)]
pub extern "C" fn meshers_output_ptr() -> usize {
    OUTPUT.with(|v| v.borrow().as_ptr() as usize)
}
#[unsafe(no_mangle)]
pub extern "C" fn meshers_clear() {
    OUTPUT.with(|v| *v.borrow_mut() = Vec::new());
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    #[ignore = "manual opt-in optimizer experiment matrix"]
    fn profile_optimizer_matrix() {
        use std::hash::{Hash, Hasher};
        for (name, shape, repeat, resolution, grade) in [
            ("gyroid", 0, 1, 16, 0.),
            ("split-p", 1, 1, 16, 0.),
            ("schwarz-p", 2, 1, 16, 0.),
            ("graded-gyroid-2", 0, 2, 12, 0.3),
            ("graded-split-p-2", 1, 2, 12, 0.3),
        ] {
            let mut times = Vec::new();
            let mut output = Value::Null;
            for _ in 0..3 {
                let start = std::time::Instant::now();
                output = match generate(shape, 2, repeat, resolution, 0.6, grade) {
                    Ok(mesh) => mesh,
                    Err(error) => {
                        println!(
                            "EXPERIMENT {}",
                            json!({"case":name,"error":error,"seconds":start.elapsed().as_secs_f64()})
                        );
                        break;
                    }
                };
                times.push(start.elapsed().as_secs_f64());
            }
            if output.is_null() {
                continue;
            }
            times.sort_by(f64::total_cmp);
            let mut hash = std::collections::hash_map::DefaultHasher::new();
            serde_json::to_vec(&output).unwrap().hash(&mut hash);
            println!(
                "EXPERIMENT {}",
                json!({"case":name,"seconds":times[1],"times":times,
                "fingerprint":format!("{:016x}",hash.finish()),"tets":output["tetrahedra"].as_array().unwrap().len(),"metrics":output["metrics"]})
            );
        }
    }
    #[test]
    #[ignore = "manual native adapter benchmark"]
    fn profile_fea_generation() {
        for shape in [0, 1] {
            let mut times = Vec::new();
            for _ in 0..4 {
                let start = std::time::Instant::now();
                let mesh = generate(shape, 2, 1, 16, 0.6, 0.).unwrap();
                let _serialized = serde_json::to_vec(&mesh).unwrap();
                times.push(start.elapsed().as_secs_f64());
                assert_eq!(mesh["metrics"]["volume"]["target_met"], true);
            }
            times.remove(0);
            times.sort_by(f64::total_cmp);
            println!("FEA shape={shape} median={:.6}", times[1]);
        }
    }
    #[test]
    #[ignore = "manual native generation benchmark"]
    fn profile_surface_extraction() {
        for shape in [0, 1] {
            for repeat in [1., 3., 7.] {
                let field = Tpms {
                    shape,
                    thickness: 0.6,
                    grade: 0.3,
                    repeat,
                };
                let band = Band {
                    field: &field,
                    bounds: [[-repeat / 2.; 3], [repeat / 2.; 3]],
                    levels: [-1., 1.],
                };
                let mut times = Vec::new();
                for run in 0..4 {
                    let start = std::time::Instant::now();
                    let mesh = triangles::extract_with_edge_refinement(
                        &band,
                        (repeat * 16.) as usize - 1,
                        false,
                    )
                    .unwrap();
                    let seconds = start.elapsed().as_secs_f64();
                    if run > 0 {
                        times.push(seconds);
                    }
                    assert!(!mesh.faces.is_empty());
                }
                times.sort_by(f64::total_cmp);
                println!(
                    "surface shape={shape} repeats={repeat} median={:.6}",
                    times[1]
                );
            }
        }
    }
    #[test]
    fn analytic_gradients_match_independent_central_differences() {
        for shape in 0..3 {
            for grade in [-0.7, 0., 0.7] {
                for repeat in [1., 3.] {
                    let field = Tpms {
                        shape,
                        thickness: 0.6,
                        grade,
                        repeat,
                    };
                    for i in 0..41 {
                        let p = std::array::from_fn(|a| {
                            (((i * (a * 12 + 7) + a * 11) % 101) as f64 / 100. - 0.5) * repeat
                        });
                        let g = field.gradient(p).unwrap();
                        for a in 0..3 {
                            let mut lo = p;
                            let mut hi = p;
                            lo[a] -= 1e-6;
                            hi[a] += 1e-6;
                            let expected = (field.value(hi) - field.value(lo)) / 2e-6;
                            assert!(
                                (g[a] - expected).abs() < 1e-6 * (1. + expected.abs()),
                                "shape={shape} grade={grade} point={p:?} axis={a}"
                            );
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn invalid_inputs_are_errors() {
        assert!(generate(0, 0, 1, 12, f64::NAN, 0.).is_err());
        assert!(generate(0, 2, 3, 32, 0.6, 0.).is_err());
    }
    #[test]
    fn surface_presets_are_closed_and_periodic() {
        for preset in [0, 1] {
            let mesh = generate(0, preset, 1, 12, 0.6, 0.).unwrap();
            assert_eq!(mesh["metrics"]["bad_edges"], 0);
            assert_eq!(
                mesh["metrics"]["periodic_matches"],
                json!([true, true, true])
            );
        }
    }
}
