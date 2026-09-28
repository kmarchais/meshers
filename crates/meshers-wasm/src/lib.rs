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
        Some(std::array::from_fn(|a| {
            let mut lo = p;
            let mut hi = p;
            lo[a] -= 1e-5;
            hi[a] += 1e-5;
            (self.value(hi) - self.value(lo)) / 2e-5
        }))
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
