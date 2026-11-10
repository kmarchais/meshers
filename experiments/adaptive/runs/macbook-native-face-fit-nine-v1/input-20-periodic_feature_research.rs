//! Explicit saved-seed research controls; never part of the production mesher.
use super::{Context, Options, ScalarField, finish};
use crate::{Mesh, Point};
use std::sync::atomic::AtomicBool;

struct Wedge {
    notched: bool,
}
impl ScalarField for Wedge {
    fn value(&self, p: Point) -> f64 {
        let base = (p[0].abs() + 0.3 * p[2].abs() - 0.35).max(p[1].abs() - 0.2);
        if self.notched {
            base.max((p[0] - 0.05).min(p[1] - 0.05))
        } else {
            base
        }
    }
    fn gradient(&self, p: Point) -> Option<Point> {
        let sign = |v: f64| {
            if v > 0. {
                1.
            } else if v < 0. {
                -1.
            } else {
                0.
            }
        };
        let a = p[0].abs() + 0.3 * p[2].abs() - 0.35;
        let b = p[1].abs() - 0.2;
        let x = [sign(p[0]), 0., 0.3 * sign(p[2])];
        let y = [0., sign(p[1]), 0.];
        let base = if a > b {
            x
        } else if b > a {
            y
        } else {
            std::array::from_fn(|i| 0.5 * x[i] + 0.5 * y[i])
        };
        if !self.notched {
            return Some(base);
        }
        let cut = if p[0] < p[1] {
            [1., 0., 0.]
        } else if p[1] < p[0] {
            [0., 1., 0.]
        } else {
            [0.5, 0.5, 0.]
        };
        Some(if a.max(b) > (p[0] - 0.05).min(p[1] - 0.05) {
            base
        } else if a.max(b) < (p[0] - 0.05).min(p[1] - 0.05) {
            cut
        } else {
            std::array::from_fn(|i| 0.5 * base[i] + 0.5 * cut[i])
        })
    }
}

#[test]
#[ignore = "bounded fixed-seed research; requires explicit MESHERS_WEDGE_SEED_JSON and MESHERS_WEDGE_OUTPUT_JSON"]
fn adapt_saved_wedge_seed() {
    let input = std::env::var("MESHERS_WEDGE_SEED_JSON").expect("explicit input path");
    let output = std::env::var("MESHERS_WEDGE_OUTPUT_JSON").expect("explicit output path");
    let data: serde_json::Value = serde_json::from_slice(&std::fs::read(input).unwrap()).unwrap();
    let bits: Vec<[u64; 3]> = serde_json::from_value(data["point_bits"].clone()).unwrap();
    let plane_bits: Vec<[u64; 4]> = serde_json::from_value(data["plane_bits"].clone()).unwrap();
    let raw_planes: Vec<_> = plane_bits
        .into_iter()
        .map(|p| p.map(f64::from_bits))
        .collect();
    let planes = crate::polyhedral::normalized(&raw_planes).unwrap();
    let regions: Option<Vec<super::intersection::ConstraintMask>> = data
        .get("region_masks")
        .map(|value| serde_json::from_value(value.clone()).unwrap());
    let field = Wedge {
        notched: regions.is_some(),
    };
    let mut mesh = Mesh {
        points: bits.into_iter().map(|p| p.map(f64::from_bits)).collect(),
        tets: serde_json::from_value(data["tetrahedra"].clone()).unwrap(),
        surface: serde_json::from_value(data["surface"].clone()).unwrap(),
    };
    let adapt = crate::tet_adapt::Options {
        tolerance: 0.002,
        min_length: 0.002,
        max_length: 0.1,
        gradation: 0.3,
        layers: 0.,
        iterations: 4,
        quality_guard: 0.3,
    };
    let options = Options {
        bounds: [[-0.5; 3], [0.5; 3]],
        cells: [24; 3],
        periodic: [false, false, true],
        minimum_quality: 0.1,
        geometry_tolerance: 0.002,
        max_tetrahedra: 250000,
        optimize_passes: 1,
        snap: 0.2,
        threads: 1,
        adapt: Some(adapt),
        ..Options::default()
    };
    let mut scale = 0f64;
    for x in [0., 0.23, 0.61, 1.] {
        for y in [0., 0.23, 0.61, 1.] {
            for z in [0., 0.23, 0.61, 1.] {
                scale = scale.max(field.value([x - 0.5, y - 0.5, z - 0.5]).abs());
            }
        }
    }
    let context = Context {
        field: &field,
        constraints: None,
        halfspaces: Some(&planes),
        halfspace_regions: regions.as_deref(),
        size_field: None,
        options,
        center: 0.,
        scale,
        optimizer_length: 1.,
        invalid: AtomicBool::new(false),
    };
    let stats = crate::tet_adapt::adapt_with_polish(
        &mut mesh,
        &context,
        0.,
        options.periodic,
        &adapt,
        options.max_tetrahedra,
        options.optimize_passes,
        None,
    )
    .unwrap();
    let point_bits: Vec<_> = mesh.points.iter().map(|p| p.map(f64::to_bits)).collect();
    let mut result = serde_json::json!({"point_bits":point_bits, "tetrahedra":mesh.tets,
        "surface":mesh.surface, "adapt":{"splits":stats.splits,"collapses":stats.collapses,
            "swaps":stats.swaps,"moves":stats.moves,"error_splits":stats.error_splits,
            "maximum_error":stats.maximum_error}, "scale":scale,
        "scope":"Fixed-seed private-core adaptation only; no public generator periodic gate was relaxed"});
    let probes: Vec<[u64; 3]> = serde_json::from_value(data["probe_bits"].clone()).unwrap();
    let probes: Vec<Point> = probes.into_iter().map(|p| p.map(f64::from_bits)).collect();
    result["probe_value_bits"] = serde_json::to_value(
        probes
            .iter()
            .map(|&p| field.value(p).to_bits())
            .collect::<Vec<_>>(),
    )
    .unwrap();
    result["probe_gradient_bits"] = serde_json::to_value(
        probes
            .iter()
            .map(|&p| field.gradient(p).unwrap().map(f64::to_bits))
            .collect::<Vec<_>>(),
    )
    .unwrap();
    match finish(mesh, &context) {
        Ok(out) => {
            result["native_finish_admitted"] = true.into();
            result["native_minimum_quality"] = out.diagnostics.minimum_mmg_quality.into();
            result["native_surface_error"] = out.diagnostics.maximum_sampled_surface_error.into();
            result["native_periodic_pair_counts"] = out.periodic_pairs.map(|p| p.len()).into();
        }
        Err(error) => {
            result["native_finish_admitted"] = false.into();
            result["native_finish_error"] = error.into();
        }
    }
    std::fs::write(output, serde_json::to_vec(&result).unwrap()).unwrap();
}
