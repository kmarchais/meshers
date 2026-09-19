use meshers_core::{
    Point,
    implicit::{Options, Region, ScalarField, generate},
};
use std::f64::consts::{PI, TAU};

fn gyroid(p: Point) -> f64 {
    let [x, y, z] = p.map(|v| v * TAU);
    x.sin() * y.cos() + y.sin() * z.cos() + z.sin() * x.cos()
}
fn schwarz_p(p: Point) -> f64 {
    p.map(|v| (TAU * v).cos()).iter().sum()
}
fn schwarz_d(p: Point) -> f64 {
    let [x, y, z] = p.map(|v| TAU * v);
    x.sin() * y.sin() * z.sin()
        + x.sin() * y.cos() * z.cos()
        + x.cos() * y.sin() * z.cos()
        + x.cos() * y.cos() * z.sin()
}
fn sphere(p: Point) -> f64 {
    ((p[0] - 0.47).powi(2) + (p[1] - 0.51).powi(2) + (p[2] - 0.49).powi(2)).sqrt() - 0.29
}

#[test]
fn empty_and_full_physical_boxes() {
    let options = Options {
        bounds: [[2., -1., 3.], [4., 2., 7.]],
        cells: [6; 3],
        ..Options::default()
    };
    let empty = generate(&|_| 1., options).unwrap();
    assert!(empty.diagnostics.empty);
    let full = generate(&|_| -1., options).unwrap();
    assert!((full.diagnostics.volume - 24.).abs() < 1e-10);
    assert!(full.boundary_tags.iter().all(|&t| t != 0));
    assert!(full.periodic_pairs.iter().all(Vec::is_empty));
}

#[test]
fn sphere_converges_and_physical_scaling_is_consistent() {
    let mut errors = Vec::new();
    for n in [10, 16, 24] {
        let o = Options {
            cells: [n; 3],
            geometry_tolerance: 0.03,
            optimize_passes: 1,
            ..Options::default()
        };
        let mesh = generate(&sphere, o).unwrap();
        errors.push((mesh.diagnostics.volume - 4. * PI * 0.29_f64.powi(3) / 3.).abs());
    }
    assert!(errors[2] < errors[1] && errors[1] < errors[0], "{errors:?}");
    let o = Options {
        cells: [12; 3],
        geometry_tolerance: 0.04,
        optimize_passes: 0,
        ..Options::default()
    };
    let a = generate(&sphere, o).unwrap();
    let b = generate(
        &|p: Point| sphere([(p[0] - 2.) / 3., (p[1] + 1.) / 3., (p[2] - 4.) / 3.]),
        Options {
            bounds: [[2., -1., 4.], [5., 2., 7.]],
            geometry_tolerance: 0.12,
            ..o
        },
    )
    .unwrap();
    assert!((b.diagnostics.volume / a.diagnostics.volume - 27.).abs() < 1e-7);
}

#[test]
fn periodic_tpms_band_and_skeletal_parts() {
    for field in [gyroid, schwarz_p, schwarz_d] {
        let o = Options {
            cells: [24; 3],
            periodic: [true; 3],
            region: Region::Band {
                lower: -0.5,
                upper: 0.5,
            },
            geometry_tolerance: 0.03,
            optimize_passes: 1,
            ..Options::default()
        };
        let result = generate(&field, o).unwrap();
        assert!(result.periodic_pairs.iter().all(|p| !p.is_empty()));
        assert!(result.diagnostics.minimum_mmg_quality > 0.);
        for sign in [-1., 1.] {
            let skeletal = generate(
                &|p| sign * field(p) + 0.25,
                Options {
                    region: Region::Negative,
                    ..o
                },
            )
            .unwrap();
            assert!(!skeletal.diagnostics.empty);
        }
    }
}

#[test]
fn custom_blend_and_variable_width() {
    let field = |p: Point| {
        let w = (1. + (TAU * p[0]).sin()) / 2.;
        let raw = (1. - w) * gyroid(p) + w * schwarz_p(p);
        raw.abs() - (0.6 + 0.1 * (TAU * p[2]).cos())
    };
    let result = generate(
        &field,
        Options {
            cells: [28; 3],
            geometry_tolerance: 0.04,
            optimize_passes: 1,
            ..Options::default()
        },
    )
    .unwrap();
    assert!(!result.diagnostics.empty);
    assert!(result.periodic_pairs.iter().all(Vec::is_empty));
}

#[test]
fn incompatible_grading_is_rejected_only_when_periodicity_is_requested() {
    let field = |p: Point| gyroid(p) - 0.1 - 0.2 * p[0];
    let o = Options {
        cells: [24; 3],
        geometry_tolerance: 0.04,
        optimize_passes: 1,
        ..Options::default()
    };
    assert!(generate(&field, o).is_ok());
    assert!(
        generate(
            &field,
            Options {
                periodic: [true, false, false],
                ..o
            }
        )
        .is_err()
    );
}

#[test]
fn explicit_axis_pairing_and_nonperiodic_coordinates() {
    // A slanted cut deliberately has different traces on the x faces.
    let field = |p: Point| p[1] - 0.2 - 0.4 * p[0];
    let o = Options {
        cells: [8; 3],
        periodic: [false, false, true],
        geometry_tolerance: 1e-6,
        ..Options::default()
    };
    let mesh = generate(&field, o).unwrap();
    assert!(mesh.periodic_pairs[0].is_empty() && mesh.periodic_pairs[1].is_empty());
    assert!(!mesh.periodic_pairs[2].is_empty());
    assert!((mesh.diagnostics.volume - 0.4).abs() < 1e-8);
}

struct AnalyticSphere;
impl ScalarField for AnalyticSphere {
    fn value(&self, p: Point) -> f64 {
        sphere(p)
    }
    fn gradient(&self, p: Point) -> Option<Point> {
        let d = [p[0] - 0.47, p[1] - 0.51, p[2] - 0.49];
        let length = d.iter().map(|v| v * v).sum::<f64>().sqrt();
        Some(d.map(|v| v / length))
    }
}
#[test]
fn supplied_and_finite_difference_gradients_agree() {
    let o = Options {
        cells: [12; 3],
        geometry_tolerance: 0.03,
        ..Options::default()
    };
    let a = generate(&sphere, o).unwrap();
    let b = generate(&AnalyticSphere, o).unwrap();
    assert!((a.diagnostics.volume - b.diagnostics.volume).abs() < 1e-6);
}

#[test]
fn invalid_fields_and_unmet_tolerance_never_succeed() {
    let o = Options {
        cells: [8; 3],
        ..Options::default()
    };
    assert!(generate(&|_| f64::NAN, o).is_err());
    assert!(
        generate(
            &sphere,
            Options {
                geometry_tolerance: 1e-10,
                ..o
            }
        )
        .is_err()
    );
    assert!(
        generate(
            &sphere,
            Options {
                max_tetrahedra: 5,
                ..o
            }
        )
        .is_err()
    );
    assert!(
        generate(
            &sphere,
            Options {
                bounds: [[0.; 3], [0.; 3]],
                ..o
            }
        )
        .is_err()
    );
}

#[test]
fn positive_field_scaling_preserves_geometry() {
    let o = Options {
        cells: [12; 3],
        geometry_tolerance: 0.03,
        optimize_passes: 0,
        ..Options::default()
    };
    let a = generate(&sphere, o).unwrap();
    let b = generate(&|p| 1e9 * sphere(p), o).unwrap();
    assert_eq!(a.mesh.tets, b.mesh.tets);
    assert!((a.diagnostics.volume - b.diagnostics.volume).abs() < 1e-10);
}

#[test]
fn constant_field_inside_band_fills_the_domain() {
    let mesh = generate(
        &|_| 0.,
        Options {
            cells: [6; 3],
            region: Region::Band {
                lower: -1.,
                upper: 1.,
            },
            ..Options::default()
        },
    )
    .unwrap();
    assert!((mesh.diagnostics.volume - 1.).abs() < 1e-10);
}

#[test]
fn two_small_interior_components_survive_sampling() {
    let field = |p: Point| {
        let sphere = |x: f64| {
            ((p[0] - x).powi(2) + (p[1] - 0.5).powi(2) + (p[2] - 0.5).powi(2)).sqrt() - 0.12
        };
        sphere(0.25).min(sphere(0.75))
    };
    let result = generate(
        &field,
        Options {
            cells: [20; 3],
            geometry_tolerance: 0.02,
            optimize_passes: 1,
            ..Options::default()
        },
    )
    .unwrap();
    let mut adjacency = vec![Vec::new(); result.mesh.points.len()];
    for t in &result.mesh.tets {
        for &i in t {
            for &j in t {
                adjacency[i].push(j);
            }
        }
    }
    let mut visited = vec![false; adjacency.len()];
    let mut components = 0;
    for i in 0..visited.len() {
        if visited[i] {
            continue;
        }
        components += 1;
        let mut pending = vec![i];
        visited[i] = true;
        while let Some(v) = pending.pop() {
            for &j in &adjacency[v] {
                if !visited[j] {
                    visited[j] = true;
                    pending.push(j);
                }
            }
        }
    }
    assert_eq!(components, 2);
}

#[test]
fn unequal_lengths_and_cell_counts_work_in_physical_space() {
    let field =
        |p: Point| ((p[0] / 2.).powi(2) + (p[1] / 3.).powi(2) + (p[2] / 4.).powi(2)).sqrt() - 0.3;
    let result = generate(
        &field,
        Options {
            bounds: [[-1., -1.5, -2.], [1., 1.5, 2.]],
            cells: [12, 16, 20],
            geometry_tolerance: 0.06,
            optimize_passes: 1,
            ..Options::default()
        },
    )
    .unwrap();
    let exact = 4. * PI * 0.6 * 0.9 * 1.2 / 3.;
    assert!((result.diagnostics.volume / exact - 1.).abs() < 0.04);
}

struct InvalidGradient;
impl ScalarField for InvalidGradient {
    fn value(&self, p: Point) -> f64 {
        sphere(p)
    }
    fn gradient(&self, _: Point) -> Option<Point> {
        Some([f64::NAN; 3])
    }
}
#[test]
fn nonfinite_gradient_never_returns_a_mesh() {
    assert!(
        generate(
            &InvalidGradient,
            Options {
                cells: [8; 3],
                ..Options::default()
            }
        )
        .is_err()
    );
}

#[test]
fn graded_blend_snapping_avoids_tiny_cut_cells() {
    let field =
        |p: Point| ((1. - p[0]) * gyroid(p) + p[0] * schwarz_p(p)).abs() - (0.5 + 0.2 * p[2]);
    let result = generate(
        &field,
        Options {
            cells: [28; 3],
            geometry_tolerance: 0.04,
            ..Options::default()
        },
    )
    .unwrap();
    assert!(result.diagnostics.minimum_mmg_quality > 0.2);
    assert!(result.diagnostics.maximum_sampled_surface_error < 0.01);
}

#[test]
fn asymmetric_band_uses_raw_endpoints() {
    let result = generate(
        &|p: Point| p[0],
        Options {
            cells: [8; 3],
            region: Region::Band {
                lower: 0.2,
                upper: 0.6,
            },
            periodic: [false, true, true],
            geometry_tolerance: 1e-6,
            ..Options::default()
        },
    )
    .unwrap();
    assert!((result.diagnostics.volume - 0.4).abs() < 1e-9);
}

#[test]
fn nonfinite_field_during_lattice_evaluation_is_rejected() {
    let field = |p: Point| if p[0] == 0.25 { f64::NAN } else { sphere(p) };
    assert!(
        generate(
            &field,
            Options {
                cells: [8; 3],
                ..Options::default()
            }
        )
        .is_err()
    );
}

#[test]
fn general_field_parallel_coloring_is_deterministic() {
    let o = Options {
        cells: [12; 3],
        threads: 1,
        geometry_tolerance: 0.03,
        ..Options::default()
    };
    let a = generate(&sphere, o).unwrap();
    let b = generate(&sphere, Options { threads: 4, ..o }).unwrap();
    assert_eq!(a.mesh.points, b.mesh.points);
    assert_eq!(a.mesh.tets, b.mesh.tets);
}

#[test]
fn twisted_sheet_keeps_a_manifold_boundary_when_snapped() {
    let field = |p: Point| {
        let s = ((p[0] + 1.5) / 3.).clamp(0., 1.);
        let angle = (-25. + 50. * s) * PI / 180.;
        let q = [
            p[0],
            angle.cos() * p[1] - angle.sin() * p[2],
            angle.sin() * p[1] + angle.cos() * p[2],
        ];
        (gyroid(q.map(|v| v + 0.125)) / 0.75_f64.sqrt()).abs() - 0.55
    };
    let result = generate(
        &field,
        Options {
            bounds: [[-1.5; 3], [1.5; 3]],
            cells: [64; 3],
            snap: 0.2,
            optimize_passes: 1,
            geometry_tolerance: 0.04,
            ..Options::default()
        },
    )
    .unwrap();
    assert!(result.diagnostics.reverted_snap_vertices > 0);
    assert!(result.diagnostics.minimum_mmg_quality > 1e-4);
    assert!(result.diagnostics.maximum_sampled_surface_error <= 0.04);
}

#[test]
fn minimum_quality_gate_rejects_flat_cells_without_changing_geometry() {
    let options = Options {
        bounds: [[0., 0., 0.], [1., 1., 0.001]],
        cells: [4; 3],
        optimize_passes: 0,
        ..Options::default()
    };
    let mesh = generate(&|_| -1., options).unwrap();
    assert!(mesh.diagnostics.minimum_mmg_quality < 0.01);
    assert!(mesh.diagnostics.elements_below_quality_01 > 0);
    assert!((mesh.diagnostics.volume - 0.001).abs() < 1e-12);
    let rejected = generate(
        &|_| -1.,
        Options {
            minimum_quality: 0.05,
            ..options
        },
    );
    assert!(
        rejected
            .err()
            .unwrap()
            .to_string()
            .contains("minimum MMG quality")
    );
    let accepted = generate(
        &|_| -1.,
        Options {
            bounds: [[0.; 3], [1.; 3]],
            minimum_quality: 0.05,
            ..options
        },
    )
    .unwrap();
    assert!(accepted.diagnostics.minimum_mmg_quality >= 0.05);
    assert!(
        generate(
            &|_| 1.,
            Options {
                minimum_quality: 1.,
                ..options
            }
        )
        .unwrap()
        .diagnostics
        .empty
    );
    for value in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            generate(
                &|_| -1.,
                Options {
                    minimum_quality: value,
                    ..options
                }
            ),
            Err(meshers_core::MeshingError::InvalidOptions(_))
        ));
    }
}
