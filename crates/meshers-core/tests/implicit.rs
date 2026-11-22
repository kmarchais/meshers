use meshers_core::{
    Point,
    implicit::{Options, Region, ScalarField, generate},
};
use std::f64::consts::{PI, TAU};

fn gyroid(p: Point) -> f64 {
    let [x, y, z] = p.map(|v| v * TAU);
    x.sin() * y.cos() + y.sin() * z.cos() + z.sin() * x.cos()
}

#[test]
fn curved_fit_preserves_corners_connectivity_and_periodic_traces() {
    let options = Options {
        bounds: [[-0.5; 3], [0.5; 3]],
        cells: [9; 3],
        region: Region::Band {
            lower: -0.55,
            upper: 0.55,
        },
        periodic: [true; 3],
        geometry_tolerance: 1.,
        optimize_passes: 1,
        order: 2,
        ..Options::default()
    };
    let before = generate(&gyroid, options).unwrap();
    let after = generate(
        &gyroid,
        Options {
            curved_fit: true,
            ..options
        },
    )
    .unwrap();
    assert_eq!(before.mesh.points, after.mesh.points);
    let a = before.quadratic.unwrap();
    let b = after.quadratic.unwrap();
    assert_eq!(a.elements.tets, b.elements.tets);
    assert_eq!(a.elements.surface, b.elements.surface);
    assert_eq!(a.elements.points.len(), b.elements.points.len());
    assert_eq!(a.periodic_pairs, b.periodic_pairs);
    let stats = b.elements.fit.unwrap();
    assert!(stats.selected_groups > 0 && stats.moved_nodes > 0);
    assert!(stats.final_distance_estimate <= stats.initial_distance_estimate);
    assert!(b.elements.minimum_jacobian_ratio > 1e-3);
    for (axis, pairs) in b.periodic_pairs.iter().enumerate() {
        for &[lo, hi] in pairs {
            for component in 0..3 {
                let shift = if component == axis { 1. } else { 0. };
                let before_error =
                    (a.elements.points[hi][component] - a.elements.points[lo][component] - shift)
                        .abs();
                let after_error =
                    (b.elements.points[hi][component] - b.elements.points[lo][component] - shift)
                        .abs();
                // Baseline projection has a 5.72e-15 residual in this FD control.
                // The fit must not increase it; moved transverse copies agree exactly.
                assert!(
                    before_error < 1e-12 && after_error <= before_error + f64::EPSILON,
                    "axis={axis} component={component} pair={lo},{hi} before={before_error:e} after={after_error:e} old={:?},{:?} new={:?},{:?}",
                    a.elements.points[lo],
                    a.elements.points[hi],
                    b.elements.points[lo],
                    b.elements.points[hi]
                );
                if component != axis
                    && (b.elements.points[hi] != a.elements.points[hi]
                        || b.elements.points[lo] != a.elements.points[lo])
                {
                    assert_eq!(
                        b.elements.points[hi][component],
                        b.elements.points[lo][component]
                    );
                }
            }
        }
    }
    assert!(
        generate(
            &gyroid,
            Options {
                order: 1,
                curved_fit: true,
                ..options
            }
        )
        .is_err()
    );
}

#[test]
fn interior_refinement_is_nested_and_preserves_periodic_boundaries() {
    let options = Options {
        cells: [4; 3],
        optimize_passes: 0,
        minimum_quality: 0.1,
        periodic: [true; 3],
        ..Options::default()
    };
    let seed = generate(&|_| -1., options).unwrap();
    let explicit_zero = generate(
        &|_| -1.,
        Options {
            interior_refine_nodes: 0,
            ..options
        },
    )
    .unwrap();
    assert_eq!(seed.mesh.points, explicit_zero.mesh.points);
    assert_eq!(seed.mesh.tets, explicit_zero.mesh.tets);
    let refined = generate(
        &|_| -1.,
        Options {
            interior_refine_nodes: 20,
            ..options
        },
    )
    .unwrap();
    let inserted = refined.diagnostics.inserted_interior_nodes;
    assert!(inserted > 0 && inserted <= 20);
    assert_eq!(
        &refined.mesh.points[..seed.mesh.points.len()],
        &seed.mesh.points
    );
    assert_eq!(refined.mesh.surface, seed.mesh.surface);
    assert_eq!(refined.boundary_tags, seed.boundary_tags);
    assert_eq!(refined.periodic_pairs, seed.periodic_pairs);
    assert!((refined.diagnostics.volume - seed.diagnostics.volume).abs() < 1e-12);
    assert!(refined.diagnostics.minimum_mmg_quality >= 0.1);
    assert!(refined.mesh.tets.len() > seed.mesh.tets.len());
    let repeat = generate(
        &|_| -1.,
        Options {
            interior_refine_nodes: 20,
            ..options
        },
    )
    .unwrap();
    assert_eq!(refined.mesh.points, repeat.mesh.points);
    assert_eq!(refined.mesh.tets, repeat.mesh.tets);
    let bounded = generate(
        &|_| -1.,
        Options {
            interior_refine_nodes: 20,
            max_tetrahedra: seed.mesh.tets.len(),
            ..options
        },
    )
    .unwrap();
    assert_eq!(bounded.diagnostics.inserted_interior_nodes, 0);
    assert_eq!(bounded.mesh.points, seed.mesh.points);
    assert_eq!(bounded.mesh.tets, seed.mesh.tets);
}

#[test]
fn adaptive_interior_refinement_preserves_the_boundary() {
    let field = |p: Point| -1. + 0.05 * (TAU * p[0]).cos();
    let options = Options {
        cells: [4; 3],
        optimize_passes: 0,
        minimum_quality: 0.1,
        periodic: [true; 3],
        adapt: Some(meshers_core::tet_adapt::Options {
            iterations: 1,
            ..Default::default()
        }),
        ..Options::default()
    };
    let seed = generate(&field, options).unwrap();
    let refined = generate(
        &field,
        Options {
            interior_refine_nodes: 20,
            ..options
        },
    )
    .unwrap();
    assert!(refined.diagnostics.inserted_interior_nodes > 0);
    assert!(refined.diagnostics.inserted_interior_nodes <= 20);
    assert_eq!(
        &refined.mesh.points[..seed.mesh.points.len()],
        &seed.mesh.points
    );
    assert_eq!(refined.mesh.surface, seed.mesh.surface);
    assert_eq!(refined.boundary_tags, seed.boundary_tags);
    assert_eq!(refined.periodic_pairs, seed.periodic_pairs);
    assert!((refined.diagnostics.volume - seed.diagnostics.volume).abs() < 1e-12);
    assert!(refined.diagnostics.minimum_mmg_quality >= 0.1);
}

#[test]
fn requested_quality_allows_tighter_curvature_refinement() {
    let mesh = generate(
        &gyroid,
        Options {
            bounds: [[-0.5; 3], [0.5; 3]],
            region: Region::Band {
                lower: -0.55,
                upper: 0.55,
            },
            cells: [48; 3],
            minimum_quality: 0.1,
            geometry_tolerance: 0.0002,
            optimize_passes: 1,
            periodic: [true; 3],
            ..Options::default()
        },
    )
    .unwrap();
    assert!(mesh.diagnostics.maximum_sampled_surface_error <= 0.0002);
    assert!(mesh.diagnostics.minimum_mmg_quality >= 0.1);
    assert!(mesh.periodic_pairs.iter().all(|pairs| !pairs.is_empty()));
}

struct BenchmarkGyroid {
    variant: usize,
    batched: bool,
}
impl ScalarField for BenchmarkGyroid {
    fn value(&self, p: Point) -> f64 {
        let width = match self.variant {
            0 => 0.25,
            1 => 0.25 * (1. + 0.4 * p[0]),
            _ => 0.125,
        };
        gyroid(p) / width
    }
    fn prefers_batches(&self) -> bool {
        self.batched
    }
}

#[test]
fn quality_optimization_preserves_the_coarse_benchmark_geometry_limit() {
    for variant in 0..3 {
        for batched in [false, true] {
            for optimize_passes in [4, 12] {
                let result = generate(
                    &BenchmarkGyroid { variant, batched },
                    Options {
                        bounds: [[-0.5; 3], [0.5; 3]],
                        cells: [8; 3],
                        region: Region::Band {
                            lower: -1.,
                            upper: 1.,
                        },
                        geometry_tolerance: 0.01,
                        minimum_quality: 0.1,
                        max_tetrahedra: 250_000,
                        optimize_passes,
                        threads: 1,
                        ..Options::default()
                    },
                )
                .unwrap();
                assert!(result.diagnostics.maximum_sampled_surface_error <= 0.01);
                assert!(result.diagnostics.minimum_mmg_quality >= 0.1);
            }
        }
    }
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
            minimum_quality: 0.1,
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
            minimum_quality: 0.1,
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
#[test]
fn optimized_mesh_is_independent_of_physical_unit_scale() {
    let mut reference: Option<(Vec<[usize; 4]>, Vec<Point>, f64)> = None;
    for scale in [1., 1e-3, 1e3] {
        let field = |p: Point| {
            let p = p.map(|x| x / scale);
            ((p[0] - 0.47).powi(2) + (p[1] - 0.51).powi(2) + (p[2] - 0.49).powi(2)).sqrt() - 0.29
        };
        let output = generate(
            &field,
            Options {
                bounds: [[0.; 3], [scale; 3]],
                cells: [10; 3],
                optimize_passes: 2,
                geometry_tolerance: 0.05 * scale,
                ..Default::default()
            },
        )
        .unwrap();
        let points = output
            .mesh
            .points
            .iter()
            .map(|p| p.map(|x| x / scale))
            .collect::<Vec<_>>();
        if let Some((tets, old_points, quality)) = &reference {
            assert!(
                output.mesh.tets == *tets,
                "Unit change altered connectivity"
            );
            for (a, b) in points.iter().zip(old_points) {
                assert!(meshers_core::norm(meshers_core::sub(*a, *b)) < 1e-6);
            }
            assert!((output.diagnostics.minimum_mmg_quality - quality).abs() < 1e-6);
        } else {
            reference = Some((
                output.mesh.tets,
                points,
                output.diagnostics.minimum_mmg_quality,
            ));
        }
    }
}

#[test]
fn snapping_does_not_erase_a_sampled_small_component() {
    let field = |p: Point| {
        let large = ((p[0] + 0.35).powi(2) + p[1].powi(2) + p[2].powi(2)).sqrt() - 0.4;
        let small =
            ((p[0] - 0.65).powi(2) + (p[1] - 0.2).powi(2) + (p[2] - 0.13).powi(2)).sqrt() - 0.075;
        large.min(small)
    };
    // Isolate clipping and snapping. An FEA quality floor may subsequently
    // reject a coarse component, but it must not disappear silently.
    let output = generate(
        &field,
        Options {
            bounds: [[-1.; 3], [1.; 3]],
            cells: [12; 3],
            optimize_passes: 0,
            minimum_quality: 0.,
            geometry_tolerance: 0.2,
            ..Default::default()
        },
    )
    .unwrap();
    let mesh = &output.mesh;
    let mut parent: Vec<usize> = (0..mesh.points.len()).collect();
    fn find(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    for tet in &mesh.tets {
        for &v in &tet[1..] {
            let a = find(&mut parent, tet[0]);
            let b = find(&mut parent, v);
            parent[a] = b;
        }
    }
    let mut groups = std::collections::BTreeSet::new();
    for tet in &mesh.tets {
        groups.insert(find(&mut parent, tet[0]));
    }
    assert_eq!(groups.len(), 2, "Snapping erased an interior component");
    assert!(
        mesh.tets
            .iter()
            .any(|tet| tet.iter().all(|&i| mesh.points[i][0] > 0.5))
    );
}

/// Thirty dominant modes of a band-limited random field with features near
/// one background cell: `threshold - sum a cos(tau k.x + phase)`.
struct Spinodoid {
    threshold: f64,
    modes: Vec<[f64; 5]>,
}

impl Spinodoid {
    fn load() -> Self {
        let mut lines = include_str!("../src/testdata/spinodoid-30-modes.txt").lines();
        let threshold = lines.next().unwrap().parse().unwrap();
        let modes = lines
            .map(|line| {
                let v: Vec<f64> = line
                    .split_whitespace()
                    .map(|x| x.parse().unwrap())
                    .collect();
                [v[0], v[1], v[2], v[3], v[4]]
            })
            .collect();
        Self { threshold, modes }
    }
}

impl ScalarField for Spinodoid {
    fn value(&self, p: Point) -> f64 {
        self.threshold
            - self
                .modes
                .iter()
                .map(|m| m[3] * (TAU * (m[0] * p[0] + m[1] * p[1] + m[2] * p[2]) + m[4]).cos())
                .sum::<f64>()
    }
    fn gradient(&self, p: Point) -> Option<Point> {
        let mut gradient = [0.; 3];
        for m in &self.modes {
            let s = TAU * m[3] * (TAU * (m[0] * p[0] + m[1] * p[1] + m[2] * p[2]) + m[4]).sin();
            for a in 0..3 {
                gradient[a] += s * m[a];
            }
        }
        Some(gradient)
    }
}

#[test]
fn fully_snapped_cap_triangles_stay_periodic() {
    // Without the cap rollback: "periodic triangle mismatch on axis 1".
    let output = generate(
        &Spinodoid::load(),
        Options {
            periodic: [true; 3],
            geometry_tolerance: 10.,
            optimize_passes: 0,
            ..Options::default()
        },
    )
    .unwrap();
    assert!(output.diagnostics.reverted_snap_vertices > 0);
    assert!(output.periodic_pairs.iter().all(|pairs| !pairs.is_empty()));
}

#[test]
fn geometry_gate_bounds_the_estimate_near_critical_points() {
    // The first-order estimate alone reports 0.22 here, above half a cell
    // diagonal, at a facet sample where the gradient nearly vanishes.
    let output = generate(
        &Spinodoid::load(),
        Options {
            cells: [16; 3],
            periodic: [true; 3],
            geometry_tolerance: 0.06,
            optimize_passes: 1,
            ..Options::default()
        },
    )
    .unwrap();
    assert!(output.diagnostics.maximum_sampled_surface_error <= 0.06);
}

#[test]
fn rejected_slivers_are_remeshed_without_breaking_periodic_caps() {
    let field = Spinodoid::load();
    let options = Options {
        periodic: [true; 3],
        geometry_tolerance: 0.1,
        ..Options::default()
    };
    let plain = generate(&field, options).unwrap();
    let gated = generate(
        &field,
        Options {
            minimum_quality: 0.1,
            ..options
        },
    )
    .unwrap();
    assert!(plain.diagnostics.minimum_mmg_quality < 0.1);
    assert!(gated.diagnostics.minimum_mmg_quality >= 0.1);
    // Generation itself rejects caps that no longer match.
    assert!(gated.periodic_pairs.iter().all(|pairs| !pairs.is_empty()));
}

#[test]
fn fibre_grazing_a_box_edge_is_defeatured_within_tolerance() {
    // A periodic fibre along x whose wall passes 0.004 from four box edges
    // leaves there a sliver of void that no conforming element can fill.
    let centre = 0.5 - 0.2 / 2_f64.sqrt() + 0.004;
    let fibre = move |p: Point| {
        let mut distance = f64::INFINITY;
        for dy in [-1., 0., 1.] {
            for dz in [-1., 0., 1.] {
                distance = distance.min((p[1] - centre - dy).hypot(p[2] - centre - dz));
            }
        }
        distance - 0.2
    };
    let options = Options {
        bounds: [[-0.5; 3], [0.5; 3]],
        cells: [16; 3],
        periodic: [true; 3],
        geometry_tolerance: 0.02,
        ..Options::default()
    };
    let plain = generate(&fibre, options).unwrap();
    let gated = generate(
        &fibre,
        Options {
            minimum_quality: 0.1,
            ..options
        },
    )
    .unwrap();
    assert!(plain.diagnostics.minimum_mmg_quality < 0.1);
    assert_eq!(plain.diagnostics.defeaturing_moves, 0);
    assert!(gated.diagnostics.minimum_mmg_quality >= 0.1);
    assert!(gated.diagnostics.defeaturing_moves > 0);
    assert!(gated.diagnostics.maximum_sampled_surface_error <= 0.02);
    // Merged crossings leave fewer paired nodes; generation itself rejects
    // caps that no longer match.
    assert!(gated.periodic_pairs.iter().all(|pairs| !pairs.is_empty()));
}

#[test]
fn a_gap_narrower_than_a_cell_is_not_bridged() {
    // Two spheres 0.03 apart: at 12 cells lattice nodes on both sides of the
    // gap are inside, and the elements between them would weld the spheres.
    let spheres = |p: Point| {
        let distance = |c: f64| (p[0] - c).hypot(p[1]).hypot(p[2]) - 0.2;
        distance(-0.215).min(distance(0.215))
    };
    let output = generate(
        &spheres,
        Options {
            bounds: [[-0.5; 3], [0.5; 3]],
            cells: [12; 3],
            geometry_tolerance: 0.01,
            minimum_quality: 0.1,
            ..Options::default()
        },
    )
    .unwrap();
    assert!(output.diagnostics.refined_lattice_edges > 0);
    assert!(output.diagnostics.minimum_mmg_quality >= 0.1);
    assert!(output.diagnostics.maximum_sampled_surface_error <= 0.01);
    // No element is left with a vertex in each sphere.
    assert!(output.mesh.tets.iter().all(|t| {
        t.iter().all(|&i| output.mesh.points[i][0] <= 0.)
            || t.iter().all(|&i| output.mesh.points[i][0] >= 0.)
    }));
}

#[test]
fn a_film_thinner_than_a_cell_is_meshed() {
    // The matrix around two spheres 0.03 apart keeps a film between them
    // that no node of a 12-cell lattice samples.
    let matrix = |p: Point| {
        let distance = |c: f64| (p[0] - c).hypot(p[1]).hypot(p[2]) - 0.2;
        -distance(-0.215).min(distance(0.215))
    };
    let output = generate(
        &matrix,
        Options {
            bounds: [[-0.5; 3], [0.5; 3]],
            cells: [12; 3],
            geometry_tolerance: 0.01,
            minimum_quality: 0.1,
            ..Options::default()
        },
    )
    .unwrap();
    assert!(output.diagnostics.refined_lattice_edges > 0);
    assert!(output.diagnostics.maximum_sampled_surface_error <= 0.01);
    assert!(output.diagnostics.minimum_mmg_quality >= 0.1);
    let in_film = |t: &[usize; 4]| {
        let c: Point =
            std::array::from_fn(|a| t.iter().map(|&i| output.mesh.points[i][a]).sum::<f64>() / 4.);
        c[0].abs() < 0.015 && c[1].hypot(c[2]) < 0.05
    };
    assert!(output.mesh.tets.iter().any(in_film));
}

#[test]
fn adaptation_keeps_a_periodic_spinodoid_unfoldable() {
    // Without the periodic guards of the adaptive operators, a collapse lays a
    // wall facet in a periodic plane or a wall vertex crosses one: "solid
    // boundary is not a closed two-manifold" or "periodic triangle mismatch".
    let tolerance = 0.05;
    let output = generate(
        &Spinodoid::load(),
        Options {
            cells: [24; 3],
            periodic: [true; 3],
            // Adaptation closes the geometry to its own estimate, not to the gate.
            geometry_tolerance: 1.2 * tolerance,
            minimum_quality: 0.1,
            adapt: Some(meshers_core::tet_adapt::Options {
                tolerance,
                min_length: tolerance,
                ..Default::default()
            }),
            ..Options::default()
        },
    )
    .unwrap();
    assert!(!output.diagnostics.adaptation_abandoned);
    assert!(output.diagnostics.minimum_mmg_quality >= 0.1);
    assert!(output.diagnostics.maximum_sampled_surface_error <= 1.2 * tolerance);
    assert!(output.periodic_pairs.iter().all(|pairs| !pairs.is_empty()));
}

#[test]
fn periodic_images_share_their_in_plane_coordinates_bit_for_bit() {
    // Images snapped onto the wall separately used to agree only to rounding.
    let gyroid = |p: [f64; 3]| {
        let [x, y, z] = p.map(|v| TAU * v);
        x.sin() * y.cos() + y.sin() * z.cos() + z.sin() * x.cos() - 0.3
    };
    let output = generate(
        &gyroid,
        Options {
            bounds: [[-0.5; 3], [0.5; 3]],
            cells: [16; 3],
            periodic: [true; 3],
            geometry_tolerance: 0.02,
            minimum_quality: 0.1,
            order: 2,
            ..Options::default()
        },
    )
    .unwrap();
    let quadratic = output.quadratic.as_ref().unwrap();
    for (p, pairs) in [
        (&output.mesh.points, &output.periodic_pairs),
        (&quadratic.elements.points, &quadratic.periodic_pairs),
    ] {
        for (axis, pairs) in pairs.iter().enumerate() {
            assert!(!pairs.is_empty());
            for &[i, j] in pairs {
                assert_eq!([p[i][axis], p[j][axis]], [-0.5, 0.5]);
                for a in (0..3).filter(|&a| a != axis) {
                    assert_eq!(
                        p[i][a].to_bits(),
                        p[j][a].to_bits(),
                        "axis {axis}: {i}, {j}"
                    );
                }
            }
        }
    }
}

#[test]
fn coverage_error_measures_the_wall_against_the_mesh() {
    // Wall points sampled on lines finer than the lattice, against the meshed
    // wall: the direction facet samples cannot see. On a sphere both
    // directions measure the same chord error, so they must agree in scale.
    let ball = |p: [f64; 3]| (p[0].powi(2) + p[1].powi(2) + p[2].powi(2)).sqrt() - 0.35;
    let options = Options {
        bounds: [[-0.5; 3], [0.5; 3]],
        cells: [16; 3],
        geometry_tolerance: 0.03,
        ..Options::default()
    };
    let fine = generate(&ball, options).unwrap().diagnostics;
    let coarse = generate(
        &ball,
        Options {
            adapt: Some(meshers_core::tet_adapt::Options {
                tolerance: 0.03,
                min_length: 0.03,
                ..Default::default()
            }),
            ..options
        },
    )
    .unwrap()
    .diagnostics;
    assert!(!coarse.adaptation_abandoned);
    for d in [&fine, &coarse] {
        let (surface, coverage) = (
            d.maximum_sampled_surface_error,
            d.maximum_sampled_coverage_error,
        );
        assert!(
            coverage > 0.3 * surface && coverage < 2. * surface,
            "{surface} {coverage}"
        );
    }
    assert!(coarse.maximum_sampled_coverage_error > 3. * fine.maximum_sampled_coverage_error);
}

#[test]
fn adaptation_recovers_fibres_grazing_the_box() {
    // Adaptation alone leaves elements pinned between the wall and a box
    // plane below any useful quality. Its result goes through the sliver
    // recovery of the lattice path, kept when the geometry gate accepts it.
    let graze = 0.5 - 0.2 / 2_f64.sqrt() + 0.004;
    for (name, centre) in [
        ("edge", [graze, graze]),
        // The fibre passes the face z = 0.5 by 0.003.
        ("face", [0., 0.5 - 0.2 + 0.003]),
    ] {
        let fibre = move |p: Point| {
            let mut distance = f64::INFINITY;
            for dy in [-1., 0., 1.] {
                for dz in [-1., 0., 1.] {
                    distance = distance.min((p[1] - centre[0] - dy).hypot(p[2] - centre[1] - dz));
                }
            }
            distance - 0.2
        };
        let tolerance = 0.02;
        let output = generate(
            &fibre,
            Options {
                bounds: [[-0.5; 3], [0.5; 3]],
                cells: [16; 3],
                periodic: [true; 3],
                geometry_tolerance: tolerance,
                minimum_quality: 0.1,
                adapt: Some(meshers_core::tet_adapt::Options {
                    tolerance,
                    min_length: tolerance,
                    ..Default::default()
                }),
                ..Options::default()
            },
        )
        .unwrap_or_else(|error| panic!("{name}: {error}"));
        let d = output.diagnostics;
        assert!(!d.adaptation_abandoned, "{name}");
        assert!(d.minimum_mmg_quality >= 0.1, "{name}");
        assert!(d.defeaturing_moves > 0, "{name}");
        assert!(d.maximum_sampled_surface_error <= tolerance, "{name}");
        assert!(d.maximum_sampled_coverage_error <= tolerance, "{name}");
        // The fibre runs along x: its caps on the x faces always match. A
        // thin strip across another face may be defeatured away.
        assert!(!output.periodic_pairs[0].is_empty(), "{name}");
    }
}

#[test]
fn a_refused_adaptation_falls_back_to_the_lattice_mesh() {
    // Adaptation works to a looser tolerance than the final checks accept and
    // coarsens the wall beyond it. The caller gets the lattice mesh of the same
    // field, flagged as such, instead of a rejection.
    let ball = |p: [f64; 3]| (p[0].powi(2) + p[1].powi(2) + p[2].powi(2)).sqrt() - 0.35;
    let tolerance = 0.004;
    let options = Options {
        bounds: [[-0.5; 3], [0.5; 3]],
        cells: [16; 3],
        periodic: [true; 3],
        geometry_tolerance: tolerance,
        minimum_quality: 0.1,
        ..Options::default()
    };
    let lattice = generate(&ball, options).unwrap();
    let adapted = generate(
        &ball,
        Options {
            adapt: Some(meshers_core::tet_adapt::Options {
                tolerance: 0.05,
                min_length: 0.1,
                max_length: 0.25,
                ..Default::default()
            }),
            ..options
        },
    )
    .unwrap();
    assert!(!lattice.diagnostics.adaptation_abandoned);
    assert!(adapted.diagnostics.adaptation_abandoned);
    assert!(adapted.diagnostics.adapt.is_none());
    assert_eq!(adapted.mesh.points, lattice.mesh.points);
    assert_eq!(adapted.mesh.tets, lattice.mesh.tets);
    assert!(adapted.diagnostics.minimum_mmg_quality >= 0.1);
    assert!(adapted.diagnostics.maximum_sampled_surface_error <= tolerance);
}
