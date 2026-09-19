use meshers_core::{
    Point,
    implicit::{self, intersection},
};
fn options() -> implicit::Options {
    implicit::Options {
        cells: [4; 3],
        snap: 0.,
        optimize_passes: 0,
        geometry_tolerance: 0.,
        periodic: [false, false, true],
        ..Default::default()
    }
}
#[test]
fn simultaneous_planes_preserve_their_intersection() {
    let a = |p: Point| p[0] - 0.61;
    let b = |p: Point| p[1] - 0.73;
    let result = intersection::generate(&[&a, &b], |p| p, options()).unwrap();
    assert!(result.constraints.contains(&3));
    for (&p, &mask) in result.mesh.points.iter().zip(&result.constraints) {
        assert!(p[0] <= 0.61 + 1e-12 && p[1] <= 0.73 + 1e-12);
        if mask & 1 != 0 {
            assert!((p[0] - 0.61).abs() < 1e-12);
        }
        if mask & 2 != 0 {
            assert!((p[1] - 0.73).abs() < 1e-12);
        }
    }
    let volume: f64 = result
        .mesh
        .tets
        .iter()
        .map(|t| meshers_core::determinant(t.map(|i| result.mesh.points[i])) / 6.)
        .sum();
    assert!((volume - 0.61 * 0.73).abs() < 1e-12);
}
#[test]
fn inverted_map_and_unsupported_options_are_rejected() {
    let a = |p: Point| p[0] - 0.61;
    let b = |p: Point| p[1] - 0.73;
    assert!(
        intersection::generate(&[&a, &b], |[x, y, z]| [-x, y, z], options())
            .err()
            .unwrap()
            .contains("inverted")
    );
    assert!(
        intersection::generate(
            &[&a, &b],
            |p| p,
            implicit::Options {
                optimize_passes: 21,
                ..options()
            }
        )
        .is_err()
    );
}

#[test]
fn checks_individual_periodic_constraints_even_when_one_is_inactive() {
    let inside = |_: Point| -2.;
    let nonperiodic = |p: Point| -3. - p[2];
    assert!(intersection::generate(&[&inside, &nonperiodic], |p| p, options()).is_err());
}

#[test]
fn joint_snapping_preserves_both_planes_and_periodic_end_triangles() {
    use std::collections::BTreeSet;
    let a = |p: Point| p[0] - 0.51;
    let b = |p: Point| p[1] - 0.74;
    let result = intersection::generate(
        &[&a, &b],
        |p| p,
        implicit::Options {
            snap: 0.2,
            ..options()
        },
    )
    .unwrap();
    assert!(result.constraints.contains(&3));
    for (&p, &mask) in result.parameters.iter().zip(&result.constraints) {
        if mask == 3 {
            assert!((p[0] - 0.51).abs() < 1e-10 && (p[1] - 0.74).abs() < 1e-10);
        }
    }
    let faces = |z: f64| {
        result
            .mesh
            .surface
            .iter()
            .filter(|f| f.iter().all(|&i| result.parameters[i][2] == z))
            .map(|f| {
                let mut keys = f.map(|i| {
                    [result.parameters[i][0], result.parameters[i][1]]
                        .map(|v| (v * 1e10).round() as i64)
                });
                keys.sort_unstable();
                keys
            })
            .collect::<BTreeSet<_>>()
    };
    assert!(!faces(0.).is_empty());
    assert_eq!(faces(0.), faces(1.));
}

#[test]
fn curved_background_retains_rotation_periodicity_and_positive_elements() {
    let a = |p: Point| p[0] - 0.61;
    let b = |p: Point| p[1] - 0.73;
    let angle = std::f64::consts::PI / 3.;
    let map = |[x, y, z]: Point| {
        [
            (2. + x) * (angle * z).cos(),
            (2. + x) * (angle * z).sin(),
            -y,
        ]
    };
    let result = intersection::generate(&[&a, &b], map, options()).unwrap();
    for t in &result.mesh.tets {
        assert!(meshers_core::determinant(t.map(|i| result.mesh.points[i])) > 0.);
    }
    let end: Vec<_> = result
        .parameters
        .iter()
        .enumerate()
        .filter(|(_, p)| p[2] == 1.)
        .collect();
    for (i, p) in result
        .parameters
        .iter()
        .enumerate()
        .filter(|(_, p)| p[2] == 0.)
    {
        let (j, _) = end
            .iter()
            .find(|(_, q)| (p[0] - q[0]).abs() < 1e-10 && (p[1] - q[1]).abs() < 1e-10)
            .unwrap();
        let [x, y, z] = result.mesh.points[i];
        let expected = [
            x * angle.cos() - y * angle.sin(),
            x * angle.sin() + y * angle.cos(),
            z,
        ];
        assert!(meshers_core::norm(meshers_core::sub(expected, result.mesh.points[*j])) < 1e-10);
    }
}

#[test]
fn nonfinite_fields_and_output_budget_are_rejected() {
    let bad = |_: Point| f64::NAN;
    let inside = |_: Point| -1.;
    assert!(intersection::generate(&[&bad, &inside], |p| p, options()).is_err());
    assert!(
        intersection::generate(
            &[&inside, &inside],
            |p| p,
            implicit::Options {
                max_tetrahedra: 1,
                ..options()
            }
        )
        .err()
        .unwrap()
        .contains("budget")
    );
}

#[test]
fn quality_improvement_keeps_plane_geometry_volume_and_periodic_faces() {
    use std::collections::BTreeSet;
    let a = |p: Point| p[0] - 0.501;
    let b = |p: Point| p[1] - 0.751;
    let o = implicit::Options {
        cells: [8; 3],
        ..options()
    };
    let baseline = intersection::generate(&[&a, &b], |p| p, o).unwrap();
    let improved = intersection::generate(
        &[&a, &b],
        |p| p,
        implicit::Options {
            optimize_passes: 4,
            ..o
        },
    )
    .unwrap();
    let quality = |r: &intersection::Output| {
        r.mesh
            .tets
            .iter()
            .map(|t| meshers_core::quality(t.map(|i| r.mesh.points[i])).powf(1.5))
            .fold(1., f64::min)
    };
    assert!(quality(&improved) > 10. * quality(&baseline));
    let volume: f64 = improved
        .mesh
        .tets
        .iter()
        .map(|t| meshers_core::determinant(t.map(|i| improved.mesh.points[i])) / 6.)
        .sum();
    assert!((volume - 0.501 * 0.751).abs() < 1e-10);
    for (&p, &mask) in improved.mesh.points.iter().zip(&improved.constraints) {
        if mask & 1 != 0 {
            assert!((p[0] - 0.501).abs() < 1e-10);
        }
        if mask & 2 != 0 {
            assert!((p[1] - 0.751).abs() < 1e-10);
        }
    }
    let faces = |z: f64| {
        improved
            .mesh
            .surface
            .iter()
            .filter(|f| f.iter().all(|&i| improved.parameters[i][2] == z))
            .map(|f| {
                let mut key = f.map(|i| {
                    [improved.parameters[i][0], improved.parameters[i][1]]
                        .map(|v| (v * 1e10).round() as i64)
                });
                key.sort_unstable();
                key
            })
            .collect::<BTreeSet<_>>()
    };
    assert_eq!(faces(0.), faces(1.));
}

#[test]
fn quality_gate_rejects_unmet_floor_and_invalid_thresholds() {
    let a = |p: Point| p[0] - 0.501;
    let b = |p: Point| p[1] - 0.751;
    for threshold in [f64::NAN, -0.1, 1.1] {
        assert!(
            intersection::generate(
                &[&a, &b],
                |p| p,
                implicit::Options {
                    minimum_quality: threshold,
                    ..options()
                }
            )
            .is_err()
        );
    }
    let error = intersection::generate(
        &[&a, &b],
        |p| p,
        implicit::Options {
            minimum_quality: 1.,
            optimize_passes: 2,
            ..options()
        },
    )
    .err()
    .unwrap();
    assert!(error.contains("minimum MMG quality"));
}

#[test]
fn roundoff_at_cylinder_grid_vertices_does_not_create_duplicate_slivers() {
    let cylinder = |p: Point| p[0].hypot(p[1]) - 1.5;
    let plane = |p: Point| p[0] - 0.9;
    let result = intersection::generate(
        &[&cylinder, &plane],
        |p| p,
        implicit::Options {
            bounds: [[-1.5, -1.5, 0.], [1.5, 1.5, 1.]],
            cells: [40, 40, 4],
            ..options()
        },
    )
    .unwrap();
    assert!(!result.mesh.tets.is_empty());
    assert!(
        result
            .mesh
            .tets
            .iter()
            .all(|t| meshers_core::determinant(t.map(|i| result.mesh.points[i])) > 0.)
    );
    for (&p, &mask) in result.mesh.points.iter().zip(&result.constraints) {
        if mask == 3 {
            assert!((p[0] - 0.9).abs() < 1e-12);
            assert!((p[0].hypot(p[1]) - 1.5).abs() < 1e-12);
        }
    }
}

#[test]
fn rescaled_constraints_preserve_snapped_planes_and_volume() {
    for scale in [1e-14, 1e-9, 1., 1e9, 1e14] {
        for optimize_passes in [0, 2] {
            let a = |p: Point| scale * (p[0] - 0.51);
            let b = |p: Point| (p[1] - 0.74) / scale;
            let result = intersection::generate(
                &[&a, &b],
                |p| p,
                implicit::Options {
                    snap: 0.2,
                    optimize_passes,
                    ..options()
                },
            )
            .unwrap();
            assert!(result.constraints.contains(&3));
            for (&p, &mask) in result.mesh.points.iter().zip(&result.constraints) {
                assert!(p[0] <= 0.51 + 1e-11 && p[1] <= 0.74 + 1e-11);
                if mask & 1 != 0 {
                    assert!((p[0] - 0.51).abs() < 1e-11, "scale={scale}, point={p:?}");
                }
                if mask & 2 != 0 {
                    assert!((p[1] - 0.74).abs() < 1e-11, "scale={scale}, point={p:?}");
                }
            }
            let volume: f64 = result
                .mesh
                .tets
                .iter()
                .map(|t| {
                    let d = meshers_core::determinant(t.map(|i| result.mesh.points[i]));
                    assert!(d > 0.);
                    d / 6.
                })
                .sum();
            assert!(
                (volume - 0.51 * 0.74).abs() < 1e-11,
                "scale={scale}, volume={volume}"
            );
        }
    }
}
