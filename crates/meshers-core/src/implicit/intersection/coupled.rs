//! Bounded coupled corner repair. Commit only an entirely passing incident cavity.
use super::{Output, ScalarField};
use crate::implicit::Options;
use crate::{Point, determinant, quality};
use slsqp_rssl::{Constraint, Slsqp, SlsqpMode};
use std::collections::{BTreeSet, HashMap};

pub(super) fn improve(
    output: &mut Output,
    o: &Options,
    fields: &[&dyn ScalarField],
    map: &dyn Fn(Point) -> Point,
) -> Result<usize, String> {
    if o.minimum_quality == 0. || o.optimize_passes == 0 || o.periodic.iter().any(|&x| x) {
        return Ok(0);
    }
    let mut bad: Vec<_> = (0..output.mesh.tets.len())
        .filter(|&i| {
            quality(output.mesh.tets[i].map(|v| output.mesh.points[v])).powf(1.5)
                < o.minimum_quality
        })
        .collect();
    bad.sort_by(|&a, &b| {
        quality(output.mesh.tets[a].map(|v| output.mesh.points[v]))
            .total_cmp(&quality(output.mesh.tets[b].map(|v| output.mesh.points[v])))
    });
    let neighborhoods: Vec<Vec<usize>> = bad
        .into_iter()
        .take(32)
        .map(|ti| output.mesh.tets[ti].to_vec())
        .collect();
    project_neighborhoods(output, o, fields, map, &neighborhoods, None)
}

pub(crate) fn project_neighborhoods(
    output: &mut Output,
    o: &Options,
    fields: &[&dyn ScalarField],
    map: &dyn Fn(Point) -> Point,
    neighborhoods: &[Vec<usize>],
    geometry: Option<f64>,
) -> Result<usize, String> {
    let target = (o.minimum_quality + 2e-7).min(1.);
    let cell: Point =
        std::array::from_fn(|a| (o.bounds[1][a] - o.bounds[0][a]) / o.cells[a] as f64);
    let h = cell.into_iter().fold(f64::INFINITY, f64::min);
    let mut incident = vec![Vec::new(); output.mesh.points.len()];
    for (i, t) in output.mesh.tets.iter().enumerate() {
        for &v in t {
            incident[v].push(i);
        }
    }
    let mut moved = 0;
    for vertices in neighborhoods.iter().take(32) {
        for f in fields {
            f.check()?;
        }
        let mut near: BTreeSet<_> = vertices.iter().copied().collect();
        for _ in 0..2 {
            near = near
                .iter()
                .flat_map(|&v| incident[v].iter().flat_map(|&i| output.mesh.tets[i]))
                .collect();
        }
        // ponytail: bound dense local work; larger cavities need a sparse repair.
        if near.len() > 64 {
            continue;
        }
        let near: Vec<_> = near.into_iter().collect();
        let ids: BTreeSet<_> = near
            .iter()
            .flat_map(|&v| incident[v].iter().copied())
            .collect();
        if ids.len() > 512 {
            continue;
        }
        let lookup: HashMap<_, _> = near.iter().enumerate().map(|(j, &v)| (v, j)).collect();
        let guarded: Vec<_> = if geometry.is_some() {
            output
                .mesh
                .surface
                .iter()
                .filter(|f| f.iter().any(|v| lookup.contains_key(v)))
                .map(|&f| {
                    let active = f.iter().fold(super::ConstraintMask::MAX, |mask, &v| {
                        mask & output.constraints[v]
                    });
                    (f, active)
                })
                .collect()
        } else {
            Vec::new()
        };
        let origins: Vec<_> = near.iter().map(|&v| output.parameters[v]).collect();
        let active: Vec<Vec<usize>> = near
            .iter()
            .map(|&v| {
                (0..fields.len())
                    .filter(|&j| output.constraints[v] & (1 << j) != 0)
                    .collect()
            })
            .collect();
        let locked: Vec<[bool; 3]> = origins
            .iter()
            .map(|p| std::array::from_fn(|a| p[a] == o.bounds[0][a] || p[a] == o.bounds[1][a]))
            .collect();
        let retract = |a: &[Point]| -> Option<Vec<Point>> {
            a.iter()
                .enumerate()
                .map(|(j, &p)| {
                    let mask = active[j].iter().fold(0u64, |mask, &i| mask | (1u64 << i));
                    let mut candidate = p;
                    for a in 0..3 {
                        if locked[j][a] {
                            candidate[a] = origins[j][a];
                        }
                    }
                    let p = crate::projection::onto_constraints(
                        fields, mask, candidate, locked[j], o.bounds,
                    )?;
                    if (0..3).any(|a| p[a] < o.bounds[0][a] || p[a] > o.bounds[1][a])
                        || fields.iter().any(|&f| {
                            let value = f.value(p);
                            value > 0.
                                && value
                                    > 1e-10 * crate::norm(super::field_gradient(f, p, h * 1e-5))
                        })
                    {
                        None
                    } else {
                        Some(p)
                    }
                })
                .collect()
        };
        let physical = |a: &[Point]| a.iter().map(|&p| map(p)).collect::<Vec<_>>();
        let passing = |phys: &[Point]| {
            ids.iter().all(|&i| {
                let tet = output.mesh.tets[i]
                    .map(|v| lookup.get(&v).map_or(output.mesh.points[v], |&j| phys[j]));
                let det = determinant(tet);
                det.is_finite() && det > 0. && quality(tet).powf(1.5) >= o.minimum_quality
            }) && guarded.iter().all(|&(face, active)| {
                let p = face.map(|v| lookup.get(&v).map_or(output.mesh.points[v], |&j| phys[j]));
                crate::projection::facet_distance(fields, u64::from(active), p, o.bounds)
                    <= geometry.unwrap()
            })
        };
        let initial = retract(&origins).filter(|a| passing(&physical(a)));
        let a = if let Some(a) = initial {
            a
        } else {
            let mut variables = vec![[None; 3]; near.len()];
            let mut bounds = Vec::new();
            for j in 0..near.len() {
                for axis in 0..3 {
                    if !locked[j][axis] {
                        variables[j][axis] = Some(bounds.len());
                        bounds.push((
                            (-1_f64).max((o.bounds[0][axis] - origins[j][axis]) / cell[axis]),
                            1_f64.min((o.bounds[1][axis] - origins[j][axis]) / cell[axis]),
                        ));
                    }
                }
            }
            if bounds.is_empty() {
                continue;
            }
            let point = |j: usize, x: &[f64]| -> Point {
                std::array::from_fn(|axis| {
                    origins[j][axis] + variables[j][axis].map_or(0., |k| x[k] * cell[axis])
                })
            };
            let point = &point;
            let mut constraints: Vec<Constraint<'_>> = Vec::new();
            for &i in &ids {
                let tet = output.mesh.tets[i];
                let lookup = &lookup;
                let reference = &output.mesh.points;
                constraints.push(Constraint::Ineq(Box::new(move |x| {
                    let p = tet.map(|v| lookup.get(&v).map_or(reference[v], |&j| map(point(j, x))));
                    let q = quality(p).powf(1.5) * determinant(p).signum();
                    if q.is_finite() { q - target } else { -1. }
                })));
            }
            for &(face, active) in &guarded {
                let lookup = &lookup;
                let reference = &output.mesh.points;
                constraints.push(Constraint::Ineq(Box::new(move |x| {
                    let p =
                        face.map(|v| lookup.get(&v).map_or(reference[v], |&j| map(point(j, x))));
                    let error =
                        crate::projection::facet_distance(fields, u64::from(active), p, o.bounds);
                    // Optimize inside the gate; solver residuals do not admit a mesh.
                    if error.is_finite() {
                        0.9999 - error / geometry.unwrap()
                    } else {
                        -1.
                    }
                })));
            }
            for j in 0..near.len() {
                let mut basis: Vec<Point> = Vec::new();
                for (i, &field) in fields.iter().enumerate() {
                    let gradient = super::field_gradient(field, origins[j], h * 1e-5);
                    let length = crate::norm(gradient);
                    let scale = if length.is_finite() && length > 0. {
                        length
                    } else {
                        let value = field.value(origins[j]).abs();
                        if value > 0. && value.is_finite() {
                            value
                        } else {
                            1.
                        }
                    };
                    let is_active = active[j].contains(&i);
                    let mut normal = gradient.map(|g| g / scale);
                    for (value, locked) in normal.iter_mut().zip(locked[j]) {
                        if locked {
                            *value = 0.;
                        }
                    }
                    for &n in &basis {
                        normal = crate::sub(normal, n.map(|v| v * crate::dot(normal, n)));
                    }
                    let rank = crate::norm(normal);
                    if is_active && rank.is_finite() && rank > 1e-10 {
                        basis.push(normal.map(|v| v / rank));
                        constraints.push(Constraint::Eq(Box::new(move |x| {
                            field.value(point(j, x)) / scale
                        })));
                    } else {
                        constraints.push(Constraint::Ineq(Box::new(move |x| {
                            -field.value(point(j, x)) / scale
                        })));
                        if is_active {
                            constraints.push(Constraint::Ineq(Box::new(move |x| {
                                field.value(point(j, x)) / scale
                            })));
                        }
                    }
                }
            }
            let iterations = 500;
            let mut solver = Slsqp::new_with_observer(
                vec![0.; bounds.len()],
                &bounds,
                Box::new(|x| 0.5 * x.iter().map(|v| v * v).sum::<f64>()),
                constraints,
                iterations,
                1e-12,
                (),
            );
            for _ in 0..(3 * iterations + 1) {
                for field in fields {
                    field.check()?;
                }
                match solver.step() {
                    Ok(SlsqpMode::Success) | Err(_) => break,
                    Ok(_) => {}
                }
            }
            for field in fields {
                field.check()?;
            }
            let parameters: Vec<_> = (0..near.len()).map(|j| point(j, &solver.x)).collect();
            let Some(a) = retract(&parameters) else {
                continue;
            };
            a
        };
        let phys = physical(&a);
        if passing(&phys) {
            for (j, &v) in near.iter().enumerate() {
                moved +=
                    usize::from(output.parameters[v] != a[j] || output.mesh.points[v] != phys[j]);
                output.parameters[v] = a[j];
                output.mesh.points[v] = phys[j];
            }
        }
    }
    Ok(moved)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn coupled_corner_motion_preserves_every_incident_quality_gate() {
        check_corner(
            include_str!("../../testdata/annular-corner-star.json"),
            [24, 12, 128],
        );
    }
    #[test]
    fn coupled_coarse_corner_motion_preserves_every_incident_quality_gate() {
        check_corner(
            include_str!("../../testdata/annular-coarse-corner-star.json"),
            [16, 8, 96],
        );
    }

    fn check_corner(input: &str, cells: [usize; 3]) {
        for scale in [1e-12, 1., 1e12] {
            check_scaled_corner(input, cells, scale);
        }
    }
    fn check_scaled_corner(input: &str, cells: [usize; 3], scale: f64) {
        #[derive(serde::Deserialize)]
        struct Star {
            points: Vec<Point>,
            parameters: Vec<Point>,
            tets: Vec<crate::Tet>,
            constraints: Vec<super::super::ConstraintMask>,
        }
        let seed: Star = serde_json::from_str(input).unwrap();
        let gyroid = |p: Point| {
            let a = p.map(|v| v * std::f64::consts::TAU);
            a[0].sin() * a[1].cos() + a[1].sin() * a[2].cos() + a[2].sin() * a[0].cos()
        };
        let upper = |p: Point| scale * (gyroid(p) - 0.55);
        let lower = |p: Point| scale * (-gyroid(p) - 0.55);
        let map = |p: Point| {
            let angle = std::f64::consts::TAU * (p[2] + 0.5);
            let r = 0.3 + 0.2 * p[1];
            [r * angle.cos(), r * angle.sin(), p[0] + 0.5]
        };
        let mut out = Output {
            mesh: crate::Mesh {
                points: seed.points,
                tets: seed.tets,
                surface: Vec::new(),
            },
            parameters: seed.parameters,
            constraints: seed.constraints,
            quality: Default::default(),
        };
        let original_points = out.mesh.points.clone();
        let original_parameters = out.parameters.clone();
        let original_masks = out.constraints.clone();
        let original_tets = out.mesh.tets.clone();
        let old = out
            .mesh
            .tets
            .iter()
            .map(|t| quality(t.map(|v| out.mesh.points[v])).powf(1.5))
            .fold(1., f64::min);
        assert!(old < 0.1);
        let o = Options {
            bounds: [[-0.5, -0.5, -1.], [0.5, 0.5, 1.]],
            cells,
            minimum_quality: 0.1,
            ..Default::default()
        };
        assert_eq!(
            improve(
                &mut out,
                &Options {
                    optimize_passes: 0,
                    ..o
                },
                &[&upper, &lower],
                &map
            )
            .unwrap(),
            0
        );
        assert_eq!(out.mesh.points, original_points);
        assert_eq!(out.parameters, original_parameters);
        assert!(
            improve(&mut out, &o, &[&upper, &lower], &map).unwrap() > 0,
            "corner grid {cells:?}, initial quality {old}"
        );
        assert!(
            out.mesh
                .tets
                .iter()
                .all(|t| determinant(t.map(|v| out.mesh.points[v])) > 0.
                    && quality(t.map(|v| out.mesh.points[v])).powf(1.5) >= 0.1)
        );
        assert_eq!(out.constraints, original_masks);
        assert_eq!(out.mesh.tets, original_tets);
        for (i, &parameter) in out.parameters.iter().enumerate() {
            for axis in 0..3 {
                if original_parameters[i][axis] == o.bounds[0][axis]
                    || original_parameters[i][axis] == o.bounds[1][axis]
                {
                    assert_eq!(parameter[axis], original_parameters[i][axis]);
                }
            }
            for (bit, field) in [&upper as &dyn ScalarField, &lower as &dyn ScalarField]
                .into_iter()
                .enumerate()
            {
                if out.constraints[i] & (1 << bit) != 0 {
                    let before = field.value(original_parameters[i]).abs() / scale;
                    let after = field.value(parameter).abs() / scale;
                    assert!(after <= before.max(1e-9));
                    if parameter != original_parameters[i]
                        || out.mesh.points[i] != original_points[i]
                    {
                        assert!(after < 1e-9);
                    }
                }
            }
            if parameter != original_parameters[i] || out.mesh.points[i] != original_points[i] {
                assert!(crate::norm(crate::sub(map(parameter), out.mesh.points[i])) < 1e-12);
            }
        }
    }
    #[test]
    fn already_passing_cavity_does_not_report_vertex_moves() {
        let points = vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
        let mut output = Output {
            mesh: crate::Mesh {
                points: points.clone(),
                tets: vec![[0, 1, 2, 3]],
                surface: vec![],
            },
            parameters: points.clone(),
            constraints: vec![0; 4],
            quality: Default::default(),
        };
        let o = Options {
            minimum_quality: 0.1,
            ..Default::default()
        };
        let field = |_: Point| -1.;
        assert_eq!(
            project_neighborhoods(&mut output, &o, &[&field], &|p| p, &[vec![0]], None).unwrap(),
            0
        );
        assert_eq!(output.mesh.points, points);
    }
}
