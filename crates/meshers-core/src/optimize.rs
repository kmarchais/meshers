//! Local vertex optimization, with translated boundary nodes moved as one group.
#[path = "optimize_batch.rs"]
mod batch;
use crate::accelerator::Group;
use crate::{Mesh, Point, determinant, dot, norm, sub};
use rayon::prelude::*;
use std::collections::BTreeMap;

fn project(
    mut p: Point,
    original: Point,
    level: Option<f64>,
    geometry: &impl crate::geometry::Geometry,
) -> Option<Point> {
    let locked = geometry.locked(original);
    for a in 0..3 {
        if locked[a] {
            p[a] = original[a];
        }
    }
    if let Some(level) = level {
        for _ in 0..8 {
            let value = geometry.value(p) - level;
            if value.abs() < 1e-13 {
                break;
            }
            let mut g = geometry.gradient(p);
            for a in 0..3 {
                if locked[a] {
                    g[a] = 0.;
                }
            }
            let g2 = dot(g, g);
            if g2 < 1e-16 {
                return None;
            }
            for a in 0..3 {
                p[a] -= value * g[a] / g2;
            }
        }
        if (geometry.value(p) - level).abs() > 1e-11 {
            return None;
        }
    }
    geometry.contains(p).then_some(p)
}

#[allow(clippy::too_many_arguments)]
pub fn optimize(
    mesh: &mut Mesh,
    threshold: f64,
    h: f64,
    passes: usize,
    surface_weight: f64,
    threads: usize,
    accelerator: Option<&dyn crate::accelerator::Factory>,
) -> Result<(), String> {
    optimize_for(
        mesh,
        threshold,
        h,
        passes,
        surface_weight,
        threads,
        accelerator,
        &crate::geometry::Legacy,
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn optimize_for(
    mesh: &mut Mesh,
    threshold: f64,
    h: f64,
    passes: usize,
    surface_weight: f64,
    threads: usize,
    accelerator: Option<&dyn crate::accelerator::Factory>,
    geometry: &impl crate::geometry::Geometry,
) -> Result<(), String> {
    if passes == 0 {
        return Ok(());
    }
    let mut cpu_profile = crate::profile::Profile::new("optimizer");
    cpu_profile.count("passes", passes);
    cpu_profile.count("threads", threads);
    let profiling = accelerator.is_some() && std::env::var_os("MESHER_GPU_PROFILE").is_some();
    let stage_start = profiling.then(std::time::Instant::now);
    let mut proposal_seconds = 0.;
    let mut acceptance_seconds = 0.;
    let mut apply_seconds = 0.;
    let mut groups: BTreeMap<[i64; 3], Vec<usize>> = BTreeMap::new();
    for (i, &p) in mesh.points.iter().enumerate() {
        groups.entry(geometry.rank(p)).or_default().push(i);
    }
    cpu_profile.mark("periodic_groups");
    let mut incident = vec![Vec::new(); mesh.points.len()];
    for (i, t) in mesh.tets.iter().enumerate() {
        for &v in t {
            incident[v].push(i);
        }
    }
    cpu_profile.mark("tetrahedron_incidence");
    let mut surface_incident = vec![Vec::new(); mesh.points.len()];
    if surface_weight > 0. {
        for f in &mesh.surface {
            if geometry.cap(f.map(|i| mesh.points[i])) {
                continue;
            }
            for &v in f {
                surface_incident[v].push(*f);
            }
        }
    }
    cpu_profile.mark("surface_incidence");
    let groups: Vec<_> = groups
        .into_values()
        .map(|members| {
            let mut cells: Vec<_> = members
                .iter()
                .flat_map(|&i| incident[i].iter().copied())
                .collect();
            cells.sort_unstable();
            cells.dedup();
            let mut faces: Vec<_> = members
                .iter()
                .flat_map(|&i| surface_incident[i].iter().copied())
                .collect();
            faces.sort_unstable();
            faces.dedup();
            (members, cells, faces)
        })
        .collect();
    cpu_profile.mark("group_adjacency");
    let apply = |mesh: &mut Mesh, group: &Group, delta: Point| {
        for &i in &group.0 {
            for (a, &d) in delta.iter().enumerate() {
                mesh.points[i][a] += d;
            }
        }
    };
    if threads == 0 && accelerator.is_none() {
        for pass in 0..passes {
            for (index, group) in groups.iter().enumerate() {
                if index % 128 == 0 {
                    geometry.check()?;
                }
                let delta = if geometry.batched() {
                    batch::proposals(
                        mesh,
                        &[group],
                        &incident,
                        threshold,
                        h,
                        pass,
                        surface_weight,
                        geometry,
                    )?[0]
                } else {
                    proposal(
                        mesh,
                        group,
                        &incident,
                        threshold,
                        h,
                        pass,
                        surface_weight,
                        geometry,
                    )
                };
                apply(mesh, group, delta);
            }
        }
        cpu_profile.mark("serial_iterations");
        return Ok(());
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .map_err(|e| format!("CPU worker pool: {e}"))?;
    cpu_profile.mark("thread_pool");
    let mut node_group = vec![0; mesh.points.len()];
    for (gi, g) in groups.iter().enumerate() {
        for &v in &g.0 {
            node_group[v] = gi;
        }
    }
    let mut assigned = vec![usize::MAX; groups.len()];
    let mut colors: Vec<Vec<usize>> = Vec::new();
    for (gi, g) in groups.iter().enumerate() {
        let mut used = vec![false; colors.len()];
        for &ci in &g.1 {
            for &v in &mesh.tets[ci] {
                let color = assigned[node_group[v]];
                if color < used.len() {
                    used[color] = true;
                }
            }
        }
        let color = used.iter().position(|&v| !v).unwrap_or(colors.len());
        if color == colors.len() {
            colors.push(Vec::new());
        }
        colors[color].push(gi);
        assigned[gi] = color;
    }
    cpu_profile.mark("coloring");
    cpu_profile.count("colors", colors.len());
    let topology_seconds = stage_start.map_or(0., |start| start.elapsed().as_secs_f64());
    let mut gpu = accelerator
        .map(|factory| factory.create(mesh, &groups, &incident, &node_group))
        .transpose()?;
    cpu_profile.mark("backend_setup");
    for pass in 0..passes {
        for color in &colors {
            geometry.check()?;
            let updates: Vec<_> = if let Some(engine) = gpu.as_mut() {
                let start = std::time::Instant::now();
                let device_acceptance = engine.device_acceptance();
                let candidates =
                    engine.proposals(mesh, color, threshold, h, surface_weight, pass)?;
                if candidates.len() != color.len()
                    || candidates.iter().flatten().any(|v| !v.is_finite())
                {
                    return Err(
                        "accelerator returned invalid candidate dimensions or values".into(),
                    );
                }
                proposal_seconds += start.elapsed().as_secs_f64();
                let start = std::time::Instant::now();
                let accepted = if device_acceptance {
                    candidates
                } else {
                    pool.install(|| {
                        color
                            .par_iter()
                            .zip(candidates.par_iter())
                            .map(|(&gi, &delta)| {
                                accept_gpu(
                                    mesh,
                                    &groups[gi],
                                    threshold,
                                    h,
                                    surface_weight,
                                    delta,
                                    geometry,
                                )
                            })
                            .collect()
                    })
                };
                acceptance_seconds += start.elapsed().as_secs_f64();
                accepted
            } else if geometry.batched() {
                let mut updates = Vec::with_capacity(color.len());
                for chunk in color.chunks(128) {
                    let selected: Vec<_> = chunk.iter().map(|&i| &groups[i]).collect();
                    updates.extend(batch::proposals(
                        mesh,
                        &selected,
                        &incident,
                        threshold,
                        h,
                        pass,
                        surface_weight,
                        geometry,
                    )?);
                }
                updates
            } else {
                pool.install(|| {
                    color
                        .par_iter()
                        .map(|&gi| {
                            proposal(
                                mesh,
                                &groups[gi],
                                &incident,
                                threshold,
                                h,
                                pass,
                                surface_weight,
                                geometry,
                            )
                        })
                        .collect()
                })
            };
            let apply_start = std::time::Instant::now();
            let mut changed_nodes = Vec::new();
            for (&gi, delta) in color.iter().zip(updates) {
                apply(mesh, &groups[gi], delta);
                if gpu.is_some() && delta != [0.; 3] {
                    changed_nodes.extend_from_slice(&groups[gi].0);
                }
            }
            if let Some(engine) = gpu.as_mut() {
                engine.synchronize_updates(mesh, &changed_nodes)?;
            }
            apply_seconds += apply_start.elapsed().as_secs_f64();
        }
    }
    cpu_profile.mark("colored_iterations");
    if profiling {
        eprintln!(
            "RUST_GPU_PROFILE {}",
            serde_json::json!({"topology_seconds":topology_seconds,"proposal_seconds":proposal_seconds,"acceptance_seconds":acceptance_seconds,"apply_seconds":apply_seconds,"stage_seconds":stage_start.map_or(0., |start| start.elapsed().as_secs_f64())})
        );
    }
    Ok(())
}
fn accept_gpu(
    mesh: &Mesh,
    group: &Group,
    threshold: f64,
    h: f64,
    weight: f64,
    delta: Point,
    geometry: &impl crate::geometry::Geometry,
) -> Point {
    if !delta.iter().all(|v| v.is_finite()) || delta == [0.; 3] {
        return [0.; 3];
    }
    let origin = mesh.points[group.0[0]];
    let value = geometry.value(origin);
    let level = (if geometry.band() {
        (value.abs() - threshold).abs()
    } else {
        value.abs()
    } < 1e-8)
        .then_some(geometry.level(value, threshold));
    let Some(candidate) = project(
        std::array::from_fn(|a| origin[a] + delta[a]),
        origin,
        level,
        geometry,
    ) else {
        return [0.; 3];
    };
    let (old_min, old_cost) = evaluate(mesh, group, threshold, h, weight, origin, geometry);
    let (new_min, new_cost) = evaluate(mesh, group, threshold, h, weight, candidate, geometry);
    if new_min >= old_min.min(0.65_f64.powi(3)) - 1e-12 && new_cost < old_cost {
        sub(candidate, origin)
    } else {
        [0.; 3]
    }
}

#[allow(clippy::too_many_arguments)]
fn proposal(
    mesh: &Mesh,
    group: &Group,
    incident: &[Vec<usize>],
    threshold: f64,
    h: f64,
    pass: usize,
    surface_weight: f64,
    geometry: &impl crate::geometry::Geometry,
) -> Point {
    let (members, _, faces) = group;
    let original = mesh.points[members[0]];
    let value = geometry.value(original);
    let level = (if geometry.band() {
        (value.abs() - threshold).abs()
    } else {
        value.abs()
    } < 1e-8)
        .then_some(geometry.level(value, threshold));
    let evaluate = |candidate| {
        evaluate(
            mesh,
            group,
            threshold,
            h,
            surface_weight,
            candidate,
            geometry,
        )
    };
    let (old_min, old_cost) = evaluate(original);
    if faces.is_empty() && old_min > 0.8_f64.powi(3) {
        return [0.; 3];
    }
    let mut best = original;
    let mut best_cost = old_cost;
    let floor = old_min.min(0.65_f64.powi(3));
    let mut try_candidate = |candidate: Point| {
        if let Some(p) = project(candidate, original, level, geometry) {
            let (worst, cost) = evaluate_bounded(
                mesh,
                group,
                threshold,
                h,
                surface_weight,
                p,
                geometry,
                floor,
                best_cost,
            );
            if worst >= floor - 1e-12 && cost < best_cost {
                best = p;
                best_cost = cost;
            }
        }
    };
    // A projected Laplacian proposal often fixes the first large defects.
    let mut avg = [0.; 3];
    let mut count = 0;
    for &member in members {
        let shift = sub(mesh.points[member], original);
        for &ci in &incident[member] {
            for &j in &mesh.tets[ci] {
                if j != member {
                    for a in 0..3 {
                        avg[a] += mesh.points[j][a] - shift[a];
                    }
                    count += 1;
                }
            }
        }
    }
    if count > 0 {
        let avg = avg.map(|x| x / count as f64);
        for weight in [0.5, 0.25] {
            try_candidate(std::array::from_fn(|a| {
                original[a] + weight * (avg[a] - original[a])
            }));
        }
    }
    // Finite differences of a smooth penalty focus work on poor tets.
    let eps = h * 0.002;
    let mut direction = [0.; 3];
    for a in 0..3 {
        let mut lo = original;
        lo[a] -= eps;
        let mut hi = original;
        hi[a] += eps;
        let c0 = project(lo, original, level, geometry).map_or(old_cost, |p| evaluate(p).1);
        let c1 = project(hi, original, level, geometry).map_or(old_cost, |p| evaluate(p).1);
        if c0.is_finite() && c1.is_finite() {
            direction[a] = c0 - c1;
        }
    }
    let length = norm(direction);
    if length > 1e-15 {
        let scale = if pass < 3 { 0.2 } else { 0.1 };
        for step in [scale, scale / 2., scale / 4., scale / 8.] {
            try_candidate(std::array::from_fn(|a| {
                original[a] + h * step * direction[a] / length
            }));
        }
    }
    sub(best, original)
}

fn evaluate(
    mesh: &Mesh,
    group: &Group,
    threshold: f64,
    h: f64,
    surface_weight: f64,
    candidate: Point,
    geometry: &impl crate::geometry::Geometry,
) -> (f64, f64) {
    evaluate_bounded(
        mesh,
        group,
        threshold,
        h,
        surface_weight,
        candidate,
        geometry,
        0.,
        f64::INFINITY,
    )
}

// A candidate below this floor will be rejected regardless of its cost. Avoid
// evaluating its remaining elements and implicit surface samples. Cost-gradient
// probes still use the unbounded evaluator above.
#[allow(clippy::too_many_arguments)]
fn evaluate_bounded(
    mesh: &Mesh,
    group: &Group,
    threshold: f64,
    h: f64,
    surface_weight: f64,
    candidate: Point,
    geometry: &impl crate::geometry::Geometry,
    floor: f64,
    cost_limit: f64,
) -> (f64, f64) {
    let (members, cells, faces) = group;
    let original = mesh.points[members[0]];
    let delta = sub(candidate, original);
    let mut worst = 1_f64;
    let mut cost = 0.;
    for &ci in cells {
        let p = mesh.tets[ci].map(|i| {
            if members.contains(&i) {
                std::array::from_fn(|a| mesh.points[i][a] + delta[a])
            } else {
                mesh.points[i]
            }
        });
        if determinant(p) <= h * h * h * 1e-8 {
            return (0., f64::INFINITY);
        }
        let l2: f64 = (0..4)
            .flat_map(|i| (i + 1..4).map(move |j| dot(sub(p[i], p[j]), sub(p[i], p[j]))))
            .sum();
        // Squared MMG quality avoids fractional powers in the hot loop.
        let q = 432. * determinant(p).powi(2) / l2.powi(3);
        if q < floor - 1e-12 {
            return (q, f64::INFINITY);
        }
        worst = worst.min(q);
        cost += 1. / (q * q);
        if surface_weight >= 0. && cost >= cost_limit {
            return (worst, f64::INFINITY);
        }
    }
    for f in faces {
        let p = f.map(|i| {
            if members.contains(&i) {
                std::array::from_fn(|a| mesh.points[i][a] + delta[a])
            } else {
                mesh.points[i]
            }
        });
        let center: Point = std::array::from_fn(|a| (p[0][a] + p[1][a] + p[2][a]) / 3.);
        let denom = norm(geometry.gradient(center)) * h * h * 0.3;
        for sample in [
            center,
            std::array::from_fn(|a| (p[0][a] + p[1][a]) * 0.5),
            std::array::from_fn(|a| (p[1][a] + p[2][a]) * 0.5),
            std::array::from_fn(|a| (p[2][a] + p[0][a]) * 0.5),
        ] {
            let error = geometry.residual(sample, threshold) / denom;
            cost += surface_weight * error.powi(2);
            if surface_weight >= 0. && cost >= cost_limit {
                return (worst, f64::INFINITY);
            }
        }
    }
    (worst, cost)
}

#[cfg(test)]
mod gpu_acceptance_tests {
    use super::*;
    #[test]
    fn bounded_cost_preserves_accepted_candidates_and_skips_rejected_surfaces() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        struct Field(AtomicUsize);
        impl crate::geometry::Geometry for Field {
            fn value(&self, p: Point) -> f64 {
                p[0]
            }
            fn gradient(&self, _: Point) -> Point {
                self.0.fetch_add(1, Ordering::Relaxed);
                [1., 0., 0.]
            }
            fn rank(&self, _: Point) -> [i64; 3] {
                [0; 3]
            }
            fn bounds(&self) -> [Point; 2] {
                [[0.; 3], [1.; 3]]
            }
        }
        let field = Field(AtomicUsize::new(0));
        let mesh = Mesh {
            points: vec![
                [0.2, 0.2, 0.2],
                [0.8, 0.2, 0.2],
                [0.2, 0.8, 0.2],
                [0.2, 0.2, 0.8],
            ],
            tets: vec![[0, 1, 2, 3]],
            surface: vec![[0, 1, 2]],
        };
        let group = (vec![0], vec![0], mesh.surface.clone());
        let original = mesh.points[0];
        let full = evaluate(&mesh, &group, 0.5, 0.6, 3., original, &field);
        assert_eq!(
            evaluate_bounded(
                &mesh,
                &group,
                0.5,
                0.6,
                3.,
                original,
                &field,
                full.0,
                f64::INFINITY
            ),
            full
        );
        field.0.store(0, Ordering::Relaxed);
        let bad = evaluate_bounded(
            &mesh,
            &group,
            0.5,
            0.6,
            3.,
            [0.39; 3],
            &field,
            full.0,
            f64::INFINITY,
        );
        assert!(bad.0 < full.0 && bad.1.is_infinite());
        assert_eq!(field.0.load(Ordering::Relaxed), 0);
        let too_costly = evaluate_bounded(&mesh, &group, 0.5, 0.6, 3., original, &field, 0., 0.);
        assert!(too_costly.1.is_infinite());
        assert_eq!(field.0.load(Ordering::Relaxed), 0);
    }
    #[test]
    fn rejects_nonfinite_outside_and_inverting_gpu_moves() {
        let mesh = Mesh {
            points: vec![
                [0.2, 0.2, 0.2],
                [0.8, 0.2, 0.2],
                [0.2, 0.8, 0.2],
                [0.2, 0.2, 0.8],
            ],
            tets: vec![[0, 1, 2, 3]],
            surface: vec![],
        };
        let group = (vec![0], vec![0], vec![]);
        assert!(determinant(mesh.tets[0].map(|i| mesh.points[i])) > 0.);
        for delta in [[f64::NAN, 0., 0.], [-2., 0., 0.], [0.6, 0.6, 0.6]] {
            assert_eq!(
                accept_gpu(&mesh, &group, 0.5, 0.6, 0., delta, &crate::geometry::Legacy),
                [0.; 3]
            );
        }
    }
}
