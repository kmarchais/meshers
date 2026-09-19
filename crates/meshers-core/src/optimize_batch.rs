//! Vector evaluation of independent proposals, preserving scalar candidate order.
use crate::{Mesh, Point, accelerator::Group, determinant, dot, geometry::Geometry, norm, sub};

fn project(
    points: &[Point],
    origins: &[Point],
    levels: &[Option<f64>],
    g: &impl Geometry,
) -> Result<Vec<Option<Point>>, String> {
    let mut result: Vec<_> = points
        .iter()
        .zip(origins)
        .map(|(&p, &o)| {
            Some(std::array::from_fn(|a| {
                if g.locked(o)[a] { o[a] } else { p[a] }
            }))
        })
        .collect();
    let mut active: Vec<_> = (0..points.len()).filter(|&i| levels[i].is_some()).collect();
    for _ in 0..8 {
        if active.is_empty() {
            break;
        }
        g.check()?;
        let samples: Vec<_> = active.iter().map(|&i| result[i].unwrap()).collect();
        let values = g.values(&samples);
        let remaining: Vec<_> = active
            .iter()
            .zip(&values)
            .filter_map(|(&i, &v)| {
                ((v - levels[i].unwrap()).abs() >= 1e-13).then_some((i, v - levels[i].unwrap()))
            })
            .collect();
        let samples: Vec<_> = remaining.iter().map(|&(i, _)| result[i].unwrap()).collect();
        let gradients = g.gradients(&samples);
        active.clear();
        for ((i, v), mut grad) in remaining.into_iter().zip(gradients) {
            let locked = g.locked(origins[i]);
            for a in 0..3 {
                if locked[a] {
                    grad[a] = 0.;
                }
            }
            let g2 = dot(grad, grad);
            if g2 < 1e-16 {
                result[i] = None;
                continue;
            }
            let p = result[i].unwrap();
            result[i] = Some(std::array::from_fn(|a| p[a] - v * grad[a] / g2));
            active.push(i);
        }
    }
    let ids: Vec<_> = (0..points.len())
        .filter(|&i| levels[i].is_some() && result[i].is_some())
        .collect();
    let samples: Vec<_> = ids.iter().map(|&i| result[i].unwrap()).collect();
    for (i, v) in ids.into_iter().zip(g.values(&samples)) {
        if (v - levels[i].unwrap()).abs() > 1e-11 {
            result[i] = None;
        }
    }
    for p in &mut result {
        if p.is_some_and(|v| !g.contains(v)) {
            *p = None;
        }
    }
    g.check()?;
    Ok(result)
}

fn evaluate(
    mesh: &Mesh,
    queries: &[(&Group, Point)],
    threshold: f64,
    h: f64,
    weight: f64,
    g: &impl Geometry,
) -> Result<Vec<(f64, f64)>, String> {
    let mut results = Vec::with_capacity(queries.len());
    // Bounded by query count, not by the total mesh size.
    for chunk in queries.chunks(128) {
        g.check()?;
        let mut current = Vec::with_capacity(chunk.len());
        let mut centers = Vec::new();
        let mut samples = Vec::new();
        let mut owners = Vec::new();
        for (qi, &(group, candidate)) in chunk.iter().enumerate() {
            let (members, cells, faces) = group;
            let delta = sub(candidate, mesh.points[members[0]]);
            let move_point = |i| {
                if members.contains(&i) {
                    std::array::from_fn(|a| mesh.points[i][a] + delta[a])
                } else {
                    mesh.points[i]
                }
            };
            let mut worst = 1_f64;
            let mut cost = 0.;
            for &ci in cells {
                let p = mesh.tets[ci].map(move_point);
                if determinant(p) <= h * h * h * 1e-8 {
                    worst = 0.;
                    cost = f64::INFINITY;
                    break;
                }
                let l2: f64 = (0..4)
                    .flat_map(|i| (i + 1..4).map(move |j| dot(sub(p[i], p[j]), sub(p[i], p[j]))))
                    .sum();
                let q = 432. * determinant(p).powi(2) / l2.powi(3);
                worst = worst.min(q);
                cost += 1. / (q * q);
            }
            current.push((worst, cost));
            if !cost.is_finite() {
                continue;
            }
            for f in faces {
                let p = f.map(move_point);
                let c = std::array::from_fn(|a| (p[0][a] + p[1][a] + p[2][a]) / 3.);
                centers.push(c);
                owners.push(qi);
                samples.extend([
                    c,
                    std::array::from_fn(|a| (p[0][a] + p[1][a]) * 0.5),
                    std::array::from_fn(|a| (p[1][a] + p[2][a]) * 0.5),
                    std::array::from_fn(|a| (p[2][a] + p[0][a]) * 0.5),
                ]);
            }
        }
        let gradients = g.gradients(&centers);
        let values = g.values(&samples);
        for ((owner, grad), values) in owners
            .into_iter()
            .zip(gradients)
            .zip(values.as_chunks::<4>().0)
        {
            let denom = norm(grad) * h * h * 0.3;
            for &value in values {
                let residual = if g.band() {
                    value.abs() - threshold
                } else {
                    value
                };
                current[owner].1 += weight * (residual / denom).powi(2);
            }
        }
        results.extend(current);
    }
    g.check()?;
    Ok(results)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn proposals(
    mesh: &Mesh,
    groups: &[&Group],
    incident: &[Vec<usize>],
    threshold: f64,
    h: f64,
    pass: usize,
    weight: f64,
    g: &impl Geometry,
) -> Result<Vec<Point>, String> {
    let origins: Vec<_> = groups.iter().map(|g| mesh.points[g.0[0]]).collect();
    let values = g.values(&origins);
    let levels: Vec<_> = values
        .into_iter()
        .map(|v| {
            (if g.band() {
                (v.abs() - threshold).abs()
            } else {
                v.abs()
            } < 1e-8)
                .then_some(g.level(v, threshold))
        })
        .collect();
    let old = evaluate(
        mesh,
        &groups
            .iter()
            .zip(&origins)
            .map(|(&g, &p)| (g, p))
            .collect::<Vec<_>>(),
        threshold,
        h,
        weight,
        g,
    )?;
    let mut best = origins.clone();
    let mut best_cost: Vec<_> = old.iter().map(|v| v.1).collect();
    let mut candidates = Vec::new();
    let mut owners = Vec::new();
    let mut slots = Vec::new();
    for (i, group) in groups.iter().enumerate() {
        if group.2.is_empty() && old[i].0 > 0.8_f64.powi(3) {
            continue;
        }
        let original = origins[i];
        let mut avg = [0.; 3];
        let mut count = 0;
        for &member in &group.0 {
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
            for w in [0.5, 0.25] {
                candidates.push(std::array::from_fn(|a| {
                    original[a] + w * (avg[a] - original[a])
                }));
                owners.push(i);
                slots.push(None);
            }
        }
        for a in 0..3 {
            for sign in [-1., 1.] {
                let mut p = original;
                p[a] += sign * (h * 0.002);
                candidates.push(p);
                owners.push(i);
                slots.push(Some(2 * a + usize::from(sign > 0.)));
            }
        }
    }
    let projected = project(
        &candidates,
        &owners.iter().map(|&i| origins[i]).collect::<Vec<_>>(),
        &owners.iter().map(|&i| levels[i]).collect::<Vec<_>>(),
        g,
    )?;
    let valid: Vec<_> = projected
        .iter()
        .enumerate()
        .filter_map(|(j, p)| p.map(|p| (j, p)))
        .collect();
    let scores = evaluate(
        mesh,
        &valid
            .iter()
            .map(|&(j, p)| (groups[owners[j]], p))
            .collect::<Vec<_>>(),
        threshold,
        h,
        weight,
        g,
    )?;
    let mut differences: Vec<_> = old.iter().map(|v| [v.1; 6]).collect();
    for ((j, p), (worst, cost)) in valid.into_iter().zip(scores) {
        let i = owners[j];
        if let Some(slot) = slots[j] {
            differences[i][slot] = cost;
        } else if worst >= old[i].0.min(0.65_f64.powi(3)) - 1e-12 && cost < best_cost[i] {
            best[i] = p;
            best_cost[i] = cost;
        }
    }
    candidates.clear();
    owners.clear();
    for i in 0..groups.len() {
        let d: Point = std::array::from_fn(|a| {
            let (lo, hi) = (differences[i][2 * a], differences[i][2 * a + 1]);
            if lo.is_finite() && hi.is_finite() {
                lo - hi
            } else {
                0.
            }
        });
        let length = norm(d);
        if length > 1e-15 {
            let scale = if pass < 3 { 0.2 } else { 0.1 };
            for step in [scale, scale / 2., scale / 4., scale / 8.] {
                candidates.push(std::array::from_fn(|a| {
                    origins[i][a] + h * step * d[a] / length
                }));
                owners.push(i);
            }
        }
    }
    let projected = project(
        &candidates,
        &owners.iter().map(|&i| origins[i]).collect::<Vec<_>>(),
        &owners.iter().map(|&i| levels[i]).collect::<Vec<_>>(),
        g,
    )?;
    let valid: Vec<_> = projected
        .iter()
        .enumerate()
        .filter_map(|(j, p)| p.map(|p| (j, p)))
        .collect();
    let scores = evaluate(
        mesh,
        &valid
            .iter()
            .map(|&(j, p)| (groups[owners[j]], p))
            .collect::<Vec<_>>(),
        threshold,
        h,
        weight,
        g,
    )?;
    for ((j, p), (worst, cost)) in valid.into_iter().zip(scores) {
        let i = owners[j];
        if worst >= old[i].0.min(0.65_f64.powi(3)) - 1e-12 && cost < best_cost[i] {
            best[i] = p;
            best_cost[i] = cost;
        }
    }
    Ok(best
        .into_iter()
        .zip(origins)
        .map(|(p, o)| sub(p, o))
        .collect())
}
