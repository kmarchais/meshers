//! Conservative short-edge removal, preserving constraint strata and periodic orbits.
use super::{MAX_CONSTRAINTS, Output};
use crate::implicit::Options;
use crate::{Point, Tet, determinant, quality};
use std::collections::{BTreeMap, BTreeSet};

type Simplex = Vec<usize>;
#[cfg(feature = "wide-constraints")]
type StratumMask = u64;
#[cfg(not(feature = "wide-constraints"))]
type StratumMask = u16;
fn subsets(vertices: &[usize], out: &mut BTreeSet<Simplex>) {
    for bits in 1..(1 << vertices.len()) {
        let mut s: Vec<_> = vertices
            .iter()
            .enumerate()
            .filter_map(|(j, &v)| (bits & (1 << j) != 0).then_some(v))
            .collect();
        s.sort_unstable();
        out.insert(s);
    }
}
fn link_ok(a: usize, b: usize, ids: &BTreeSet<usize>, tets: &[Tet], masks: &[StratumMask]) -> bool {
    let mut la = BTreeSet::new();
    let mut lb = BTreeSet::new();
    let mut le = BTreeSet::new();
    let mut faces = BTreeMap::<[usize; 3], usize>::new();
    for &i in ids {
        let t = tets[i];
        if t.contains(&a) {
            subsets(
                &t.into_iter().filter(|&v| v != a).collect::<Vec<_>>(),
                &mut la,
            );
        }
        if t.contains(&b) {
            subsets(
                &t.into_iter().filter(|&v| v != b).collect::<Vec<_>>(),
                &mut lb,
            );
        }
        if t.contains(&a) && t.contains(&b) {
            subsets(
                &t.into_iter()
                    .filter(|&v| v != a && v != b)
                    .collect::<Vec<_>>(),
                &mut le,
            );
        }
        for mut f in [
            [t[0], t[1], t[2]],
            [t[0], t[1], t[3]],
            [t[0], t[2], t[3]],
            [t[1], t[2], t[3]],
        ] {
            if f.contains(&a) || f.contains(&b) {
                f.sort_unstable();
                *faces.entry(f).or_default() += 1;
            }
        }
    }
    if la.intersection(&lb).cloned().collect::<BTreeSet<_>>() != le {
        return false;
    }
    la.clear();
    lb.clear();
    le.clear();
    let boundary: Vec<_> = faces
        .iter()
        .filter_map(|(&f, &count)| (count == 1).then_some(f))
        .collect();
    for (f, count) in faces {
        if count == 1 {
            if f.contains(&a) {
                subsets(
                    &f.into_iter().filter(|&v| v != a).collect::<Vec<_>>(),
                    &mut la,
                );
            }
            if f.contains(&b) {
                subsets(
                    &f.into_iter().filter(|&v| v != b).collect::<Vec<_>>(),
                    &mut lb,
                );
            }
            if f.contains(&a) && f.contains(&b) {
                subsets(
                    &f.into_iter()
                        .filter(|&v| v != a && v != b)
                        .collect::<Vec<_>>(),
                    &mut le,
                );
            }
        }
    }
    if la.intersection(&lb).cloned().collect::<BTreeSet<_>>() != le {
        return false;
    }
    feature_link_ok(a, b, &boundary, masks)
}

// A topologically valid solid collapse can still connect two separate rims.
// Preserve the feature subcomplex too: no new labelled curve edge, no lost
// curve edge except the collapsed one, and the 1D link condition on each rim.
fn feature_link_ok(a: usize, b: usize, faces: &[[usize; 3]], masks: &[StratumMask]) -> bool {
    let features = |faces: &[[usize; 3]]| {
        let mut edges = BTreeSet::new();
        for f in faces {
            for [u, v] in [[f[0], f[1]], [f[1], f[2]], [f[2], f[0]]] {
                let common = masks[u] & masks[v];
                if common.count_ones() < 2 {
                    continue;
                }
                for x in 0..MAX_CONSTRAINTS + 6 {
                    if common & (1 << x) == 0 {
                        continue;
                    }
                    for y in x + 1..MAX_CONSTRAINTS + 6 {
                        let bits = (1 << x) | (1 << y);
                        if common & bits == bits {
                            edges.insert((bits, u.min(v), u.max(v)));
                        }
                    }
                }
            }
        }
        edges
    };
    let old = features(faces);
    let mapped: Vec<_> = faces
        .iter()
        .filter_map(|f| {
            let f = f.map(|v| if v == a { b } else { v });
            (f[0] != f[1] && f[1] != f[2] && f[0] != f[2]).then_some(f)
        })
        .collect();
    let new = features(&mapped);
    let expected: BTreeSet<_> = old
        .iter()
        .filter_map(|&(bits, u, v)| {
            let u = if u == a { b } else { u };
            let v = if v == a { b } else { v };
            (u != v).then_some((bits, u.min(v), u.max(v)))
        })
        .collect();
    if new != expected {
        return false;
    }
    for &(bits, u, v) in &old {
        if u == a.min(b) && v == a.max(b) {
            let neighbors = |node| {
                old.iter()
                    .filter_map(|&(m, x, y)| {
                        if m != bits {
                            None
                        } else if x == node {
                            Some(y)
                        } else if y == node {
                            Some(x)
                        } else {
                            None
                        }
                    })
                    .collect::<BTreeSet<_>>()
            };
            if neighbors(a).intersection(&neighbors(b)).next().is_some() {
                return false;
            }
        }
    }
    true
}
fn q(points: &[Point], t: Tet) -> f64 {
    quality(t.map(|i| points[i])).powf(1.5)
}
fn cap(p: Point, o: &Options) -> StratumMask {
    (0..3).fold(0, |bits, a| {
        bits | if p[a] == o.bounds[0][a] {
            1 << (MAX_CONSTRAINTS + 2 * a)
        } else if p[a] == o.bounds[1][a] {
            1 << (MAX_CONSTRAINTS + 1 + 2 * a)
        } else {
            0
        }
    })
}

pub(super) fn improve(output: &mut Output, o: &Options) -> Result<usize, String> {
    let mesh = &mut output.mesh;
    let n = mesh.points.len();
    let mut incident = vec![Vec::new(); n];
    for (i, t) in mesh.tets.iter().enumerate() {
        for &v in t {
            incident[v].push(i);
        }
    }
    let mut alive = vec![true; mesh.tets.len()];
    let mut active = vec![true; n];
    let rank = |p: Point| -> [i64; 3] {
        std::array::from_fn(|a| {
            let [lo, hi] = o.bounds;
            let x = if o.periodic[a] && p[a] == hi[a] {
                lo[a]
            } else {
                p[a]
            };
            ((x - lo[a]) / (hi[a] - lo[a]) * 1e11).round() as i64
        })
    };
    let mut groups = BTreeMap::<[i64; 3], Vec<usize>>::new();
    for (i, &p) in output.parameters.iter().enumerate() {
        groups.entry(rank(p)).or_default().push(i);
    }
    let h = (0..3)
        .map(|a| (o.bounds[1][a] - o.bounds[0][a]) / o.cells[a] as f64)
        .fold(f64::INFINITY, f64::min);
    let masks: Vec<_> = output
        .constraints
        .iter()
        .zip(&output.parameters)
        .map(|(&m, &p)| StratumMask::from(m) | cap(p, o))
        .collect();
    let mut accepted = 0;
    for _ in 0..o.optimize_passes {
        let mut edges = BTreeSet::new();
        for (i, &t) in mesh.tets.iter().enumerate() {
            if alive[i] && q(&mesh.points, t) < 0.1 {
                for a in 0..4 {
                    for b in a + 1..4 {
                        edges.insert([t[a].min(t[b]), t[a].max(t[b])]);
                    }
                }
            }
        }
        let length = |[a, b]: [usize; 2]| crate::norm(crate::sub(mesh.points[a], mesh.points[b]));
        let mut edges: Vec<_> = edges
            .into_iter()
            .filter(|&e| length(e) < 0.75 * h)
            .collect();
        edges.sort_by(|&a, &b| length(a).total_cmp(&length(b)).then(a.cmp(&b)));
        let before = accepted;
        for [x, y] in edges {
            if !active[x] || !active[y] {
                continue;
            }
            for (a, b) in [(x, y), (y, x)] {
                if masks[a] & !masks[b] != 0 || masks[a].count_ones() >= 3 {
                    continue;
                }
                let ga = &groups[&rank(output.parameters[a])];
                let gb = &groups[&rank(output.parameters[b])];
                let aa: Vec<_> = ga.iter().copied().filter(|&v| active[v]).collect();
                let bb: Vec<_> = gb.iter().copied().filter(|&v| active[v]).collect();
                if aa.len() != bb.len() {
                    continue;
                }
                let mut replace = BTreeMap::new();
                for &v in &aa {
                    if let Some(&w) = bb
                        .iter()
                        .find(|&&w| cap(output.parameters[v], o) == cap(output.parameters[w], o))
                    {
                        replace.insert(v, w);
                    }
                }
                if replace.len() != aa.len()
                    || replace.values().copied().collect::<BTreeSet<_>>().len() != bb.len()
                {
                    continue;
                }
                let mut all = BTreeSet::new();
                let mut valid = true;
                for (&v, &w) in &replace {
                    if masks[v] & !masks[w] != 0 || masks[v].count_ones() >= 3 {
                        valid = false;
                        break;
                    }
                    let ids: BTreeSet<_> = incident[v]
                        .iter()
                        .chain(&incident[w])
                        .copied()
                        .filter(|&i| {
                            alive[i] && (mesh.tets[i].contains(&v) || mesh.tets[i].contains(&w))
                        })
                        .collect();
                    if !ids
                        .iter()
                        .any(|&i| mesh.tets[i].contains(&v) && mesh.tets[i].contains(&w))
                        || !link_ok(v, w, &ids, &mesh.tets, &masks)
                    {
                        valid = false;
                        break;
                    }
                    all.extend(ids);
                }
                if !valid {
                    continue;
                }
                let old_min = all
                    .iter()
                    .map(|&i| q(&mesh.points, mesh.tets[i]))
                    .fold(1., f64::min);
                let old_bad = all
                    .iter()
                    .filter(|&&i| q(&mesh.points, mesh.tets[i]) < 0.1)
                    .count();
                let mut new_min: f64 = 1.;
                let mut new_bad = 0;
                let mut updates = Vec::new();
                for &i in &all {
                    let t = mesh.tets[i].map(|v| *replace.get(&v).unwrap_or(&v));
                    if t.into_iter().collect::<BTreeSet<_>>().len() < 4 {
                        updates.push((i, None));
                        continue;
                    }
                    let det = determinant(t.map(|v| mesh.points[v]));
                    let quality = q(&mesh.points, t);
                    if !det.is_finite() || det <= 0. || quality + 1e-12 < old_min.min(0.1) {
                        valid = false;
                        break;
                    }
                    new_min = new_min.min(quality);
                    new_bad += usize::from(quality < 0.1);
                    updates.push((i, Some(t)));
                }
                if !valid || (new_bad >= old_bad && new_min <= old_min * 1.001) {
                    continue;
                }
                for (i, t) in updates {
                    if let Some(t) = t {
                        mesh.tets[i] = t;
                    } else {
                        alive[i] = false;
                    }
                }
                for (&v, &w) in &replace {
                    active[v] = false;
                    let extra = incident[v].clone();
                    incident[w].extend(extra);
                    incident[w].sort_unstable();
                    incident[w].dedup();
                }
                accepted += replace.len();
                break;
            }
        }
        if before == accepted {
            break;
        }
    }
    mesh.tets = mesh
        .tets
        .iter()
        .enumerate()
        .filter_map(|(i, &t)| alive[i].then_some(t))
        .collect();
    let mut remap = vec![usize::MAX; n];
    let mut points = Vec::new();
    let mut parameters = Vec::new();
    let mut constraints = Vec::new();
    for t in &mut mesh.tets {
        for v in t {
            if remap[*v] == usize::MAX {
                remap[*v] = points.len();
                points.push(mesh.points[*v]);
                parameters.push(output.parameters[*v]);
                constraints.push(output.constraints[*v]);
            }
            *v = remap[*v];
        }
    }
    mesh.points = points;
    output.parameters = parameters;
    output.constraints = constraints;
    mesh.surface = crate::boundary(&mesh.tets)?;
    if !super::super::nonmanifold_vertices(&mesh.surface, mesh.points.len()).is_empty() {
        return Err("quality repair produced nonmanifold boundary".into());
    }
    Ok(accepted)
}

pub(super) fn smooth(
    output: &mut Output,
    o: &Options,
    fields: &[&dyn super::ScalarField],
    map: &impl Fn(Point) -> Point,
    regions: Option<&[super::ConstraintMask]>,
) -> Result<usize, String> {
    let [lo, hi] = o.bounds;
    let h = (0..3)
        .map(|a| (hi[a] - lo[a]) / o.cells[a] as f64)
        .fold(f64::INFINITY, f64::min);
    let canonical = |p: Point| -> Point {
        std::array::from_fn(|a| {
            if o.periodic[a] && p[a] == hi[a] {
                lo[a]
            } else {
                p[a]
            }
        })
    };
    let rank = |p: Point| -> [i64; 3] {
        let p = canonical(p);
        std::array::from_fn(|a| ((p[a] - lo[a]) / (hi[a] - lo[a]) * 1e11).round() as i64)
    };
    let mut groups = BTreeMap::<[i64; 3], Vec<usize>>::new();
    let mut incident = vec![Vec::new(); output.mesh.points.len()];
    for (i, t) in output.mesh.tets.iter().enumerate() {
        for &v in t {
            incident[v].push(i);
        }
    }
    for (i, &p) in output.parameters.iter().enumerate() {
        groups.entry(rank(p)).or_default().push(i);
    }
    let mut moved = 0;
    for (iteration, group) in groups.values().enumerate() {
        if iteration % 4096 == 0 {
            for f in fields {
                f.check()?;
            }
        }
        let ids: BTreeSet<_> = group
            .iter()
            .flat_map(|&v| incident[v].iter().copied())
            .collect();
        let old_min = ids
            .iter()
            .map(|&i| q(&output.mesh.points, output.mesh.tets[i]))
            .fold(1., f64::min);
        if old_min >= 0.15 {
            continue;
        }
        let first = group[0];
        let original = output.parameters[first];
        let origin = canonical(original);
        let neighbors: BTreeSet<_> = incident[first]
            .iter()
            .flat_map(|&i| output.mesh.tets[i])
            .filter(|&v| v != first)
            .collect();
        let mean: Point = std::array::from_fn(|a| {
            neighbors
                .iter()
                .map(|&v| output.parameters[v][a])
                .sum::<f64>()
                / neighbors.len() as f64
        });
        let direction = crate::sub(mean, original);
        let mut candidates = Vec::new();
        for alpha in [1., 0.5, 0.25] {
            candidates.push(std::array::from_fn(|a| origin[a] + alpha * direction[a]));
        }
        for axis in 0..3 {
            for sign in [-1., 1.] {
                let mut p = origin;
                p[axis] += sign * 0.15 * h;
                candidates.push(p);
            }
        }
        let active: Vec<_> = (0..fields.len())
            .filter(|&j| output.constraints[first] & (1 << j) != 0)
            .collect();
        let mut best_score = old_min;
        let mut best = None;
        for candidate in candidates {
            let Some(p) = super::project_candidate(fields, &active, origin, candidate, o, h * 1e-5)
            else {
                continue;
            };
            if (0..3).any(|a| p[a] < lo[a] || p[a] > hi[a])
                || crate::norm(crate::sub(p, origin)) > 0.5 * h
            {
                continue;
            }
            let outside = |mask: super::ConstraintMask| {
                fields
                    .iter()
                    .enumerate()
                    .any(|(i, f)| mask & (1 << i) != 0 && f.value(p) > 1e-9)
            };
            if !active.is_empty()
                && regions.map_or_else(
                    || outside(super::ConstraintMask::MAX),
                    |regions| regions.iter().all(|&mask| outside(mask)),
                )
            {
                continue;
            }
            let parameters: Vec<Point> = group
                .iter()
                .map(|&v| {
                    std::array::from_fn(|a| {
                        if o.periodic[a] && output.parameters[v][a] == hi[a] {
                            hi[a]
                        } else {
                            p[a]
                        }
                    })
                })
                .collect();
            let physical: Vec<Point> = parameters.iter().map(|&p| map(p)).collect();
            if physical.iter().flatten().any(|x| !x.is_finite()) {
                continue;
            }
            let mut minimum: f64 = 1.;
            let mut valid = true;
            for &i in &ids {
                let tet = output.mesh.tets[i].map(|v| {
                    group
                        .iter()
                        .position(|&g| g == v)
                        .map_or(output.mesh.points[v], |j| physical[j])
                });
                if determinant(tet) <= 0. {
                    valid = false;
                    break;
                }
                minimum = minimum.min(quality(tet).powf(1.5));
                if minimum <= best_score + 1e-10 {
                    valid = false;
                    break;
                }
            }
            if valid && minimum > best_score + 1e-10 {
                best_score = minimum;
                best = Some((parameters, physical));
            }
        }
        if let Some((parameters, physical)) = best {
            for (j, &v) in group.iter().enumerate() {
                output.parameters[v] = parameters[j];
                output.mesh.points[v] = physical[j];
            }
            moved += group.len();
        }
    }
    Ok(moved)
}

fn orient_candidates(mut candidates: Vec<Tet>, old: &[Tet], points: &[Point]) -> Option<Vec<Tet>> {
    let old_min = old.iter().map(|&t| q(points, t)).fold(1., f64::min);
    let old_volume: f64 = old.iter().map(|t| determinant(t.map(|v| points[v]))).sum();
    let mut volume = 0.;
    let mut minimum: f64 = 1.;
    for t in &mut candidates {
        let det = determinant(t.map(|v| points[v]));
        if !det.is_finite() || det == 0. {
            return None;
        }
        if det < 0. {
            t.swap(2, 3);
        }
        volume += det.abs();
        minimum = minimum.min(q(points, *t));
    }
    let old_bad = old.iter().filter(|&&t| q(points, t) < 0.1).count();
    let new_bad = candidates.iter().filter(|&&t| q(points, t) < 0.1).count();
    (minimum > old_min + 1e-10
        && new_bad <= old_bad
        && (volume - old_volume).abs() <= old_volume * 1e-10)
        .then_some(candidates)
}

// Enumerate only the small (degree <= 7) edge rings, as in adaptive edge removal.
fn ring_candidates(a: usize, b: usize, old: &[Tet]) -> Vec<Vec<Tet>> {
    let mut link = BTreeMap::<usize, Vec<usize>>::new();
    for t in old {
        let rest: Vec<_> = t.iter().copied().filter(|&v| v != a && v != b).collect();
        if rest.len() != 2 {
            return Vec::new();
        }
        link.entry(rest[0]).or_default().push(rest[1]);
        link.entry(rest[1]).or_default().push(rest[0]);
    }
    if link.len() != old.len() || link.values().any(|v| v.len() != 2) {
        return Vec::new();
    }
    let mut ring = vec![*link.keys().next().unwrap()];
    while ring.len() < link.len() {
        let next = link[ring.last().unwrap()]
            .iter()
            .filter(|v| !ring.contains(v))
            .min();
        let Some(&next) = next else {
            return Vec::new();
        };
        ring.push(next);
    }
    if !link[ring.last().unwrap()].contains(&ring[0]) {
        return Vec::new();
    }
    fn triangles(r: &[usize]) -> Vec<Vec<[usize; 3]>> {
        if r.len() < 3 {
            return vec![Vec::new()];
        }
        let mut out = Vec::new();
        for k in 1..r.len() - 1 {
            for left in triangles(&r[..=k]) {
                for right in triangles(&r[k..]) {
                    let mut t = vec![[r[0], r[k], r[r.len() - 1]]];
                    t.extend_from_slice(&left);
                    t.extend_from_slice(&right);
                    out.push(t);
                }
            }
        }
        out
    }
    triangles(&ring)
        .into_iter()
        .map(|triangles| {
            triangles
                .into_iter()
                .flat_map(|[u, v, w]| [[u, v, w, a], [u, v, w, b]])
                .collect()
        })
        .collect()
}

#[derive(Clone)]
struct BoundaryPlan {
    faces: [[usize; 3]; 2],
    ids: [usize; 2],
    old: [Tet; 2],
    new: [Tet; 2],
    diagonal: [usize; 2],
    label: StratumMask,
}
fn flip_boundary(output: &mut Output, o: &Options) -> Result<usize, String> {
    let masks: Vec<_> = output
        .constraints
        .iter()
        .zip(&output.parameters)
        .map(|(&m, &p)| StratumMask::from(m) | cap(p, o))
        .collect();
    let mesh = &mut output.mesh;
    let n = mesh.points.len();
    let mut incident = vec![Vec::new(); n];
    for (i, t) in mesh.tets.iter().enumerate() {
        for &v in t {
            incident[v].push(i);
        }
    }
    let mut edges = BTreeMap::<[usize; 2], Vec<([usize; 3], usize)>>::new();
    for &f in &mesh.surface {
        let Some(&i) = incident[f[0]]
            .iter()
            .find(|&&i| f.iter().all(|v| mesh.tets[i].contains(v)))
        else {
            return Err("boundary face has no owner".into());
        };
        for [a, b] in [[f[0], f[1]], [f[0], f[2]], [f[1], f[2]]] {
            edges.entry([a.min(b), a.max(b)]).or_default().push((f, i));
        }
    }
    let normface = |mut f: [usize; 3]| {
        f.sort_unstable();
        f
    };
    let mut plans = BTreeMap::new();
    for ([a, b], adjacent) in edges {
        if adjacent.len() != 2 || (masks[a] & masks[b]).count_ones() != 1 {
            continue;
        }
        let [(f, i), (g, j)] = [adjacent[0], adjacent[1]];
        if i == j {
            continue;
        }
        let old = [mesh.tets[i], mesh.tets[j]];
        let oldmin = old.iter().map(|&t| q(&mesh.points, t)).fold(1., f64::min);
        let bit = (masks[a] & masks[b]).trailing_zeros() as usize;
        let periodic = bit >= MAX_CONSTRAINTS && o.periodic[(bit - MAX_CONSTRAINTS) / 2];
        if oldmin >= 0.1 || periodic {
            continue;
        }
        let star: BTreeSet<_> = incident[a]
            .iter()
            .copied()
            .filter(|&k| mesh.tets[k].contains(&b))
            .collect();
        if star != BTreeSet::from([i, j]) {
            continue;
        }
        let c = *f.iter().find(|&&v| v != a && v != b).unwrap();
        let d = *g.iter().find(|&&v| v != a && v != b).unwrap();
        let label = masks[a] & masks[b] & masks[c] & masks[d];
        if label.count_ones() != 1 || (masks[c] & masks[d]).count_ones() != 1 {
            continue;
        }
        if incident[c].iter().any(|&k| mesh.tets[k].contains(&d)) {
            continue;
        }
        let v = *old[0].iter().find(|v| !f.contains(v)).unwrap();
        if !old[1].contains(&v) || g.contains(&v) {
            continue;
        }
        let mut new = [[v, c, d, a], [v, c, d, b]];
        let det = new.map(|t| determinant(t.map(|v| mesh.points[v])));
        if det[0] * det[1] >= 0. || det.iter().any(|x| !x.is_finite()) {
            continue;
        }
        for (t, det) in new.iter_mut().zip(det) {
            if det < 0. {
                t.swap(0, 1);
            }
        }
        plans.insert(
            [a, b],
            BoundaryPlan {
                faces: [f, g],
                ids: [i, j],
                old,
                new,
                diagonal: [c, d],
                label,
            },
        );
    }
    let mut expected: BTreeSet<_> = mesh.surface.iter().copied().map(normface).collect();
    let mut used = BTreeSet::new();
    let mut accepted = 0;
    for plan in plans.values() {
        let bit = plan.label.trailing_zeros() as usize;
        let periodic = bit >= MAX_CONSTRAINTS && o.periodic[(bit - MAX_CONSTRAINTS) / 2];
        if periodic {
            continue;
        }
        let selected = vec![plan];
        let ids: BTreeSet<_> = selected.iter().flat_map(|p| p.ids).collect();
        if ids.len() != 2 * selected.len() || ids.iter().any(|i| used.contains(i)) {
            continue;
        }
        if selected.iter().any(|plan| {
            let [c, d] = plan.diagonal;
            incident[c].iter().any(|&k| mesh.tets[k].contains(&d))
        }) {
            continue;
        }
        let oldmin = selected
            .iter()
            .flat_map(|p| p.old)
            .map(|t| q(&mesh.points, t))
            .fold(1., f64::min);
        let minimum = selected
            .iter()
            .flat_map(|p| p.new)
            .map(|t| q(&mesh.points, t))
            .fold(1., f64::min);
        let oldbad = selected
            .iter()
            .flat_map(|p| p.old)
            .filter(|&t| q(&mesh.points, t) < 0.1)
            .count();
        let newbad = selected
            .iter()
            .flat_map(|p| p.new)
            .filter(|&t| q(&mesh.points, t) < 0.1)
            .count();
        if oldmin >= 0.1 || minimum <= oldmin + 1e-10 || newbad > oldbad {
            continue;
        }
        for plan in selected {
            for (&k, t) in plan.ids.iter().zip(plan.new) {
                mesh.tets[k] = t;
                for v in t {
                    incident[v].push(k);
                    incident[v].sort_unstable();
                    incident[v].dedup();
                }
            }
            for f in plan.faces {
                expected.remove(&normface(f));
            }
            let [c, d] = plan.diagonal;
            let edge = plan.faces[0]
                .iter()
                .copied()
                .filter(|v| plan.faces[1].contains(v))
                .collect::<Vec<_>>();
            expected.insert(normface([edge[0], c, d]));
            expected.insert(normface([edge[1], c, d]));
            used.extend(plan.ids);
            accepted += 1;
        }
    }
    let surface = crate::boundary(&mesh.tets)?;
    if surface
        .iter()
        .copied()
        .map(normface)
        .collect::<BTreeSet<_>>()
        != expected
        || !super::super::nonmanifold_vertices(&surface, n).is_empty()
    {
        return Err("boundary reconnection changed cavity interfaces or manifold".into());
    }
    mesh.surface = surface;
    Ok(accepted)
}

pub(super) fn flip(output: &mut Output, o: &Options) -> Result<usize, String> {
    let boundary = flip_boundary(output, o)?;
    let mesh = &mut output.mesh;
    let mut incident = vec![Vec::new(); mesh.points.len()];
    let faces = |t: Tet| -> [[usize; 3]; 4] {
        let mut f = [
            [t[0], t[1], t[2]],
            [t[0], t[1], t[3]],
            [t[0], t[2], t[3]],
            [t[1], t[2], t[3]],
        ];
        for tri in &mut f {
            tri.sort_unstable();
        }
        f
    };
    let mut face_ids = BTreeMap::<[usize; 3], Vec<usize>>::new();
    let mut edges = BTreeSet::new();
    for &t in &mesh.tets {
        if q(&mesh.points, t) < 0.1 {
            for f in faces(t) {
                face_ids.entry(f).or_default();
            }
            for a in 0..4 {
                for b in a + 1..4 {
                    edges.insert([t[a].min(t[b]), t[a].max(t[b])]);
                }
            }
        }
    }
    for (i, &t) in mesh.tets.iter().enumerate() {
        for v in t {
            incident[v].push(i);
        }
        for f in faces(t) {
            if let Some(ids) = face_ids.get_mut(&f) {
                ids.push(i);
            }
        }
    }
    let mut alive = vec![true; mesh.tets.len()];
    let mut accepted = 0;
    let mut live_count = mesh.tets.len();
    for (f, ids) in face_ids {
        if ids.len() != 2 || ids.iter().any(|&i| !alive[i]) {
            continue;
        }
        let a = *mesh.tets[ids[0]].iter().find(|v| !f.contains(v)).unwrap();
        let b = *mesh.tets[ids[1]].iter().find(|v| !f.contains(v)).unwrap();
        if incident[a]
            .iter()
            .any(|&i| alive[i] && mesh.tets[i].contains(&b))
        {
            continue;
        }
        let [u, v, w] = f;
        let candidates = vec![[a, b, u, v], [a, b, v, w], [a, b, w, u]];
        let Some(new) = orient_candidates(
            candidates,
            &[mesh.tets[ids[0]], mesh.tets[ids[1]]],
            &mesh.points,
        ) else {
            continue;
        };
        if live_count
            .saturating_add(new.len())
            .saturating_sub(ids.len())
            > o.max_tetrahedra
        {
            continue;
        }
        live_count = live_count + new.len() - ids.len();
        for i in ids {
            alive[i] = false;
        }
        for t in new {
            let i = mesh.tets.len();
            for v in t {
                incident[v].push(i);
            }
            mesh.tets.push(t);
            alive.push(true);
        }
        accepted += 1;
    }
    for [a, b] in edges {
        let ids: Vec<_> = incident[a]
            .iter()
            .copied()
            .filter(|&i| alive[i] && mesh.tets[i].contains(&b))
            .collect();
        if !(3..=7).contains(&ids.len()) {
            continue;
        }
        let old: Vec<_> = ids.iter().map(|&i| mesh.tets[i]).collect();
        let boundary = |ts: &[Tet]| {
            let mut counts = BTreeMap::<[usize; 3], usize>::new();
            for &t in ts {
                for f in faces(t) {
                    *counts.entry(f).or_default() += 1;
                }
            }
            counts
                .into_iter()
                .filter_map(|(f, n)| (n == 1).then_some(f))
                .collect::<BTreeSet<_>>()
        };
        let old_boundary = boundary(&old);
        let mut best = None;
        let mut best_min = old.iter().map(|&t| q(&mesh.points, t)).fold(1., f64::min);
        for candidate in ring_candidates(a, b, &old) {
            let conflict = candidate.iter().any(|t| {
                (0..4).any(|j| {
                    (j + 1..4).any(|k| {
                        let (u, v) = (t[j], t[k]);
                        !old.iter().any(|t| t.contains(&u) && t.contains(&v))
                            && incident[u]
                                .iter()
                                .any(|&i| alive[i] && mesh.tets[i].contains(&v))
                    })
                })
            });
            if conflict {
                continue;
            }
            let Some(new) = orient_candidates(candidate, &old, &mesh.points) else {
                continue;
            };
            let minimum = new.iter().map(|&t| q(&mesh.points, t)).fold(1., f64::min);
            if minimum <= best_min || boundary(&new) != old_boundary {
                continue;
            }
            let conflict = new.iter().flat_map(|&t| faces(t)).any(|face| {
                !old_boundary.contains(&face)
                    && incident[face[0]].iter().any(|&i| {
                        alive[i]
                            && !ids.contains(&i)
                            && face.iter().all(|v| mesh.tets[i].contains(v))
                    })
            });
            if conflict || crate::boundary(&new).is_err() {
                continue;
            }
            best_min = minimum;
            best = Some(new);
        }
        let Some(new) = best else {
            continue;
        };
        if live_count
            .saturating_add(new.len())
            .saturating_sub(old.len())
            > o.max_tetrahedra
        {
            continue;
        }
        live_count = live_count + new.len() - old.len();
        for i in ids {
            alive[i] = false;
        }
        for t in new {
            let i = mesh.tets.len();
            for v in t {
                incident[v].push(i);
            }
            mesh.tets.push(t);
            alive.push(true);
        }
        accepted += 1;
    }
    mesh.tets = mesh
        .tets
        .iter()
        .enumerate()
        .filter_map(|(i, &t)| alive[i].then_some(t))
        .collect();
    let normalize = |mut faces: Vec<[usize; 3]>| {
        for f in &mut faces {
            f.sort_unstable();
        }
        faces.sort_unstable();
        faces
    };
    let surface = crate::boundary(&mesh.tets)?;
    if normalize(surface.clone()) != normalize(mesh.surface.clone()) {
        return Err("local reconnection changed boundary triangulation".into());
    }
    mesh.surface = surface;
    Ok(accepted + boundary)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rejects_creating_a_shortcut_between_separate_rims() {
        let masks = [1, 5, 5, 1, 1];
        assert!(!feature_link_ok(0, 1, &[[0, 2, 3], [1, 3, 4]], &masks));
    }
    #[test]
    fn permits_removing_an_interior_point_of_one_rim() {
        let masks = [5, 5, 5, 5, 1, 4];
        let faces = [
            [2, 0, 4],
            [2, 5, 0],
            [0, 1, 4],
            [0, 5, 1],
            [1, 3, 4],
            [1, 5, 3],
        ];
        assert!(feature_link_ok(0, 1, &faces, &masks));
    }
    #[test]
    fn forbids_collapsing_a_three_edge_feature_loop() {
        let masks = [5, 5, 5, 1, 4];
        let faces = [
            [0, 1, 3],
            [1, 2, 3],
            [2, 0, 3],
            [1, 0, 4],
            [2, 1, 4],
            [0, 2, 4],
        ];
        assert!(!feature_link_ok(0, 1, &faces, &masks));
    }
    #[cfg(feature = "wide-constraints")]
    #[test]
    fn high_constraint_and_cap_bits_preserve_feature_links() {
        let wall = 1 << 31;
        let cap = 1 << (MAX_CONSTRAINTS + 5);
        let masks = [wall | cap, wall | cap, wall | cap, wall | cap, wall, cap];
        let faces = [
            [2, 0, 4],
            [2, 5, 0],
            [0, 1, 4],
            [0, 5, 1],
            [1, 3, 4],
            [1, 5, 3],
        ];
        assert!(feature_link_ok(0, 1, &faces, &masks));
        let masks = [wall, wall | cap, wall | cap, wall, wall];
        assert!(!feature_link_ok(0, 1, &[[0, 2, 3], [1, 3, 4]], &masks));
    }
}
