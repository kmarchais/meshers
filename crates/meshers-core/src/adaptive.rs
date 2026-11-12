//! Conforming surface-edge bisection, applied to whole periodic edge orbits.
use crate::{Mesh, Point, boundary, determinant, dot, norm, quality, sub};
use std::collections::{BTreeMap, HashMap, HashSet};
type Edge = [usize; 2];
fn edge(a: usize, b: usize) -> Edge {
    if a < b { [a, b] } else { [b, a] }
}

/// One nested pass over disjoint interior edge stars, ranked by squared length.
/// Boundary coordinates and triangles remain unchanged. The node budget is an upper limit.
pub(crate) fn refine_interior(
    mesh: &mut Mesh,
    node_budget: usize,
    minimum_quality: f64,
    max_tets: usize,
    geometry: &impl crate::geometry::Geometry,
) -> Result<usize, String> {
    if node_budget == 0 || mesh.tets.is_empty() {
        return Ok(0);
    }
    geometry.check()?;
    let boundary_edges: HashSet<_> = mesh
        .surface
        .iter()
        .flat_map(|f| [edge(f[0], f[1]), edge(f[1], f[2]), edge(f[2], f[0])])
        .collect();
    let mut incidence = Vec::with_capacity(mesh.tets.len() * 6);
    for (ci, t) in mesh.tets.iter().enumerate() {
        for i in 0..4 {
            for j in i + 1..4 {
                incidence.push((edge(t[i], t[j]), ci));
            }
        }
    }
    incidence.sort_unstable();
    geometry.check()?;
    let [lo, hi] = geometry.bounds();
    let midpoint = |e: Edge| {
        std::array::from_fn::<_, 3, _>(|a| (mesh.points[e[0]][a] + mesh.points[e[1]][a]) * 0.5)
    };
    let mut candidates = Vec::new();
    let mut start = 0;
    while start < incidence.len() {
        let e = incidence[start].0;
        let mut end = start + 1;
        while end < incidence.len() && incidence[end].0 == e {
            end += 1;
        }
        let p = midpoint(e);
        let delta = sub(mesh.points[e[0]], mesh.points[e[1]]);
        let length2 = dot(delta, delta);
        if !boundary_edges.contains(&e)
            && (0..3).all(|a| {
                p[a] > lo[a] + (hi[a] - lo[a]) * 1e-9 && p[a] < hi[a] - (hi[a] - lo[a]) * 1e-9
            })
            && length2 > 1e-20
        {
            candidates.push((e, start, end, length2));
        }
        start = end;
    }
    candidates.sort_unstable_by(|a, b| b.3.total_cmp(&a.3).then(a.0.cmp(&b.0)));
    let mut removed = vec![false; mesh.tets.len()];
    let mut additions = Vec::new();
    let original_nodes = mesh.points.len();
    let mut output_tets = mesh.tets.len();
    for (index, (e, start, end, _)) in candidates.into_iter().enumerate() {
        if mesh.points.len() - original_nodes >= node_budget {
            break;
        }
        if index % 1024 == 0 {
            geometry.check()?;
        }
        let star = &incidence[start..end];
        if star.iter().any(|&(_, ci)| removed[ci])
            || output_tets.saturating_add(star.len()) > max_tets
        {
            continue;
        }
        let p = std::array::from_fn(|a| (mesh.points[e[0]][a] + mesh.points[e[1]][a]) * 0.5);
        let id = mesh.points.len();
        let mut children = Vec::with_capacity(star.len() * 2);
        for &(_, ci) in star {
            for endpoint in e {
                let mut t = mesh.tets[ci];
                *t.iter_mut().find(|v| **v == endpoint).unwrap() = id;
                let xyz = t.map(|i| if i == id { p } else { mesh.points[i] });
                let det = determinant(xyz);
                let edge_sum: f64 = (0..4)
                    .flat_map(|i| {
                        (i + 1..4).map(move |j| {
                            let delta = sub(xyz[i], xyz[j]);
                            dot(delta, delta)
                        })
                    })
                    .sum();
                let q = (432. * det * det / edge_sum.powi(3)).sqrt();
                if !det.is_finite() || det <= 0. || !q.is_finite() || q < minimum_quality {
                    children.clear();
                    break;
                }
                children.push(t);
            }
            if children.is_empty() {
                break;
            }
        }
        if children.len() == star.len() * 2 {
            mesh.points.push(p);
            for &(_, ci) in star {
                removed[ci] = true;
            }
            additions.extend(children);
            output_tets += star.len();
        }
    }
    mesh.tets = mesh
        .tets
        .iter()
        .enumerate()
        .filter_map(|(i, &t)| (!removed[i]).then_some(t))
        .chain(additions)
        .collect();
    geometry.check()?;
    if mesh.surface != boundary(&mesh.tets)? {
        return Err("interior refinement changed the oriented boundary".into());
    }
    Ok(mesh.points.len() - original_nodes)
}
fn orbit(mesh: &Mesh, e: Edge, geometry: &impl crate::geometry::Geometry) -> [[i64; 3]; 2] {
    let mut k = e.map(|i| geometry.rank(mesh.points[i]));
    k.sort_unstable();
    k
}
fn error(p: Point, t: f64, geometry: &impl crate::geometry::Geometry) -> f64 {
    geometry.residual(p, t).abs() / norm(geometry.gradient(p))
}
fn midpoint(
    mesh: &Mesh,
    e: Edge,
    t: f64,
    geometry: &impl crate::geometry::Geometry,
) -> Option<Point> {
    let [a, b] = e.map(|i| mesh.points[i]);
    let mut p = std::array::from_fn(|i| (a[i] + b[i]) * 0.5);
    let level = geometry.level(geometry.value(a), t);
    for _ in 0..12 {
        let value = geometry.value(p) - level;
        if value.abs() < 1e-13 {
            break;
        }
        let mut g = geometry.gradient(p);
        for i in 0..3 {
            if a[i] == b[i] && geometry.locked(a)[i] {
                g[i] = 0.;
            }
        }
        let g2 = dot(g, g);
        if g2 < 1e-16 {
            return None;
        }
        for i in 0..3 {
            p[i] -= value * g[i] / g2;
        }
    }
    ((geometry.value(p) - level).abs() < 1e-11 && geometry.contains(p)).then_some(p)
}
pub fn refine(
    mesh: &mut Mesh,
    threshold: f64,
    tolerance: f64,
    max_tets: usize,
) -> Result<usize, String> {
    refine_for(
        mesh,
        threshold,
        tolerance,
        max_tets,
        &crate::geometry::Legacy,
    )
}
pub(crate) fn refine_for(
    mesh: &mut Mesh,
    threshold: f64,
    tolerance: f64,
    max_tets: usize,
    geometry: &impl crate::geometry::Geometry,
) -> Result<usize, String> {
    if tolerance == 0. {
        return Ok(0);
    }
    let mut accepted = 0;
    for _ in 0..3 {
        geometry.check()?;
        let mut profile = crate::profile::Profile::new("adaptive_round");
        profile.count("input_tetrahedra", mesh.tets.len());
        let mut scores = BTreeMap::<_, f64>::new();
        let mut surface_edges = BTreeMap::<_, Vec<Edge>>::new();
        let mut batched_errors = Vec::new();
        if geometry.batched() {
            let mut samples = Vec::new();
            for f in &mesh.surface {
                let p = f.map(|i| mesh.points[i]);
                if geometry.cap(p) {
                    continue;
                }
                samples.push(std::array::from_fn(|a| (p[0][a] + p[1][a] + p[2][a]) / 3.));
                for [i, j] in [[0, 1], [1, 2], [2, 0]] {
                    // Match the scalar path's index-sorted edge arithmetic.
                    let e = edge(f[i], f[j]);
                    samples.push(std::array::from_fn(|a| {
                        (mesh.points[e[0]][a] + mesh.points[e[1]][a]) * 0.5
                    }));
                }
            }
            for chunk in samples.chunks(4096) {
                geometry.check()?;
                let values = geometry.values(chunk);
                let gradients = geometry.gradients(chunk);
                batched_errors.extend(values.into_iter().zip(gradients).map(|(v, g)| {
                    (if geometry.band() {
                        v.abs() - threshold
                    } else {
                        v
                    })
                    .abs()
                        / norm(g)
                }));
            }
        }
        let mut batched_errors = batched_errors.into_iter();
        for f in &mesh.surface {
            if geometry.cap(f.map(|i| mesh.points[i])) {
                continue;
            }
            let p = f.map(|i| mesh.points[i]);
            let center = std::array::from_fn(|a| (p[0][a] + p[1][a] + p[2][a]) / 3.);
            let center_error = if geometry.batched() {
                batched_errors.next().unwrap()
            } else {
                error(center, threshold, geometry)
            };
            let face_edges = [edge(f[0], f[1]), edge(f[1], f[2]), edge(f[2], f[0])];
            let longest = *face_edges
                .iter()
                .max_by(|a, b| {
                    norm(sub(mesh.points[a[0]], mesh.points[a[1]]))
                        .total_cmp(&norm(sub(mesh.points[b[0]], mesh.points[b[1]])))
                })
                .unwrap();
            for e in face_edges {
                let key = orbit(mesh, e, geometry);
                surface_edges.entry(key).or_default().push(e);
                let mid =
                    std::array::from_fn(|a| (mesh.points[e[0]][a] + mesh.points[e[1]][a]) * 0.5);
                let score = (if geometry.batched() {
                    batched_errors.next().unwrap()
                } else {
                    error(mid, threshold, geometry)
                })
                .max(if e == longest { center_error * 0.8 } else { 0. });
                if score > tolerance {
                    scores
                        .entry(key)
                        .and_modify(|v| *v = v.max(score))
                        .or_insert(score);
                }
            }
        }
        profile.mark("surface_scoring");
        profile.count("candidate_orbits", scores.len());
        let mut candidates: Vec<_> = scores.into_iter().collect();
        candidates.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
        profile.mark("candidate_sort");
        let mut edges: HashMap<Edge, Vec<usize>> = HashMap::new();
        for (key, _) in &candidates {
            for &e in &surface_edges[key] {
                edges.entry(e).or_default();
            }
        }
        for (ci, t) in mesh.tets.iter().enumerate() {
            for i in 0..4 {
                for j in i + 1..4 {
                    if let Some(cells) = edges.get_mut(&edge(t[i], t[j])) {
                        cells.push(ci);
                    }
                }
            }
        }
        profile.mark("edge_incidence");
        profile.count("unique_edges", edges.len());

        let mut removed = HashSet::new();
        let mut additions = Vec::new();
        let mut changed = 0;
        for (key, _) in candidates {
            let mut copies = surface_edges[&key].clone();
            copies.sort_unstable();
            copies.dedup();
            let cells: Vec<_> = copies
                .iter()
                .flat_map(|e| edges[e].iter().copied())
                .collect();
            if cells.iter().any(|ci| removed.contains(ci))
                || mesh.tets.len() + additions.len() - removed.len() + cells.len() > max_tets
            {
                continue;
            }
            let mut new_points = Vec::new();
            let mut new_tets = Vec::new();
            let mut valid = true;
            for &e in &copies {
                let Some(p) = midpoint(mesh, e, threshold, geometry) else {
                    valid = false;
                    break;
                };
                let id = mesh.points.len() + new_points.len();
                new_points.push(p);
                for &ci in &edges[&e] {
                    for endpoint in e {
                        let mut t = mesh.tets[ci];
                        *t.iter_mut().find(|v| **v == endpoint).unwrap() = id;
                        let xyz = t.map(|i| if i == id { p } else { mesh.points[i] });
                        if determinant(xyz) <= 0.
                            || quality(xyz).powf(1.5) < geometry.refinement_quality_floor()
                        {
                            valid = false;
                            break;
                        }
                        new_tets.push(t);
                    }
                    if !valid {
                        break;
                    }
                }
                if !valid {
                    break;
                }
            }
            if valid {
                removed.extend(cells);
                mesh.points.extend(new_points);
                additions.extend(new_tets);
                changed += 1;
            }
        }
        profile.mark("split_validation");
        profile.count("accepted_orbits", changed);
        if changed == 0 {
            break;
        }
        let deleted: Vec<_> = mesh
            .tets
            .iter()
            .enumerate()
            .filter_map(|(i, t)| removed.contains(&i).then_some(*t))
            .collect();
        let surface = update_boundary(&mesh.surface, &deleted, &additions)?;
        mesh.tets = mesh
            .tets
            .iter()
            .enumerate()
            .filter_map(|(ci, &t)| (!removed.contains(&ci)).then_some(t))
            .chain(additions)
            .collect();
        profile.mark("tetrahedron_compaction");
        mesh.surface = surface;
        if std::env::var_os("MESHER_CHECK_TOPOLOGY").is_some()
            && mesh.surface != boundary(&mesh.tets)?
        {
            return Err("incremental boundary disagrees with full reconstruction".into());
        }
        profile.mark("boundary_rebuild");
        accepted += changed;
    }
    // Retain a full final manifold/orientation audit of the volume topology.
    if mesh.surface != boundary(&mesh.tets)? {
        return Err("incremental boundary disagrees with final volume topology".into());
    }
    Ok(accepted)
}

// Boundary(new volume) = boundary(old volume) - boundary(deleted) + boundary(added).
fn update_boundary(
    surface: &[[usize; 3]],
    deleted: &[crate::Tet],
    added: &[crate::Tet],
) -> Result<Vec<[usize; 3]>, String> {
    let mut faces: Vec<_> = surface
        .iter()
        .copied()
        .map(|f| (crate::face_key(f), f, 0u8))
        .chain(
            deleted
                .iter()
                .flat_map(crate::tet_faces)
                .map(|f| (crate::face_key(f), f, 1)),
        )
        .chain(
            added
                .iter()
                .flat_map(crate::tet_faces)
                .map(|f| (crate::face_key(f), f, 2)),
        )
        .collect();
    faces.sort_unstable_by_key(|v| v.0);
    let mut result = Vec::new();
    let mut i = 0;
    while i < faces.len() {
        let mut j = i + 1;
        while j < faces.len() && faces[j].0 == faces[i].0 {
            j += 1;
        }
        let signed = |f: [usize; 3]| if crate::face_parity(f) == 0 { 1i32 } else { -1 };
        let sum: i32 = faces[i..j]
            .iter()
            .map(|(_, f, kind)| signed(*f) * if *kind == 1 { -1 } else { 1 })
            .sum();
        if sum.abs() > 1 {
            return Err("invalid incremental boundary multiplicity".into());
        }
        if sum != 0 {
            let chosen = faces[i..j]
                .iter()
                .filter(|(_, f, k)| *k != 1 && signed(*f) == sum)
                .max_by_key(|(_, _, k)| *k)
                .ok_or("incremental boundary has no surviving face")?;
            result.push(chosen.1);
        }
        i = j;
    }
    result.sort_unstable();
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn interior_quality_rejection_leaves_mesh_unchanged() {
        let mut mesh = crate::implicit::generate(
            &|_| -1.,
            crate::implicit::Options {
                cells: [4; 3],
                optimize_passes: 0,
                ..Default::default()
            },
        )
        .unwrap()
        .mesh;
        let points = mesh.points.clone();
        let tets = mesh.tets.clone();
        assert_eq!(
            refine_interior(&mut mesh, 20, 1., 10000, &crate::geometry::Legacy).unwrap(),
            0
        );
        assert_eq!(mesh.points, points);
        assert_eq!(mesh.tets, tets);
    }
    #[test]
    fn incremental_boundary_matches_full_edge_bisection() {
        let old = [[0, 1, 2, 3]];
        let new = [[0, 4, 2, 3], [4, 1, 2, 3]];
        assert_eq!(
            update_boundary(&boundary(&old).unwrap(), &old, &new).unwrap(),
            boundary(&new).unwrap()
        );
    }
}
