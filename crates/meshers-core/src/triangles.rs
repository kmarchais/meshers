//! Experimental direct triangle surface extraction on a Freudenthal sampling grid.
//! Only crossing edges and boundary polygons are retained. No volume mesh is built.
use crate::{MeshingError, Point, cross, dot, implicit::ScalarField, sub, surface_band::Band};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Serialize, Deserialize)]
pub struct TriangleMesh {
    pub points: Vec<Point>,
    pub faces: Vec<[usize; 3]>,
    pub labels: Vec<u8>,
}
type Result<T> = std::result::Result<T, MeshingError>;
fn fail(message: &str) -> MeshingError {
    MeshingError::GenerationFailed(message.into())
}
struct Extractor<'a, 'b, F: ScalarField> {
    band: &'a Band<'b, F>,
    grid: Vec<Point>,
    values: Vec<f64>,
    mesh: TriangleMesh,
    cache: HashMap<(usize, usize, u8), usize>,
}
impl<F: ScalarField> Extractor<'_, '_, F> {
    fn original(&mut self, id: usize) -> usize {
        *self.cache.entry((id, id, 8)).or_insert_with(|| {
            let out = self.mesh.points.len();
            self.mesh.points.push(self.grid[id]);
            out
        })
    }
    fn root(&mut self, a: usize, b: usize, label: u8) -> Result<usize> {
        let (a, b) = if a < b { (a, b) } else { (b, a) };
        if let Some(&id) = self.cache.get(&(a, b, label)) {
            return Ok(id);
        }
        let level = self.band.levels[if label == 0 { 1 } else { 0 }];
        let va = self.values[a] - level;
        let vb = self.values[b] - level;
        if va == 0.0 {
            return Ok(self.original(a));
        }
        if vb == 0.0 {
            return Ok(self.original(b));
        }
        if va.signum() == vb.signum() {
            return Err(fail("No bracket on surface edge"));
        }
        let mut lo = 0.;
        let mut hi = 1.;
        let mut f_lo = va;
        let mut f_hi = vb;
        let mut t = 0.5;
        for iteration in 0..44 {
            let midpoint = (lo + hi) * 0.5;
            t = if iteration < 12 {
                let secant = (lo * f_hi - hi * f_lo) / (f_hi - f_lo);
                if secant.is_finite()
                    && secant > lo + 1e-6 * (hi - lo)
                    && secant < hi - 1e-6 * (hi - lo)
                {
                    secant
                } else {
                    midpoint
                }
            } else {
                midpoint
            };
            let p = std::array::from_fn(|k| self.grid[a][k] * (1. - t) + self.grid[b][k] * t);
            let v = self.band.field.value(p) - level;
            if !v.is_finite() {
                return Err(fail("Non-finite field during edge solve"));
            }
            if v.abs() <= 1e-12 * (1. + va.abs().max(vb.abs())) {
                break;
            }
            if v.signum() == f_lo.signum() {
                lo = t;
                f_lo = v;
            } else {
                hi = t;
                f_hi = v;
            }
        }
        let p = std::array::from_fn(|k| self.grid[a][k] * (1. - t) + self.grid[b][k] * t);
        let id = self.mesh.points.len();
        self.mesh.points.push(p);
        self.cache.insert((a, b, label), id);
        Ok(id)
    }
    fn polygon(&mut self, mut ids: Vec<usize>, label: u8) -> Result<()> {
        ids.sort_unstable();
        ids.dedup();
        if ids.len() < 3 {
            return Ok(());
        }
        let center = std::array::from_fn(|k| {
            ids.iter().map(|&i| self.mesh.points[i][k]).sum::<f64>() / ids.len() as f64
        });
        let normal = self.band.normal(center, label);
        let axis = sub(self.mesh.points[ids[0]], center);
        let perpendicular = cross(normal, axis);
        ids.sort_by(|&a, &b| {
            let a = sub(self.mesh.points[a], center);
            let b = sub(self.mesh.points[b], center);
            dot(a, perpendicular)
                .atan2(dot(a, axis))
                .total_cmp(&dot(b, perpendicular).atan2(dot(b, axis)))
        });
        // Choose the fan with the best worst triangle, avoiding a fixed diagonal bias.
        let canonical_faces = |faces: &[[usize; 3]]| {
            let mut triangles: Vec<_> = faces
                .iter()
                .map(|face| {
                    let mut vertices = face.map(|i| {
                        let mut p = self.mesh.points[i];
                        p[((label - 2) / 2) as usize] = 0.;
                        p.map(|x| (x * 1e10).round() as i64)
                    });
                    vertices.sort_unstable();
                    vertices
                })
                .collect();
            triangles.sort_unstable();
            triangles
        };
        let mut best = None;
        for start in 0..ids.len() {
            let mut faces = Vec::new();
            let mut worst = f64::INFINITY;
            for j in 1..ids.len() - 1 {
                let face = [
                    ids[start],
                    ids[(start + j) % ids.len()],
                    ids[(start + j + 1) % ids.len()],
                ];
                let [a, b, c] = face.map(|i| self.mesh.points[i]);
                let n = self
                    .band
                    .normal(std::array::from_fn(|k| (a[k] + b[k] + c[k]) / 3.), label);
                let e = sub(b, a);
                let d = sub(c, a);
                let edge = sub(c, b);
                let score = dot(cross(e, d), n) / (dot(e, e) + dot(d, d) + dot(edge, edge));
                worst = worst.min(score);
                faces.push(face);
            }
            let key = if label >= 2 {
                canonical_faces(&faces)
            } else {
                Vec::new()
            };
            if best
                .as_ref()
                .is_none_or(|(score, old_key, _): &(f64, Vec<[[i64; 3]; 3]>, _)| {
                    if label >= 2 {
                        worst > *score + 1e-12
                            || ((worst - *score).abs() <= 1e-12 && key < *old_key)
                    } else {
                        worst > *score
                    }
                })
            {
                best = Some((worst, key, faces))
            }
        }
        let (score, _, faces) = best.unwrap();
        if score <= 0. {
            return Err(fail("Inverted triangle after analytic edge intersection"));
        }
        for face in faces {
            self.mesh.faces.push(face);
            self.mesh.labels.push(label)
        }
        Ok(())
    }
    fn cap(&mut self, face: [usize; 3], label: u8) -> Result<()> {
        // Each clipped point retains its original grid-edge support so both
        // implicit walls and clipping caps share exactly the same vertex IDs.
        let mut poly = face.map(|id| ([id, id], self.values[id], None)).to_vec();
        for wall in [1u8, 0] {
            let level = self.band.levels[if wall == 0 { 1 } else { 0 }];
            let inside = |v: f64| if wall == 0 { v <= level } else { v >= level };
            let mut next = Vec::new();
            for j in 0..poly.len() {
                let a = poly[j];
                let b = poly[(j + 1) % poly.len()];
                if inside(a.1) {
                    next.push(a)
                }
                if inside(a.1) != inside(b.1) {
                    let mut support = vec![a.0[0], a.0[1], b.0[0], b.0[1]];
                    support.sort_unstable();
                    support.dedup();
                    if support.len() != 2 {
                        return Err(fail("Invalid clipping edge support"));
                    }
                    let id = self.root(support[0], support[1], wall)?;
                    next.push(([support[0], support[1]], level, Some(id)));
                }
            }
            poly = next;
        }
        let ids = poly
            .into_iter()
            .map(|(support, _, id)| id.unwrap_or_else(|| self.original(support[0])))
            .collect();
        self.polygon(ids, label)
    }
}

pub fn extract<F: ScalarField>(band: &Band<'_, F>, cells: usize) -> Result<TriangleMesh> {
    band.validate()?;
    if !(4..=128).contains(&cells) {
        return Err(MeshingError::InvalidOptions(
            "Surface cells must be in 4..=128".into(),
        ));
    }
    let n = cells + 1;
    let mut grid = Vec::with_capacity(n * n * n);
    let mut values = Vec::with_capacity(n * n * n);
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let p = std::array::from_fn(|k| {
                    band.bounds[0][k]
                        + (band.bounds[1][k] - band.bounds[0][k]) * [x, y, z][k] as f64
                            / cells as f64
                });
                let value = band.field.value(p);
                if !value.is_finite() {
                    return Err(fail("Non-finite sampled field"));
                }
                grid.push(p);
                values.push(value);
            }
        }
    }
    let mut e = Extractor {
        band,
        grid,
        values,
        mesh: TriangleMesh {
            points: Vec::new(),
            faces: Vec::new(),
            labels: Vec::new(),
        },
        cache: HashMap::new(),
    };
    let stride = [1, n, n * n];
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                let base = x + n * (y + n * z);
                for perm in [
                    [0, 1, 2],
                    [0, 2, 1],
                    [1, 0, 2],
                    [1, 2, 0],
                    [2, 0, 1],
                    [2, 1, 0],
                ] {
                    let tet = [
                        base,
                        base + stride[perm[0]],
                        base + stride[perm[0]] + stride[perm[1]],
                        base + 1 + n + n * n,
                    ];
                    for label in 0..2 {
                        let level = band.levels[if label == 0 { 1 } else { 0 }];
                        let mut ids = Vec::new();
                        for i in 0..4 {
                            for j in i + 1..4 {
                                if (e.values[tet[i]] < level) != (e.values[tet[j]] < level) {
                                    ids.push(e.root(tet[i], tet[j], label)?)
                                }
                            }
                        }
                        e.polygon(ids, label)?;
                    }
                    for omit in 0..4 {
                        let face: Vec<_> = (0..4).filter(|&j| j != omit).map(|j| tet[j]).collect();
                        for label in 2..8 {
                            let axis = (label - 2) / 2;
                            let boundary_index = if (label - 2) % 2 == 0 { 0 } else { cells };
                            if face
                                .iter()
                                .all(|&i| (i / stride[axis]) % n == boundary_index)
                            {
                                e.cap([face[0], face[1], face[2]], label as u8)?;
                            }
                        }
                    }
                }
            }
        }
    }
    // Avoid a silent success for unresolved or empty geometry in the benchmark.
    if e.mesh.faces.is_empty() {
        return Err(fail("No resolved triangle surface"));
    }
    band.field.check().map_err(MeshingError::GenerationFailed)?;
    Ok(e.mesh)
}

/// Constrained neighbor averaging with local orientation checks. Connectivity
/// stays fixed; a rejected move leaves the original vertex unchanged.
pub fn smooth<F: ScalarField>(
    mesh: &mut TriangleMesh,
    band: &Band<'_, F>,
    iterations: usize,
) -> Result<()> {
    smooth_impl(mesh, band, iterations, false)
}

fn shape<F: ScalarField>(points: [Point; 3], label: u8, band: &Band<'_, F>) -> f64 {
    let [a, b, c] = points;
    let normal = band.normal(std::array::from_fn(|k| (a[k] + b[k] + c[k]) / 3.0), label);
    shape_with_normal(points, normal)
}

fn shape_with_normal(points: [Point; 3], normal: Point) -> f64 {
    let [a, b, c] = points;
    let e = sub(b, a);
    let d = sub(c, a);
    let other = sub(c, b);
    2.0 * 3.0f64.sqrt() * dot(cross(e, d), normal)
        / (dot(e, e) + dot(d, d) + dot(other, other)).max(1e-300)
}

/// Collapse short sliver edges and improve patch connectivity without adding nodes.
/// Local flips improve the worst triangle shape, or reduce excess valence while
/// preserving a shape floor. Constrained relaxation then preserves weak cells.
pub fn improve<F: ScalarField>(
    mesh: &mut TriangleMesh,
    band: &Band<'_, F>,
    rounds: usize,
) -> Result<usize> {
    smooth_impl(mesh, band, 0, true)?;
    let mut edits = 0;
    for _ in 0..rounds {
        collapse_slivers(mesh, band)?;
        let mut edges = std::collections::BTreeMap::<(usize, usize), Vec<(usize, usize)>>::new();
        let mut degree = vec![0i32; mesh.points.len()];
        for (fi, f) in mesh.faces.iter().enumerate() {
            for j in 0..3 {
                let a = f[j];
                let b = f[(j + 1) % 3];
                edges.entry((a.min(b), a.max(b))).or_default().push((fi, j));
            }
        }
        for &(a, b) in edges.keys() {
            degree[a] += 1;
            degree[b] += 1;
        }
        let mut touched = vec![false; mesh.points.len()];
        for (&(a, b), adjacent) in &edges {
            if adjacent.len() != 2 {
                return Err(fail("Expected closed triangle mesh"));
            }
            let (fi, j) = adjacent[0];
            let (gi, k) = adjacent[1];
            if mesh.labels[fi] != mesh.labels[gi] {
                continue;
            }
            let f = mesh.faces[fi];
            let g = mesh.faces[gi];
            let u = f[j];
            let v = f[(j + 1) % 3];
            let c = f[(j + 2) % 3];
            let d = g[(k + 2) % 3];
            if [a, b, c, d].iter().any(|&i| touched[i])
                || c == d
                || edges.contains_key(&(c.min(d), c.max(d)))
            {
                continue;
            }
            let next = [[c, d, v], [d, c, u]];
            let before = shape(f.map(|i| mesh.points[i]), mesh.labels[fi], band).min(shape(
                g.map(|i| mesh.points[i]),
                mesh.labels[gi],
                band,
            ));
            let after = next
                .iter()
                .map(|f| shape(f.map(|i| mesh.points[i]), mesh.labels[fi], band))
                .fold(f64::INFINITY, f64::min);
            let gain = [a, b, c, d]
                .iter()
                .map(|&i| (degree[i] - 6).pow(2))
                .sum::<i32>()
                - [(a, -1), (b, -1), (c, 1), (d, 1)]
                    .iter()
                    .map(|&(i, delta)| (degree[i] + delta - 6).pow(2))
                    .sum::<i32>();
            if after > before + 0.01 || (gain > 0 && after >= before.max(0.5)) {
                mesh.faces[fi] = next[0];
                mesh.faces[gi] = next[1];
                for i in [a, b, c, d] {
                    touched[i] = true;
                }
                edits += 1;
            }
        }
        smooth_impl(mesh, band, 2, true)?;
    }
    Ok(edits)
}

/// Redistribute existing feature vertices using only their two curve neighbours.
/// Triple intersections stay fixed. Projection and a local shape floor protect
/// both incident patches; connectivity and the number of vertices are unchanged.
pub fn redistribute_features<F: ScalarField>(
    mesh: &mut TriangleMesh,
    band: &Band<'_, F>,
    passes: usize,
) -> Result<usize> {
    smooth_impl(mesh, band, 0, true)?;
    if passes == 0 {
        return Ok(0);
    }
    let mut incident = vec![Vec::new(); mesh.points.len()];
    let mut masks = vec![0u8; mesh.points.len()];
    let mut edges = std::collections::BTreeMap::<(usize, usize), u8>::new();
    for (fi, f) in mesh.faces.iter().enumerate() {
        for j in 0..3 {
            incident[f[j]].push(fi);
            masks[f[j]] |= 1 << mesh.labels[fi];
            let (a, b) = (f[j], f[(j + 1) % 3]);
            *edges.entry((a.min(b), a.max(b))).or_default() |= 1 << mesh.labels[fi];
        }
    }
    let mut neighbors = vec![Vec::new(); mesh.points.len()];
    for ((a, b), mask) in edges {
        if mask.count_ones() == 2 {
            neighbors[a].push(b);
            neighbors[b].push(a);
        }
    }
    let length = |a: Point, b: Point| {
        let d = sub(a, b);
        dot(d, d).sqrt()
    };
    let mut moved = 0;
    for _ in 0..passes {
        let mut changed = false;
        for i in 0..mesh.points.len() {
            if masks[i].count_ones() != 2 || neighbors[i].len() != 2 {
                continue;
            }
            let a = mesh.points[neighbors[i][0]];
            let b = mesh.points[neighbors[i][1]];
            let old = mesh.points[i];
            let before = (length(old, a) - length(old, b)).abs();
            let floor = incident[i]
                .iter()
                .map(|&fi| {
                    shape(
                        mesh.faces[fi].map(|v| mesh.points[v]),
                        mesh.labels[fi],
                        band,
                    )
                })
                .fold(0.45, f64::min);
            let target: Point = std::array::from_fn(|k| 0.5 * (a[k] + b[k]));
            for attempt in 0..10 {
                let weight = 0.5f64.powi(attempt);
                let trial = std::array::from_fn(|k| old[k] + weight * (target[k] - old[k]));
                let Ok(trial) = band.project(trial, masks[i]) else {
                    continue;
                };
                let value = band.field.value(trial);
                if masks[i] & 3 == 0 && (value < band.levels[0] || value > band.levels[1]) {
                    continue;
                }
                if (length(trial, a) - length(trial, b)).abs() >= before * 0.999 {
                    continue;
                }
                if incident[i].iter().all(|&fi| {
                    shape(
                        mesh.faces[fi].map(|v| if v == i { trial } else { mesh.points[v] }),
                        mesh.labels[fi],
                        band,
                    ) >= floor.max(1e-10)
                }) {
                    mesh.points[i] = trial;
                    moved += 1;
                    changed = true;
                    break;
                }
            }
        }
        if !changed {
            break;
        }
    }
    band.field.check().map_err(MeshingError::GenerationFailed)?;
    Ok(moved)
}

/// Optimize the weak triangle patches while holding analytic feature constraints.
/// Numerical local derivatives keep this research objective independent of the
/// field's optional gradient and never alter connectivity.
pub fn polish<F: ScalarField>(
    mesh: &mut TriangleMesh,
    band: &Band<'_, F>,
    passes: usize,
) -> Result<()> {
    smooth_impl(mesh, band, 0, true)?;
    let mut incident = vec![Vec::new(); mesh.points.len()];
    let mut masks = vec![0u8; mesh.points.len()];
    for (fi, f) in mesh.faces.iter().enumerate() {
        for &v in f {
            incident[v].push(fi);
            masks[v] |= 1 << mesh.labels[fi];
        }
    }
    for _ in 0..passes {
        let mut moved = 0;
        for i in 0..mesh.points.len() {
            if incident[i].is_empty() {
                continue;
            }
            let old = mesh.points[i];
            let eval = |candidate: Point| {
                let mut loss = 0.0;
                let mut worst = 1.0f64;
                for &fi in &incident[i] {
                    let q = shape(
                        mesh.faces[fi].map(|v| if v == i { candidate } else { mesh.points[v] }),
                        mesh.labels[fi],
                        band,
                    );
                    if q <= 0.0 || !q.is_finite() {
                        return (f64::INFINITY, q);
                    }
                    worst = worst.min(q);
                    loss += -q.ln() + 100.0 * (0.75 - q).max(0.0).powi(2);
                }
                (loss, worst)
            };
            let (before, worst) = eval(old);
            if worst >= 0.76 {
                continue;
            }
            let edge = incident[i]
                .iter()
                .flat_map(|&fi| mesh.faces[fi])
                .filter(|&v| v != i)
                .map(|v| {
                    let d = sub(mesh.points[v], old);
                    dot(d, d).sqrt()
                })
                .sum::<f64>()
                / (2 * incident[i].len()) as f64;
            let h = edge * 1e-4;
            let gradient: Point = std::array::from_fn(|k| {
                let mut a = old;
                let mut b = old;
                a[k] += h;
                b[k] -= h;
                (eval(a).0 - eval(b).0) / (2.0 * h)
            });
            let length = dot(gradient, gradient).sqrt();
            if !length.is_finite() || length < 1e-14 {
                continue;
            }
            for attempt in 0..16 {
                let step = 0.2 * edge / length * 0.5f64.powi(attempt);
                let candidate = std::array::from_fn(|k| old[k] - step * gradient[k]);
                let Ok(candidate) = band.project(candidate, masks[i]) else {
                    continue;
                };
                let value = band.field.value(candidate);
                if masks[i] & 3 == 0 && (value < band.levels[0] || value > band.levels[1]) {
                    continue;
                }
                let (loss, next_worst) = eval(candidate);
                if loss < before - 1e-12 && next_worst >= worst.min(0.65) - 1e-12 {
                    mesh.points[i] = candidate;
                    moved += 1;
                    break;
                }
            }
        }
        if moved == 0 {
            break;
        }
    }
    band.field.check().map_err(MeshingError::GenerationFailed)
}

/// Move periodic copies as one vertex while keeping the extracted cap connectivity.
/// The input must be a periodic field on the requested bounds. This deliberately
/// avoids topology edits, which would need paired operations on opposite caps.
pub fn polish_periodic<F: ScalarField>(
    mesh: &mut TriangleMesh,
    band: &Band<'_, F>,
    periodic: [bool; 3],
    passes: usize,
) -> Result<()> {
    smooth_impl(mesh, band, 0, true)?;
    let mut incident = vec![Vec::new(); mesh.points.len()];
    let mut masks = vec![0u8; mesh.points.len()];
    for (fi, face) in mesh.faces.iter().enumerate() {
        for &v in face {
            incident[v].push(fi);
            masks[v] |= 1 << mesh.labels[fi];
        }
    }
    let mut orbits = std::collections::BTreeMap::<[i64; 3], Vec<usize>>::new();
    for (i, p) in mesh.points.iter().enumerate() {
        let key = std::array::from_fn(|axis| {
            let width = band.bounds[1][axis] - band.bounds[0][axis];
            let mut value = (p[axis] - band.bounds[0][axis]) / width;
            if periodic[axis] && (value - 1.0).abs() < 1e-9 {
                value = 0.0;
            }
            (value * 1e9).round() as i64
        });
        orbits.entry(key).or_default().push(i);
    }
    let orbits: Vec<_> = orbits.into_values().collect();
    if orbits.iter().any(|orbit| orbit.len() > 8) {
        return Err(fail("Coincident periodic surface vertices"));
    }
    let mut owner = vec![0; mesh.points.len()];
    let mut slot = vec![0; mesh.points.len()];
    let mut patches = Vec::with_capacity(orbits.len());
    for (oi, orbit) in orbits.iter().enumerate() {
        let mut patch: Vec<_> = orbit
            .iter()
            .flat_map(|&i| incident[i].iter().copied())
            .collect();
        patch.sort_unstable();
        patch.dedup();
        patches.push(patch);
        for (j, &i) in orbit.iter().enumerate() {
            owner[i] = oi;
            slot[i] = j;
        }
    }
    let mut neighbors = vec![Vec::new(); orbits.len()];
    for face in &mesh.faces {
        let ids = face.map(|i| owner[i]);
        for &a in &ids {
            for &b in &ids {
                if a != b {
                    neighbors[a].push(b);
                }
            }
        }
    }
    for row in &mut neighbors {
        row.sort_unstable();
        row.dedup();
    }
    let mut active: Vec<_> = (0..orbits.len()).collect();
    for _ in 0..passes {
        let mut moved = 0;
        let mut next_active = vec![false; orbits.len()];
        for &oi in &active {
            let orbit = &orbits[oi];
            let base = orbit[0];
            let old = mesh.points[base];
            let patch = &patches[oi];
            let normals: Vec<_> = patch
                .iter()
                .map(|&fi| {
                    let points = mesh.faces[fi].map(|v| mesh.points[v]);
                    band.normal(
                        std::array::from_fn(|axis| {
                            (points[0][axis] + points[1][axis] + points[2][axis]) / 3.0
                        }),
                        mesh.labels[fi],
                    )
                })
                .collect();
            let translated = |candidate: Point| -> [Point; 8] {
                let mut points = [[0.; 3]; 8];
                for (j, &i) in orbit.iter().enumerate() {
                    points[j] = std::array::from_fn(|axis| {
                        let low = band.bounds[0][axis];
                        let high = band.bounds[1][axis];
                        if periodic[axis]
                            && (mesh.points[i][axis] - low).abs() < 1e-9
                            && (old[axis] - high).abs() < 1e-9
                        {
                            candidate[axis] - (high - low)
                        } else if periodic[axis]
                            && (mesh.points[i][axis] - high).abs() < 1e-9
                            && (old[axis] - low).abs() < 1e-9
                        {
                            candidate[axis] + (high - low)
                        } else {
                            candidate[axis]
                        }
                    });
                }
                points
            };
            let eval = |candidate: Point| {
                let replacements = translated(candidate);
                let mut loss = 0.0;
                let mut worst = 1.0f64;
                for (&fi, &normal) in patch.iter().zip(&normals) {
                    let points = mesh.faces[fi].map(|v| {
                        if owner[v] == oi {
                            replacements[slot[v]]
                        } else {
                            mesh.points[v]
                        }
                    });
                    let q = shape_with_normal(points, normal);
                    if q <= 0.0 || !q.is_finite() {
                        return (f64::INFINITY, q);
                    }
                    worst = worst.min(q);
                    loss += -q.ln() + 100.0 * (0.75 - q).max(0.0).powi(2);
                }
                (loss, worst)
            };
            let (before, worst) = eval(old);
            if worst >= 0.76 {
                continue;
            }
            let edge = incident[base]
                .iter()
                .flat_map(|&fi| mesh.faces[fi])
                .filter(|&v| v != base)
                .map(|v| dot(sub(mesh.points[v], old), sub(mesh.points[v], old)).sqrt())
                .sum::<f64>()
                / (2 * incident[base].len()) as f64;
            let target: Point = std::array::from_fn(|axis| {
                incident[base]
                    .iter()
                    .flat_map(|&fi| mesh.faces[fi])
                    .filter(|&v| v != base)
                    .map(|v| mesh.points[v][axis])
                    .sum::<f64>()
                    / (2 * incident[base].len()) as f64
            });
            let mut gradient = None;
            let mut accepted = None;
            for attempt in 0..18 {
                let candidate = if attempt < 2 {
                    let weight = 0.5f64.powi(attempt + 1);
                    std::array::from_fn(|axis| old[axis] + weight * (target[axis] - old[axis]))
                } else {
                    let (direction, length) = gradient.get_or_insert_with(|| {
                        let h = edge * 1e-4;
                        let direction: Point = std::array::from_fn(|axis| {
                            let mut a = old;
                            let mut b = old;
                            a[axis] += h;
                            b[axis] -= h;
                            (eval(a).0 - eval(b).0) / (2.0 * h)
                        });
                        (direction, dot(direction, direction).sqrt())
                    });
                    if !length.is_finite() || *length < 1e-14 {
                        break;
                    }
                    let step = 0.2 * edge / *length * 0.5f64.powi(attempt - 2);
                    std::array::from_fn(|axis| old[axis] - step * direction[axis])
                };
                let Ok(candidate) = band.project(candidate, masks[base]) else {
                    continue;
                };
                let replacements = translated(candidate);
                if orbit.iter().zip(&replacements).any(|(&i, &p)| {
                    let value = band.field.value(p);
                    !value.is_finite()
                        || (masks[i] & 3 == 0 && (value < band.levels[0] || value > band.levels[1]))
                }) {
                    continue;
                }
                let (loss, next_worst) = eval(candidate);
                if loss < before - 1e-12 && next_worst >= worst.min(0.65) - 1e-12 {
                    let exact = translated(candidate);
                    let mut exact_loss = 0.0;
                    let mut exact_worst = 1.0f64;
                    for &fi in patch {
                        let points = mesh.faces[fi].map(|v| {
                            if owner[v] == oi {
                                exact[slot[v]]
                            } else {
                                mesh.points[v]
                            }
                        });
                        let q = shape(points, mesh.labels[fi], band);
                        if q <= 0.0 || !q.is_finite() {
                            exact_loss = f64::INFINITY;
                            break;
                        }
                        exact_worst = exact_worst.min(q);
                        exact_loss += -q.ln() + 100.0 * (0.75 - q).max(0.0).powi(2);
                    }
                    if exact_loss >= before - 1e-12 || exact_worst < worst.min(0.65) - 1e-12 {
                        continue;
                    }
                    accepted = Some(replacements);
                    break;
                }
            }
            if let Some(replacements) = accepted {
                for (&i, &p) in orbit.iter().zip(&replacements) {
                    mesh.points[i] = p;
                }
                moved += 1;
                next_active[oi] = true;
                for &j in &neighbors[oi] {
                    next_active[j] = true;
                }
            }
        }
        if moved == 0 {
            break;
        }
        active = next_active
            .into_iter()
            .enumerate()
            .filter_map(|(i, yes)| yes.then_some(i))
            .collect();
    }
    band.field.check().map_err(MeshingError::GenerationFailed)
}

fn collapse_slivers<F: ScalarField>(mesh: &mut TriangleMesh, band: &Band<'_, F>) -> Result<()> {
    use std::collections::{BTreeMap, BTreeSet};
    let mut edges = BTreeMap::<(usize, usize), Vec<usize>>::new();
    let mut adjacent = vec![BTreeSet::new(); mesh.points.len()];
    let mut incident = vec![Vec::new(); mesh.points.len()];
    let mut masks = vec![0u8; mesh.points.len()];
    let scores: Vec<_> = mesh
        .faces
        .iter()
        .zip(&mesh.labels)
        .map(|(f, &l)| shape(f.map(|i| mesh.points[i]), l, band))
        .collect();
    for (fi, f) in mesh.faces.iter().enumerate() {
        for j in 0..3 {
            let a = f[j];
            let b = f[(j + 1) % 3];
            edges.entry((a.min(b), a.max(b))).or_default().push(fi);
            adjacent[a].insert(b);
            adjacent[b].insert(a);
            incident[a].push(fi);
            masks[a] |= 1 << mesh.labels[fi];
        }
    }
    let length = |a: usize, b: usize| {
        let d = sub(mesh.points[a], mesh.points[b]);
        dot(d, d)
    };
    let mean = edges.keys().map(|&(a, b)| length(a, b)).sum::<f64>() / edges.len() as f64;
    let mut candidates: Vec<_> = edges
        .iter()
        .filter(|((a, b), faces)| {
            length(*a, *b) < 0.36 * mean && faces.iter().any(|&i| scores[i] < 0.4)
        })
        .map(|(&edge, _)| edge)
        .collect();
    candidates.sort_by(|&(a, b), &(c, d)| length(a, b).total_cmp(&length(c, d)));
    let mut touched = vec![false; mesh.points.len()];
    let mut removed = vec![false; mesh.faces.len()];
    for (a, b) in candidates {
        if touched[a] || touched[b] {
            continue;
        }
        let shared = &edges[&(a, b)];
        if shared.len() != 2 {
            continue;
        }
        let opposite: BTreeSet<_> = shared
            .iter()
            .flat_map(|&fi| mesh.faces[fi])
            .filter(|&v| v != a && v != b)
            .collect();
        if adjacent[a]
            .intersection(&adjacent[b])
            .copied()
            .collect::<BTreeSet<_>>()
            != opposite
        {
            continue;
        }
        let mut best = None;
        for (keep, remove) in [(a, b), (b, a)] {
            if masks[keep] | masks[remove] != masks[keep] {
                continue;
            }
            let mut minimum = f64::INFINITY;
            let mut valid = true;
            for &fi in &incident[remove] {
                if shared.contains(&fi) {
                    continue;
                }
                let next = mesh.faces[fi].map(|v| if v == remove { keep } else { v });
                let score = shape(next.map(|v| mesh.points[v]), mesh.labels[fi], band);
                if score < scores[fi].clamp(0.05, 0.3) {
                    valid = false;
                    break;
                }
                minimum = minimum.min(score);
            }
            if valid && best.is_none_or(|(_, _, q)| minimum > q) {
                best = Some((keep, remove, minimum))
            }
        }
        if let Some((keep, remove, _)) = best {
            for &fi in shared {
                removed[fi] = true
            }
            for &fi in &incident[remove] {
                mesh.faces[fi] = mesh.faces[fi].map(|v| if v == remove { keep } else { v })
            }
            for &v in adjacent[a].iter().chain(&adjacent[b]) {
                touched[v] = true
            }
            touched[a] = true;
            touched[b] = true;
        }
    }
    let mut index = 0;
    mesh.faces.retain(|_| {
        let retain = !removed[index];
        index += 1;
        retain
    });
    let mut index = 0;
    mesh.labels.retain(|_| {
        let retain = !removed[index];
        index += 1;
        retain
    });
    let mut mapping = vec![usize::MAX; mesh.points.len()];
    let mut points = Vec::new();
    for f in &mut mesh.faces {
        for v in f {
            if mapping[*v] == usize::MAX {
                mapping[*v] = points.len();
                points.push(mesh.points[*v])
            }
            *v = mapping[*v];
        }
    }
    mesh.points = points;
    Ok(())
}

fn smooth_impl<F: ScalarField>(
    mesh: &mut TriangleMesh,
    band: &Band<'_, F>,
    iterations: usize,
    guard_quality: bool,
) -> Result<()> {
    band.validate()?;
    if mesh.faces.len() != mesh.labels.len()
        || mesh.labels.iter().any(|&label| label > 7)
        || mesh.points.iter().flatten().any(|v| !v.is_finite())
        || mesh.faces.iter().any(|f| {
            f.iter().any(|&v| v >= mesh.points.len())
                || f[0] == f[1]
                || f[1] == f[2]
                || f[0] == f[2]
        })
    {
        return Err(MeshingError::InvalidOptions("Invalid triangle mesh".into()));
    }
    let mut neighbors = vec![Vec::new(); mesh.points.len()];
    let mut incident = vec![Vec::new(); mesh.points.len()];
    let mut masks = vec![0u8; mesh.points.len()];
    for (i, face) in mesh.faces.iter().enumerate() {
        for &a in face {
            incident[a].push(i);
            masks[a] |= 1 << mesh.labels[i];
            for &b in face {
                if a != b {
                    neighbors[a].push(b)
                }
            }
        }
    }
    for adjacent in &mut neighbors {
        adjacent.sort_unstable();
        adjacent.dedup()
    }
    if masks
        .iter()
        .any(|&m| m & 3 == 3 || m & 12 == 12 || m & 48 == 48 || m & 192 == 192)
    {
        return Err(fail("Incompatible vertex constraints"));
    }
    for _ in 0..iterations {
        for i in 0..mesh.points.len() {
            if neighbors[i].is_empty() {
                continue;
            }
            let old = mesh.points[i];
            let floor = if guard_quality {
                incident[i]
                    .iter()
                    .map(|&fi| {
                        shape(
                            mesh.faces[fi].map(|j| mesh.points[j]),
                            mesh.labels[fi],
                            band,
                        )
                    })
                    .fold(0.5, f64::min)
            } else {
                0.0
            };
            let target: Point = std::array::from_fn(|k| {
                neighbors[i].iter().map(|&j| mesh.points[j][k]).sum::<f64>()
                    / neighbors[i].len() as f64
            });
            for step in 0..8 {
                let weight = 0.5f64.powi(step + 1);
                let candidate = std::array::from_fn(|k| old[k] + weight * (target[k] - old[k]));
                let Ok(candidate) = band.project(candidate, masks[i]) else {
                    continue;
                };
                if masks[i] & 3 == 0 {
                    let value = band.field.value(candidate);
                    if !value.is_finite() || value < band.levels[0] || value > band.levels[1] {
                        continue;
                    }
                }
                let valid = incident[i].iter().all(|&fi| {
                    let [a, b, c] =
                        mesh.faces[fi].map(|j| if j == i { candidate } else { mesh.points[j] });
                    let n = band.normal(
                        std::array::from_fn(|k| (a[k] + b[k] + c[k]) / 3.),
                        mesh.labels[fi],
                    );
                    let e = sub(b, a);
                    let d = sub(c, a);
                    dot(cross(e, d), n) > 1e-10 * (dot(e, e) + dot(d, d))
                        && (!guard_quality
                            || shape([a, b, c], mesh.labels[fi], band) >= floor - 1e-12)
                });
                if valid {
                    mesh.points[i] = candidate;
                    break;
                }
            }
        }
    }
    band.field.check().map_err(MeshingError::GenerationFailed)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn check<F: ScalarField>(mesh: &TriangleMesh, band: &Band<'_, F>) {
        let mut edges = HashMap::<(usize, usize), (usize, i32)>::new();
        let mut used = vec![false; mesh.points.len()];
        for (face, &label) in mesh.faces.iter().zip(&mesh.labels) {
            let [a, b, c] = face.map(|i| mesh.points[i]);
            let n = band.normal(std::array::from_fn(|k| (a[k] + b[k] + c[k]) / 3.), label);
            assert!(dot(cross(sub(b, a), sub(c, a)), n) > 0.);
            for j in 0..3 {
                let a = face[j];
                let b = face[(j + 1) % 3];
                used[a] = true;
                let item = edges.entry((a.min(b), a.max(b))).or_default();
                item.0 += 1;
                item.1 += if a < b { 1 } else { -1 };
                let p = mesh.points[a];
                if label < 2 {
                    assert!(
                        (band.field.value(p) - band.levels[if label == 0 { 1 } else { 0 }]).abs()
                            < 1e-9
                    )
                } else {
                    assert!(
                        (p[(label as usize - 2) / 2]
                            - band.bounds[(label as usize - 2) % 2][(label as usize - 2) / 2])
                            .abs()
                            < 1e-12
                    )
                }
            }
        }
        assert!(used.into_iter().all(|b| b));
        assert!(edges.values().all(|&e| e == (2, 0)));
    }
    #[test]
    fn clipped_sphere_band_stays_closed_after_smoothing() {
        let field = |p: Point| dot(p, p);
        let band = Band {
            field: &field,
            bounds: [[-0.83; 3], [0.83; 3]],
            levels: [0.2, 1.0],
        };
        let mut mesh = extract(&band, 12).unwrap();
        check(&mesh, &band);
        let faces = mesh.faces.clone();
        smooth(&mut mesh, &band, 4).unwrap();
        check(&mesh, &band);
        assert_eq!(faces, mesh.faces);
        improve(&mut mesh, &band, 4).unwrap();
        check(&mesh, &band);
        let before = mesh.points.clone();
        let faces = mesh.faces.clone();
        let mut masks = vec![0u8; mesh.points.len()];
        for (f, &l) in mesh.faces.iter().zip(&mesh.labels) {
            for &i in f {
                masks[i] |= 1 << l;
            }
        }
        assert!(redistribute_features(&mut mesh, &band, 10).unwrap() > 0);
        assert_eq!(faces, mesh.faces);
        for (i, &mask) in masks.iter().enumerate() {
            if mask.count_ones() != 2 {
                assert_eq!(before[i], mesh.points[i]);
            }
        }
        check(&mesh, &band);
        polish(&mut mesh, &band, 4).unwrap();
        check(&mesh, &band);
    }
    #[test]
    fn lattice_aligned_plane_band_has_shared_cap_vertices() {
        let field = |p: Point| p[0];
        let band = Band {
            field: &field,
            bounds: [[-1.; 3], [1.; 3]],
            levels: [-0.25, 0.25],
        };
        let mut mesh = extract(&band, 8).unwrap();
        check(&mesh, &band);
        smooth(&mut mesh, &band, 3).unwrap();
        check(&mesh, &band);
        improve(&mut mesh, &band, 4).unwrap();
        check(&mesh, &band);
    }
    #[test]
    fn split_p_cap_triangulations_match_on_opposite_faces() {
        let field = |p: Point| {
            let [x, y, z] = p.map(|v| v * std::f64::consts::TAU);
            1.1 * ((2. * x).sin() * y.cos() * z.sin()
                + (2. * y).sin() * z.cos() * x.sin()
                + (2. * z).sin() * x.cos() * y.sin())
                - 0.2
                    * ((2. * x).cos() * (2. * y).cos()
                        + (2. * y).cos() * (2. * z).cos()
                        + (2. * z).cos() * (2. * x).cos())
                - 0.4 * ((2. * x).cos() + (2. * y).cos() + (2. * z).cos())
        };
        let band = Band {
            field: &field,
            bounds: [[-0.5; 3], [0.5; 3]],
            levels: [-0.25, 0.25],
        };
        let mesh = extract(&band, 24).unwrap();
        check(&mesh, &band);
        for axis in 0..3 {
            let caps = [0, 1].map(|side| {
                mesh.faces
                    .iter()
                    .zip(&mesh.labels)
                    .filter(|&(_, &label)| label == (2 + 2 * axis + side) as u8)
                    .map(|(face, _)| {
                        let mut vertices = face.map(|i| {
                            let mut p = mesh.points[i];
                            p[axis] = 0.;
                            p.map(|x| (x * 1e9).round() as i64)
                        });
                        vertices.sort_unstable();
                        vertices
                    })
                    .collect::<std::collections::BTreeSet<_>>()
            });
            assert!(!caps[0].is_empty());
            assert_eq!(caps[0], caps[1]);
        }
    }
    #[test]
    fn invalid_input_returns_errors() {
        let field = |_: Point| f64::NAN;
        let band = Band {
            field: &field,
            bounds: [[-1.; 3], [1.; 3]],
            levels: [-0.25, 0.25],
        };
        assert!(extract(&band, 8).is_err());
        assert!(extract(&band, 0).is_err());
        let mut mesh = TriangleMesh {
            points: vec![[0.; 3]],
            faces: vec![[0, 1, 2]],
            labels: vec![0],
        };
        assert!(smooth(&mut mesh, &band, 1).is_err());
    }
}
