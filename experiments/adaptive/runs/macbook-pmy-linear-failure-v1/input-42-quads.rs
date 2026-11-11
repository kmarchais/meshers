//! Experimental native quad extraction and constrained optimization.
//!
//! This research API is feature-gated. Uniform dual contouring does not yet
//! provide a feature-aligned patch layout or guarantee resolved topology.
use crate::{MeshingError, Point, cross, dot, implicit::ScalarField, norm, sub};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

type Result<T> = std::result::Result<T, MeshingError>;
fn invalid(message: &str) -> MeshingError {
    MeshingError::InvalidOptions(message.into())
}
fn failed(message: &str) -> MeshingError {
    MeshingError::GenerationFailed(message.into())
}
fn add(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] + b[i])
}
fn mul(a: Point, s: f64) -> Point {
    a.map(|v| v * s)
}
fn unit(a: Point) -> Point {
    mul(a, 1.0 / norm(a).max(1e-300))
}

/// Labels: upper/lower implicit wall (0/1), then x-/x+, y-/y+, z-/z+.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct QuadMesh {
    pub points: Vec<Point>,
    pub quads: Vec<[usize; 4]>,
    pub labels: Vec<u8>,
}

/// An implicit band intersected with an axis-aligned box.
pub struct Band<'a, F: ScalarField> {
    pub field: &'a F,
    pub bounds: [Point; 2],
    pub levels: [f64; 2],
}

impl<F: ScalarField> Band<'_, F> {
    pub(crate) fn validate(&self) -> Result<()> {
        if !((self.levels[0].is_finite() || self.levels[0] == f64::NEG_INFINITY)
            && self.levels[1].is_finite()
            && self.levels[0] < self.levels[1])
            || (0..3).any(|i| {
                !self.bounds[0][i].is_finite()
                    || !self.bounds[1][i].is_finite()
                    || self.bounds[0][i] >= self.bounds[1][i]
            })
        {
            return Err(invalid(
                "Ordered levels with a finite upper bound and finite box bounds are required",
            ));
        }
        if !self.length().is_finite() || self.length() <= 0.0 {
            return Err(invalid("Box extent must have a finite positive length"));
        }
        Ok(())
    }
    fn length(&self) -> f64 {
        norm(sub(self.bounds[1], self.bounds[0]))
    }
    pub(crate) fn gradient(&self, p: Point) -> Point {
        self.field.gradient(p).unwrap_or_else(|| {
            let h = self.length() * 1e-6;
            std::array::from_fn(|i| {
                let mut a = p;
                let mut b = p;
                a[i] += h;
                b[i] -= h;
                (self.field.value(a) - self.field.value(b)) / (2.0 * h)
            })
        })
    }
    pub(crate) fn normal(&self, p: Point, label: u8) -> Point {
        if label < 2 {
            mul(unit(self.gradient(p)), if label == 0 { 1.0 } else { -1.0 })
        } else {
            let mut n = [0.0; 3];
            n[(label as usize - 2) / 2] = if label.is_multiple_of(2) { -1.0 } else { 1.0 };
            n
        }
    }
    fn constraints(&self, mesh: &QuadMesh) -> Result<Vec<u8>> {
        self.validate()?;
        if mesh.points.is_empty()
            || mesh.quads.is_empty()
            || mesh.labels.len() != mesh.quads.len()
            || mesh.points.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "Expected nonempty finite points and one label per quad",
            ));
        }
        let mut masks = vec![0_u8; mesh.points.len()];
        for (q, &label) in mesh.quads.iter().zip(&mesh.labels) {
            if label > 7
                || q.iter().any(|&i| i >= masks.len())
                || (0..4).any(|i| (i + 1..4).any(|j| q[i] == q[j]))
            {
                return Err(invalid("Invalid quad indices or surface label"));
            }
            for &v in q {
                masks[v] |= 1 << label;
            }
        }
        if masks
            .iter()
            .any(|&m| m == 0 || (0..4).any(|j| (m >> (2 * j)) & 3 == 3))
        {
            return Err(failed(
                "Unresolved or incompatible surface constraints; increase extraction resolution",
            ));
        }
        Ok(masks)
    }
    pub(crate) fn project(&self, mut p: Point, mask: u8) -> Result<Point> {
        let free: Point = std::array::from_fn(|i| {
            if (mask >> (2 + 2 * i)) & 3 == 0 {
                1.0
            } else {
                0.0
            }
        });
        for (i, coordinate) in p.iter_mut().enumerate() {
            for side in 0..2 {
                if mask & (1 << (2 + 2 * i + side)) != 0 {
                    *coordinate = self.bounds[side][i];
                }
            }
        }
        if mask & 3 != 0 {
            let target = self.levels[usize::from(mask & 1 != 0)];
            let mut converged = false;
            for _ in 0..20 {
                let residual = self.field.value(p) - target;
                let g = self.gradient(p);
                let g = std::array::from_fn(|i| g[i] * free[i]);
                let gn = norm(g);
                if !residual.is_finite() || !gn.is_finite() || gn <= f64::MIN_POSITIVE {
                    return Err(failed("Singular/non-finite field during projection"));
                }
                let distance = residual / gn;
                if distance.abs() <= self.length() * 1e-12 {
                    converged = true;
                    break;
                }
                p = sub(
                    p,
                    mul(
                        unit(g),
                        distance.clamp(-self.length() * 0.05, self.length() * 0.05),
                    ),
                );
            }
            if !converged {
                return Err(failed("Surface projection did not converge"));
            }
        }
        if (0..3).any(|i| {
            !p[i].is_finite()
                || p[i] < self.bounds[0][i] - self.length() * 1e-10
                || p[i] > self.bounds[1][i] + self.length() * 1e-10
        }) {
            return Err(failed("Projection left the clipping box"));
        }
        Ok(p)
    }
    fn tangent(&self, p: Point, mut g: Point, mask: u8) -> Point {
        for (i, v) in g.iter_mut().enumerate() {
            if (mask >> (2 + 2 * i)) & 3 != 0 {
                *v = 0.0;
            }
        }
        if mask & 3 != 0 {
            let mut n = self.gradient(p);
            for (i, v) in n.iter_mut().enumerate() {
                if (mask >> (2 + 2 * i)) & 3 != 0 {
                    *v = 0.0;
                }
            }
            let n = unit(n);
            g = sub(g, mul(n, dot(g, n)));
        }
        g
    }
}

/// Preserve face connectivity while inserting shared edge midpoints and centers.
pub fn refine(mesh: &QuadMesh) -> Result<QuadMesh> {
    if mesh.quads.len() != mesh.labels.len() {
        return Err(invalid("Missing quad labels"));
    }
    let mut out = QuadMesh {
        points: mesh.points.clone(),
        quads: Vec::new(),
        labels: Vec::new(),
    };
    let mut edges = BTreeMap::new();
    for (q, &label) in mesh.quads.iter().zip(&mesh.labels) {
        if q.iter().any(|&v| v >= mesh.points.len()) {
            return Err(invalid("Invalid quad index"));
        }
        let mut mid = [0; 4];
        for j in 0..4 {
            let a = q[j];
            let b = q[(j + 1) % 4];
            mid[j] = *edges.entry((a.min(b), a.max(b))).or_insert_with(|| {
                let id = out.points.len();
                out.points
                    .push(mul(add(mesh.points[a], mesh.points[b]), 0.5));
                id
            });
        }
        let center = out.points.len();
        out.points.push(mul(
            q.iter().fold([0.0; 3], |s, &v| add(s, mesh.points[v])),
            0.25,
        ));
        for j in 0..4 {
            out.quads.push([q[j], mid[j], center, mid[(j + 3) % 4]]);
            out.labels.push(label);
        }
    }
    Ok(out)
}

/// Catmull-Clark refinement with crease stencils on labelled patch boundaries.
/// All new positions are projected onto their incident analytic constraints.
pub fn subdivide<F: ScalarField>(mesh: &QuadMesh, geometry: &Band<'_, F>) -> Result<QuadMesh> {
    let masks = geometry.constraints(mesh)?;
    let mut out = refine(mesh)?;
    let centers = mesh
        .quads
        .iter()
        .map(|q| {
            mul(
                q.iter().fold([0.0; 3], |s, &v| add(s, mesh.points[v])),
                0.25,
            )
        })
        .collect::<Vec<_>>();
    let mut edges: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    let mut incident = vec![Vec::new(); mesh.points.len()];
    for (face, q) in mesh.quads.iter().enumerate() {
        for j in 0..4 {
            let a = q[j];
            let b = q[(j + 1) % 4];
            edges.entry((a.min(b), a.max(b))).or_default().push(face);
            incident[a].push(face);
        }
    }
    let mut neighbors = vec![Vec::new(); mesh.points.len()];
    let mut crease = neighbors.clone();
    for (&(a, b), faces) in &edges {
        if faces.len() != 2 {
            return Err(failed("Smooth refinement requires closed manifold edges"));
        }
        neighbors[a].push(b);
        neighbors[b].push(a);
        if mesh.labels[faces[0]] != mesh.labels[faces[1]] {
            crease[a].push(b);
            crease[b].push(a);
        }
    }
    for v in 0..mesh.points.len() {
        let p = mesh.points[v];
        out.points[v] = if masks[v].count_ones() > 2 || crease[v].len() > 2 {
            p
        } else if crease[v].len() == 2 {
            mul(
                add(
                    mul(p, 6.0),
                    add(mesh.points[crease[v][0]], mesh.points[crease[v][1]]),
                ),
                0.125,
            )
        } else if crease[v].is_empty() {
            let n = incident[v].len() as f64;
            let face = mul(
                incident[v]
                    .iter()
                    .fold([0.0; 3], |s, &f| add(s, centers[f])),
                1.0 / n,
            );
            let edge = mul(
                neighbors[v]
                    .iter()
                    .fold([0.0; 3], |s, &w| add(s, mul(add(p, mesh.points[w]), 0.5))),
                1.0 / neighbors[v].len() as f64,
            );
            mul(add(add(face, mul(edge, 2.0)), mul(p, n - 3.0)), 1.0 / n)
        } else {
            p
        };
    }
    // refine() keeps each parent quad's four children together, so child j's
    // second vertex identifies the shared midpoint of parent edge j.
    for (face, q) in mesh.quads.iter().enumerate() {
        for j in 0..4 {
            let a = q[j];
            let b = q[(j + 1) % 4];
            let faces = &edges[&(a.min(b), a.max(b))];
            if mesh.labels[faces[0]] == mesh.labels[faces[1]] {
                let midpoint = out.quads[4 * face + j][1];
                out.points[midpoint] = mul(
                    add(
                        add(mesh.points[a], mesh.points[b]),
                        add(centers[faces[0]], centers[faces[1]]),
                    ),
                    0.25,
                );
            }
        }
    }
    let masks = geometry.constraints(&out)?;
    out.points = out
        .points
        .iter()
        .zip(masks)
        .map(|(&p, m)| geometry.project(p, m))
        .collect::<Result<Vec<_>>>()?;
    Ok(out)
}

/// Uniform dual extraction. No triangles or tetrahedra are constructed.
/// Ambiguous cells and unresolved thin bands may require a finer lattice.
pub fn extract_candidate<F: ScalarField>(geometry: &Band<'_, F>, cells: usize) -> Result<QuadMesh> {
    geometry.validate()?;
    if !geometry.levels[0].is_finite() {
        return Err(invalid(
            "Dual quad extraction requires two finite band levels",
        ));
    }
    if !(4..=128).contains(&cells) {
        return Err(invalid("Extraction cells must be in 4..=128"));
    }
    let h: Point =
        std::array::from_fn(|i| (geometry.bounds[1][i] - geometry.bounds[0][i]) / cells as f64);
    let start: Point = std::array::from_fn(|i| geometry.bounds[0][i] - 0.513 * h[i]);
    let n = cells + 2;
    let m = n + 1;
    let point = |c: [usize; 3]| std::array::from_fn(|i| start[i] + c[i] as f64 * h[i]);
    let index = |c: [usize; 3]| (c[0] * m + c[1]) * m + c[2];
    let cell_index = |c: [usize; 3]| (c[0] * n + c[1]) * n + c[2];
    let sample = |p: Point| {
        let v = geometry.field.value(p);
        // Normalize field magnitude by band width, not by a fixed epsilon.
        let s = (geometry.levels[1] - geometry.levels[0]) / geometry.length();
        let mut best = ((v - geometry.levels[1]) / s, 0_u8);
        if (geometry.levels[0] - v) / s > best.0 {
            best = ((geometry.levels[0] - v) / s, 1);
        }
        for (a, &coordinate) in p.iter().enumerate() {
            for side in 0..2 {
                let d = if side == 0 {
                    geometry.bounds[0][a] - coordinate
                } else {
                    coordinate - geometry.bounds[1][a]
                };
                if d > best.0 {
                    best = (d, (2 + 2 * a + side) as u8);
                }
            }
        }
        best
    };
    let mut values = vec![0.0; m * m * m];
    for i in 0..m {
        for j in 0..m {
            for k in 0..m {
                let c = [i, j, k];
                values[index(c)] = sample(point(c)).0;
            }
        }
    }
    if values.iter().any(|v| !v.is_finite()) {
        return Err(failed("Non-finite field samples"));
    }
    let mut mesh = QuadMesh {
        points: Vec::new(),
        quads: Vec::new(),
        labels: Vec::new(),
    };
    let mut lookup = vec![usize::MAX; n * n * n * 12];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let c = [i, j, k];
                let mut intersections = [None; 12];
                let mut endpoints = [([0; 3], [0; 3]); 12];
                for axis in 0..3 {
                    for b in 0..2 {
                        for d in 0..2 {
                            let u = (axis + 1) % 3;
                            let w = (axis + 2) % 3;
                            let mut a = c;
                            a[u] += b;
                            a[w] += d;
                            let mut z = a;
                            z[axis] += 1;
                            let edge = axis * 4 + b * 2 + d;
                            endpoints[edge] = (a, z);
                            let sign = values[index(a)] < 0.0;
                            if sign == (values[index(z)] < 0.0) {
                                continue;
                            }
                            let mut pa = point(a);
                            let mut pb = point(z);
                            for _ in 0..32 {
                                let mid = mul(add(pa, pb), 0.5);
                                if (sample(mid).0 < 0.0) == sign {
                                    pa = mid;
                                } else {
                                    pb = mid;
                                }
                            }
                            intersections[edge] = Some(mul(add(pa, pb), 0.5));
                        }
                    }
                }
                // Split disconnected contour loops inside a cell. Resolve each
                // ambiguous face with the same geometric sample on both sides.
                let mut groups = std::array::from_fn::<_, 12, _>(|e| e);
                for axis in 0..3 {
                    for side in 0..2 {
                        let active = (0..12)
                            .filter(|&e| {
                                intersections[e].is_some()
                                    && endpoints[e].0[axis] == c[axis] + side
                                    && endpoints[e].1[axis] == c[axis] + side
                            })
                            .collect::<Vec<_>>();
                        let mut pairs = Vec::new();
                        if active.len() == 2 {
                            pairs.push((active[0], active[1]));
                        } else if active.len() == 4 {
                            let mut center =
                                std::array::from_fn(|i| start[i] + (c[i] as f64 + 0.5) * h[i]);
                            center[axis] = start[axis] + (c[axis] + side) as f64 * h[axis];
                            let inside = sample(center).0 < 0.0;
                            for (i, &e) in active.iter().enumerate() {
                                for &f in &active[i + 1..] {
                                    for v in [endpoints[e].0, endpoints[e].1] {
                                        if (v == endpoints[f].0 || v == endpoints[f].1)
                                            && (values[index(v)] < 0.0) != inside
                                        {
                                            pairs.push((e, f));
                                        }
                                    }
                                }
                            }
                        }
                        for (a, b) in pairs {
                            let old = groups[b];
                            let new = groups[a];
                            for g in &mut groups {
                                if *g == old {
                                    *g = new;
                                }
                            }
                        }
                    }
                }
                for group in 0..12 {
                    let members = (0..12)
                        .filter(|&e| groups[e] == group && intersections[e].is_some())
                        .collect::<Vec<_>>();
                    if members.is_empty() {
                        continue;
                    }
                    let id = mesh.points.len();
                    let sum = members
                        .iter()
                        .fold([0.0; 3], |s, &e| add(s, intersections[e].unwrap()));
                    mesh.points.push(mul(sum, 1.0 / members.len() as f64));
                    for e in members {
                        lookup[cell_index(c) * 12 + e] = id;
                    }
                }
            }
        }
    }
    for axis in 0..3 {
        let u = (axis + 1) % 3;
        let w = (axis + 2) % 3;
        for a in 0..n {
            for b in 1..n {
                for d in 1..n {
                    let mut c = [0; 3];
                    c[axis] = a;
                    c[u] = b;
                    c[w] = d;
                    let mut end = c;
                    end[axis] += 1;
                    let sign = values[index(c)] < 0.0;
                    if sign == (values[index(end)] < 0.0) {
                        continue;
                    }
                    let mut q = [0; 4];
                    for (j, (du, dw)) in [(1, 1), (0, 1), (0, 0), (1, 0)].into_iter().enumerate() {
                        let mut t = c;
                        t[u] -= du;
                        t[w] -= dw;
                        q[j] = lookup[cell_index(t) * 12 + axis * 4 + du * 2 + dw];
                    }
                    if q.contains(&usize::MAX) {
                        return Err(failed("Incomplete dual cell"));
                    }
                    if !sign {
                        q.reverse();
                    }
                    let mut pa = point(c);
                    let mut pb = point(end);
                    for _ in 0..32 {
                        let mid = mul(add(pa, pb), 0.5);
                        if (sample(mid).0 < 0.0) == sign {
                            pa = mid;
                        } else {
                            pb = mid;
                        }
                    }
                    mesh.quads.push(q);
                    mesh.labels.push(sample(mul(add(pa, pb), 0.5)).1);
                }
            }
        }
    }
    geometry
        .field
        .check()
        .map_err(MeshingError::GenerationFailed)?;
    geometry.constraints(&mesh)?;
    Ok(mesh)
}

/// Local same-patch edge rotations that reduce valence error. This is only a
/// layout experiment: local edits do not replace a global feature-aligned field.
pub fn regularize<F: ScalarField>(
    mesh: &mut QuadMesh,
    geometry: &Band<'_, F>,
    passes: usize,
) -> Result<usize> {
    regularize_impl(mesh, geometry, passes, false)
}

/// Also permit valence-neutral rotations when they improve the pair's worst shape.
pub fn improve_layout<F: ScalarField>(
    mesh: &mut QuadMesh,
    geometry: &Band<'_, F>,
    passes: usize,
) -> Result<usize> {
    regularize_impl(mesh, geometry, passes, true)
}

fn regularize_impl<F: ScalarField>(
    mesh: &mut QuadMesh,
    geometry: &Band<'_, F>,
    passes: usize,
    allow_neutral: bool,
) -> Result<usize> {
    geometry.constraints(mesh)?;
    let mut changed = 0;
    for _ in 0..passes {
        let mut adjacency: BTreeMap<(usize, usize), Vec<(usize, usize)>> = BTreeMap::new();
        let mut valence = vec![0_i32; mesh.points.len()];
        for (i, q) in mesh.quads.iter().enumerate() {
            for j in 0..4 {
                valence[q[j]] += 1;
                let a = q[j];
                let b = q[(j + 1) % 4];
                adjacency
                    .entry((a.min(b), a.max(b)))
                    .or_default()
                    .push((i, j));
            }
        }
        let mut touched = vec![false; mesh.quads.len()];
        let mut step = 0;
        for (&(a, b), faces) in &adjacency {
            if faces.len() != 2 {
                continue;
            }
            let (i, ei) = faces[0];
            let (j, ej) = faces[1];
            if touched[i] || touched[j] || mesh.labels[i] != mesh.labels[j] {
                continue;
            }
            let q = mesh.quads[i];
            let r = mesh.quads[j];
            if q[ei] != r[(ej + 1) % 4] || q[(ei + 1) % 4] != r[ej] {
                continue;
            }
            let hex = [
                q[ei],
                r[(ej + 2) % 4],
                r[(ej + 3) % 4],
                q[(ei + 1) % 4],
                q[(ei + 2) % 4],
                q[(ei + 3) % 4],
            ];
            let mut unique = hex;
            unique.sort_unstable();
            if unique.windows(2).any(|w| w[0] == w[1]) {
                continue;
            }
            let center = mul(
                hex.iter().fold([0.0; 3], |p, &v| add(p, mesh.points[v])),
                1.0 / 6.0,
            );
            let n = geometry.normal(center, mesh.labels[i]);
            let shape = |face: [usize; 4]| {
                (0..4)
                    .map(|k| {
                        let p = mesh.points[face[k]];
                        let e = sub(mesh.points[face[(k + 1) % 4]], p);
                        let b = sub(mesh.points[face[(k + 3) % 4]], p);
                        2.0 * dot(cross(e, b), n) / (dot(e, e) + dot(b, b)).max(1e-30)
                    })
                    .fold(1.0, f64::min)
            };
            let old_shape = shape(q).min(shape(r));
            let mut choice = None;
            for k in 1..=2 {
                let c = hex[k];
                let d = hex[k + 3];
                if adjacency.contains_key(&(c.min(d), c.max(d))) {
                    continue;
                }
                let gain = [a, b]
                    .iter()
                    .map(|&v| (valence[v] - 4).pow(2) - (valence[v] - 5).pow(2))
                    .sum::<i32>()
                    + [c, d]
                        .iter()
                        .map(|&v| (valence[v] - 4).pow(2) - (valence[v] - 3).pow(2))
                        .sum::<i32>();
                if gain < 0 || (gain == 0 && !allow_neutral) {
                    continue;
                }
                let first = std::array::from_fn(|v| hex[(k + v) % 6]);
                let second = [hex[k], hex[(k + 3) % 6], hex[(k + 4) % 6], hex[(k + 5) % 6]];
                let next_shape = shape(first).min(shape(second));
                if gain == 0 && next_shape < old_shape + 0.02 {
                    continue;
                }
                if next_shape < old_shape.min(0.25) {
                    continue;
                }
                if choice
                    .as_ref()
                    .is_none_or(|(best, _, _, _, _)| gain > *best)
                {
                    choice = Some((gain, c, d, first, second));
                }
            }
            if let Some((_, c, d, first, second)) = choice {
                mesh.quads[i] = first;
                mesh.quads[j] = second;
                valence[a] -= 1;
                valence[b] -= 1;
                valence[c] += 1;
                valence[d] += 1;
                touched[i] = true;
                touched[j] = true;
                step += 1;
                // Prevent use of stale adjacency at any neighboring face.
                for &v in &hex {
                    for (&(u, w), neighbors) in &adjacency {
                        if u == v || w == v {
                            for &(face, _) in neighbors {
                                touched[face] = true;
                            }
                        }
                    }
                }
            }
        }
        changed += step;
        if step == 0 {
            break;
        }
    }
    Ok(changed)
}

/// Experimental diagonal quad collapse with local link and geometry checks.
/// Reduces coarse lattice singularities before refinement; validate afterwards.
pub fn coarsen<F: ScalarField>(
    mesh: &mut QuadMesh,
    geometry: &Band<'_, F>,
    passes: usize,
) -> Result<usize> {
    let mut total = 0;
    for _ in 0..passes {
        let masks = geometry.constraints(mesh)?;
        let mut incident = vec![Vec::new(); mesh.points.len()];
        let mut neighbors = vec![std::collections::BTreeSet::new(); mesh.points.len()];
        for (i, q) in mesh.quads.iter().enumerate() {
            for k in 0..4 {
                incident[q[k]].push(i);
                neighbors[q[k]].insert(q[(k + 1) % 4]);
                neighbors[q[k]].insert(q[(k + 3) % 4]);
            }
        }
        let mut used = vec![false; mesh.points.len()];
        let mut removed = vec![false; mesh.quads.len()];
        let mut count = 0;
        for (face, remove) in removed.iter_mut().enumerate() {
            let q = mesh.quads[face];
            if q.iter().any(|&v| used[v]) {
                continue;
            }
            let mut choice = None;
            for k in 0..2 {
                let a = q[k];
                let c = q[k + 2];
                let b = q[(k + 1) % 4];
                let d = q[(k + 3) % 4];
                let merged = masks[a] | masks[c];
                if (0..4).any(|j| (merged >> (2 * j)) & 3 == 3)
                    || (masks[a] != merged && masks[c] != merged)
                {
                    continue;
                }
                if incident[b].len() < 4 || incident[d].len() < 4 {
                    continue;
                }
                let common = neighbors[a]
                    .intersection(&neighbors[c])
                    .copied()
                    .collect::<Vec<_>>();
                if common.len() != 2 || !common.contains(&b) || !common.contains(&d) {
                    continue;
                }
                let val = |v: usize| incident[v].len() as i32;
                let gain = (val(a) - 4).pow(2)
                    + (val(c) - 4).pow(2)
                    + (val(b) - 4).pow(2)
                    + (val(d) - 4).pow(2)
                    - (val(a) + val(c) - 6).pow(2)
                    - (val(b) - 5).pow(2)
                    - (val(d) - 5).pow(2);
                if gain < 0 {
                    continue;
                }
                let Ok(position) =
                    geometry.project(mul(add(mesh.points[a], mesh.points[c]), 0.5), merged)
                else {
                    continue;
                };
                let mut affected = incident[a].clone();
                affected.extend(&incident[c]);
                affected.sort_unstable();
                affected.dedup();
                if affected
                    .iter()
                    .any(|&f| mesh.quads[f].iter().any(|&v| used[v]))
                {
                    continue;
                }
                let mut safe = true;
                for &f in &affected {
                    if f == face {
                        continue;
                    }
                    let old = mesh.quads[f];
                    if old.contains(&a) && old.contains(&c) {
                        safe = false;
                        break;
                    }
                    let normal = geometry.normal(
                        mul(
                            old.iter().fold([0.0; 3], |p, &v| add(p, mesh.points[v])),
                            0.25,
                        ),
                        mesh.labels[f],
                    );
                    let shape = |points: [Point; 4]| {
                        (0..4)
                            .map(|i| {
                                let e = sub(points[(i + 1) % 4], points[i]);
                                let b = sub(points[(i + 3) % 4], points[i]);
                                2.0 * dot(cross(e, b), normal) / (dot(e, e) + dot(b, b)).max(1e-30)
                            })
                            .fold(1.0, f64::min)
                    };
                    let before = shape(old.map(|v| mesh.points[v]));
                    let after = shape(old.map(|v| {
                        if v == a || v == c {
                            position
                        } else {
                            mesh.points[v]
                        }
                    }));
                    if after < before.min(0.25) {
                        safe = false;
                        break;
                    }
                }
                if safe
                    && choice
                        .as_ref()
                        .is_none_or(|(best, _, _, _, _)| gain > *best)
                {
                    choice = Some((gain, a, c, position, affected));
                }
            }
            if let Some((_, a, c, position, affected)) = choice {
                mesh.points[a] = position;
                *remove = true;
                count += 1;
                for f in affected {
                    for v in &mut mesh.quads[f] {
                        used[*v] = true;
                        if *v == c {
                            *v = a;
                        }
                    }
                }
            }
        }
        if count == 0 {
            break;
        }
        let mut faces = Vec::new();
        let mut labels = Vec::new();
        let mut points = Vec::new();
        let mut remap = vec![usize::MAX; mesh.points.len()];
        for (i, &q) in mesh.quads.iter().enumerate() {
            if !removed[i] {
                faces.push(q.map(|v| {
                    if remap[v] == usize::MAX {
                        remap[v] = points.len();
                        points.push(mesh.points[v]);
                    }
                    remap[v]
                }));
                labels.push(mesh.labels[i]);
            }
        }
        *mesh = QuadMesh {
            points,
            quads: faces,
            labels,
        };
        total += count;
    }
    Ok(total)
}

/// Insert one conforming quad row inside each feature-bounded patch. A concave
/// boundary vertex must have an inward edge; one incident patch quad cannot
/// represent a reflex boundary angle with positive corner Jacobians.
pub fn boundary_rows<F: ScalarField>(mesh: &QuadMesh, geometry: &Band<'_, F>) -> Result<QuadMesh> {
    let masks = geometry.constraints(mesh)?;
    let mut out = QuadMesh {
        points: mesh.points.clone(),
        quads: Vec::new(),
        labels: Vec::new(),
    };
    let mut centers: BTreeMap<(usize, u8), (Point, usize)> = BTreeMap::new();
    type DirectedPatchEdge = (usize, usize, u8);
    let mut edges: BTreeMap<(usize, usize), Vec<DirectedPatchEdge>> = BTreeMap::new();
    for (q, &label) in mesh.quads.iter().zip(&mesh.labels) {
        let center = mul(
            q.iter().fold([0.0; 3], |p, &v| add(p, mesh.points[v])),
            0.25,
        );
        for j in 0..4 {
            let a = q[j];
            let b = q[(j + 1) % 4];
            edges
                .entry((a.min(b), a.max(b)))
                .or_default()
                .push((a, b, label));
            if masks[a].count_ones() > 1 {
                let entry = centers.entry((a, label)).or_insert(([0.0; 3], 0));
                entry.0 = add(entry.0, center);
                entry.1 += 1;
            }
        }
    }
    let mut duplicates = BTreeMap::new();
    for ((v, label), (sum, count)) in centers {
        let delta = sub(mul(sum, 1.0 / count as f64), mesh.points[v]);
        let mut step = 0.7;
        let mut candidate = None;
        for _ in 0..20 {
            if let Ok(p) = geometry.project(add(mesh.points[v], mul(delta, step)), 1 << label) {
                candidate = Some(p);
                break;
            }
            step *= 0.5;
        }
        let point = candidate
            .ok_or_else(|| failed("Cannot initialize a boundary row inside the geometry"))?;
        duplicates.insert((v, label), out.points.len());
        out.points.push(point);
    }
    for (q, &label) in mesh.quads.iter().zip(&mesh.labels) {
        out.quads
            .push(q.map(|v| *duplicates.get(&(v, label)).unwrap_or(&v)));
        out.labels.push(label);
    }
    for entries in edges.values() {
        if entries.len() != 2 {
            return Err(failed(
                "Boundary rows require closed manifold edge incidence",
            ));
        }
        if entries[0].2 == entries[1].2 {
            continue;
        }
        for &(a, b, label) in entries {
            out.quads
                .push([a, b, duplicates[&(b, label)], duplicates[&(a, label)]]);
            out.labels.push(label);
        }
    }
    geometry.constraints(&out)?;
    Ok(out)
}

#[derive(Clone, Copy)]
pub struct OptimizeOptions {
    pub iterations: usize,
    pub size_weight: f64,
    pub minimum_shape: f64,
}
impl Default for OptimizeOptions {
    fn default() -> Self {
        Self {
            iterations: 320,
            size_weight: 0.35,
            minimum_shape: 0.0,
        }
    }
}
#[derive(Debug, Serialize)]
pub struct Optimization {
    pub iterations: usize,
    pub initial_energy: f64,
    pub final_energy: f64,
    pub line_search_stalled: bool,
    pub objective_evaluations: usize,
    pub objective_seconds: f64,
    pub projection_seconds: f64,
    pub total_seconds: f64,
}

/// Quality of a closed quad surface. Positive sampled Jacobians are a local
/// check, not a proof of global non-intersection or geometric accuracy.
#[derive(Debug, Serialize)]
pub struct Quality {
    pub quads: usize,
    pub irregular_vertices: usize,
    pub boundary_edges: usize,
    pub nonmanifold_edges: usize,
    pub orientation_errors: usize,
    pub bad_vertex_links: usize,
    pub duplicate_faces: usize,
    pub minimum_shape: f64,
    pub shape_size_p05: f64,
    pub shape_size_mean: f64,
    pub minimum_sampled_jacobian: f64,
    pub invalid_cells: usize,
}

/// Evaluate topology and bilinear cell orientation against the intended surface.
pub fn quality<F: ScalarField>(mesh: &QuadMesh, geometry: &Band<'_, F>) -> Result<Quality> {
    geometry.constraints(mesh)?;
    let mut edges: BTreeMap<(usize, usize), Vec<(usize, bool)>> = BTreeMap::new();
    let mut links = vec![Vec::new(); mesh.points.len()];
    let mut unique = std::collections::BTreeSet::new();
    let mut duplicates = 0;
    let mut shape = Vec::new();
    let mut areas = Vec::new();
    let mut minimum_jacobian = f64::INFINITY;
    let mut invalid_cells = 0;
    for (face, (q, &label)) in mesh.quads.iter().zip(&mesh.labels).enumerate() {
        let mut key = *q;
        key.sort_unstable();
        if !unique.insert(key) {
            duplicates += 1;
        }
        let p = q.map(|v| mesh.points[v]);
        let center = mul(p.iter().fold([0.0; 3], |a, &b| add(a, b)), 0.25);
        let n = geometry.normal(center, label);
        if n.iter().any(|v| !v.is_finite()) || norm(n) < 0.5 {
            return Err(failed("Invalid surface normal during quality evaluation"));
        }
        let mut worst = 1.0_f64;
        let mut cell_jacobian = 1.0_f64;
        for i in 0..4 {
            let a = q[i];
            let b = q[(i + 1) % 4];
            edges
                .entry((a.min(b), a.max(b)))
                .or_default()
                .push((face, a < b));
            links[a].push((q[(i + 3) % 4], b));
            let e = sub(p[(i + 1) % 4], p[i]);
            let b = sub(p[(i + 3) % 4], p[i]);
            let denominator = dot(e, e) + dot(b, b);
            worst = worst.min(if denominator > 0.0 {
                2.0 * dot(cross(e, b), n) / denominator
            } else {
                0.0
            });
        }
        for u in 0..5 {
            for v in 0..5 {
                let u = f64::from(u) / 4.0;
                let v = f64::from(v) / 4.0;
                let du = add(mul(sub(p[1], p[0]), 1.0 - v), mul(sub(p[2], p[3]), v));
                let dv = add(mul(sub(p[3], p[0]), 1.0 - u), mul(sub(p[2], p[1]), u));
                let denominator = norm(du) * norm(dv);
                cell_jacobian = cell_jacobian.min(if denominator > 0.0 {
                    dot(cross(du, dv), n) / denominator
                } else {
                    0.0
                });
            }
        }
        if cell_jacobian <= 1e-8 {
            invalid_cells += 1;
        }
        minimum_jacobian = minimum_jacobian.min(cell_jacobian);
        shape.push(worst);
        areas.push(
            0.5 * (norm(cross(sub(p[1], p[0]), sub(p[2], p[0])))
                + norm(cross(sub(p[2], p[0]), sub(p[3], p[0])))),
        );
    }
    let mut bad_links = 0;
    for link in &links {
        let mut graph: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for &(a, b) in link {
            graph.entry(a).or_default().push(b);
            graph.entry(b).or_default().push(a);
        }
        let mut visited = std::collections::BTreeSet::new();
        let mut stack = graph.keys().next().copied().into_iter().collect::<Vec<_>>();
        while let Some(v) = stack.pop() {
            if visited.insert(v) {
                stack.extend(graph[&v].iter().copied());
            }
        }
        if graph.is_empty()
            || visited.len() != graph.len()
            || graph.values().any(|adj| adj.len() != 2)
        {
            bad_links += 1;
        }
    }
    let mean = areas.iter().sum::<f64>() / areas.len() as f64;
    let mut scores = shape
        .iter()
        .zip(&areas)
        .map(|(&s, &a)| {
            s * if mean > 0.0 && a > 0.0 {
                (a / mean).min(mean / a)
            } else {
                0.0
            }
        })
        .collect::<Vec<_>>();
    scores.sort_by(f64::total_cmp);
    let index = 0.05 * (scores.len() - 1) as f64;
    let lo = index.floor() as usize;
    let hi = index.ceil() as usize;
    Ok(Quality {
        quads: mesh.quads.len(),
        irregular_vertices: links.iter().filter(|v| v.len() != 4).count(),
        boundary_edges: edges.values().filter(|v| v.len() == 1).count(),
        nonmanifold_edges: edges.values().filter(|v| v.len() > 2).count(),
        orientation_errors: edges
            .values()
            .filter(|v| v.len() == 2 && v[0].1 == v[1].1)
            .count(),
        bad_vertex_links: bad_links,
        duplicate_faces: duplicates,
        minimum_shape: shape.into_iter().fold(f64::INFINITY, f64::min),
        shape_size_p05: scores[lo] * (1.0 - index.fract()) + scores[hi] * index.fract(),
        shape_size_mean: scores.iter().sum::<f64>() / scores.len() as f64,
        minimum_sampled_jacobian: minimum_jacobian,
        invalid_cells,
    })
}

impl Quality {
    /// Closed manifold connectivity and positive sampled cell Jacobians.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.boundary_edges == 0
            && self.nonmanifold_edges == 0
            && self.orientation_errors == 0
            && self.bad_vertex_links == 0
            && self.duplicate_faces == 0
            && self.invalid_cells == 0
            && self.minimum_sampled_jacobian.is_finite()
    }
}

// Positive extension of determinant permits untangling. Shape/size objective is
// independently implemented from its algebraic definition, without external code.
fn corner(
    e: Point,
    b: Point,
    n: Point,
    area: f64,
    options: OptimizeOptions,
) -> (f64, Point, Point) {
    let d = dot(cross(e, b), n) / area;
    let t = (dot(e, e) + dot(b, b)) / area;
    let r = (d * d + 4e-6).sqrt();
    let p = if d >= 0.0 {
        (d + r) * 0.5
    } else {
        2e-6 / (r - d)
    };
    let pd = p / r;
    let mut energy = t / (2.0 * p) + options.size_weight * (d + 1.0 / p - 2.0);
    let mut dt = 1.0 / (2.0 * p);
    let mut dd = options.size_weight - (t * 0.5 + options.size_weight) * pd / (p * p);
    let deficit = (options.minimum_shape - 2.0 * p / t.max(1e-30)).max(0.0);
    if deficit > 0.0 {
        energy += 10000.0 * deficit * deficit;
        dt += 40000.0 * deficit * p / t.powi(2);
        dd -= 40000.0 * deficit * pd / t;
    }
    (
        energy,
        add(mul(e, 2.0 * dt / area), mul(cross(b, n), dd / area)),
        add(mul(b, 2.0 * dt / area), mul(cross(n, e), dd / area)),
    )
}
fn energy(
    points: &[Point],
    mesh: &QuadMesh,
    normals: &[Point],
    area: f64,
    options: OptimizeOptions,
    direct_quality: bool,
) -> (f64, Vec<Point>) {
    if direct_quality {
        return quality_energy(points, mesh, normals, area, options);
    }
    let mut g = vec![[0.0; 3]; points.len()];
    let mut value = 0.0;
    let scale = 1.0 / (4 * mesh.quads.len()) as f64;
    for (q, &n) in mesh.quads.iter().zip(normals) {
        for j in 0..4 {
            let a = q[j];
            let next = q[(j + 1) % 4];
            let prev = q[(j + 3) % 4];
            let (v, ge, gb) = corner(
                sub(points[next], points[a]),
                sub(points[prev], points[a]),
                n,
                area,
                options,
            );
            value += v * scale;
            g[next] = add(g[next], mul(ge, scale));
            g[prev] = add(g[prev], mul(gb, scale));
            g[a] = sub(g[a], mul(add(ge, gb), scale));
        }
    }
    (value, g)
}

fn quality_energy(
    points: &[Point],
    mesh: &QuadMesh,
    normals: &[Point],
    area: f64,
    options: OptimizeOptions,
) -> (f64, Vec<Point>) {
    let mut gradient = vec![[0.0; 3]; points.len()];
    let mut loss = 0.0;
    let scale = 1.0 / mesh.quads.len() as f64;
    for (face, &normal) in mesh.quads.iter().zip(normals) {
        let e: [Point; 4] =
            std::array::from_fn(|j| sub(points[face[(j + 1) % 4]], points[face[j]]));
        let b: [Point; 4] =
            std::array::from_fn(|j| sub(points[face[(j + 3) % 4]], points[face[j]]));
        let d: [f64; 4] = std::array::from_fn(|j| dot(cross(e[j], b[j]), normal) / area);
        let t: [f64; 4] = std::array::from_fn(|j| (dot(e[j], e[j]) + dot(b[j], b[j])) / area);
        let shape: [f64; 4] = std::array::from_fn(|j| 2.0 * d[j] / t[j].max(1e-30));
        let minimum = shape.into_iter().fold(f64::INFINITY, f64::min);
        let weights = shape.map(|q| (-64.0 * (q - minimum)).exp());
        let sum = weights.iter().sum::<f64>();
        let soft_shape = minimum - (sum * 0.25).ln() / 64.0;
        let signed_area = d.iter().sum::<f64>() * 0.25;
        let root = (signed_area * signed_area + 4e-12).sqrt();
        let positive = if signed_area >= 0.0 {
            (signed_area + root) * 0.5
        } else {
            2e-12 / (root - signed_area)
        };
        let log = positive.ln();
        let smooth_abs = (log * log + 1e-4).sqrt();
        let relative = (-smooth_abs).exp();
        let relative_derivative = -relative * log / smooth_abs / root;
        let score = soft_shape * relative;
        let deficit = (0.48 - score).max(0.0);
        let importance = 1.0 + 80.0 * deficit;
        loss += (1.0 - score + 40.0 * deficit * deficit) * scale;
        for j in 0..4 {
            let mut dt =
                importance * relative * weights[j] / sum * 2.0 * d[j] / t[j].max(1e-30).powi(2);
            let mut dd = importance
                * (-relative * weights[j] / sum * 2.0 / t[j].max(1e-30)
                    - soft_shape * relative_derivative * 0.25);
            let shape_deficit = (options.minimum_shape - shape[j]).max(0.0);
            loss += 1000.0 * shape_deficit * shape_deficit * scale;
            dt += 4000.0 * shape_deficit * d[j] / t[j].max(1e-30).powi(2);
            dd -= 4000.0 * shape_deficit / t[j].max(1e-30);
            let (base, base_e, base_b) = corner(e[j], b[j], normal, area, options);
            let ge = add(
                add(
                    mul(e[j], 2.0 * dt / area),
                    mul(cross(b[j], normal), dd / area),
                ),
                mul(base_e, 0.0125),
            );
            let gb = add(
                add(
                    mul(b[j], 2.0 * dt / area),
                    mul(cross(normal, e[j]), dd / area),
                ),
                mul(base_b, 0.0125),
            );
            loss += base * 0.0125 * scale;
            let a = face[j];
            let next = face[(j + 1) % 4];
            let prev = face[(j + 3) % 4];
            gradient[next] = add(gradient[next], mul(ge, scale));
            gradient[prev] = add(gradient[prev], mul(gb, scale));
            gradient[a] = sub(gradient[a], mul(add(ge, gb), scale));
        }
    }
    (loss, gradient)
}
fn vdot(a: &[Point], b: &[Point]) -> f64 {
    a.iter().zip(b).map(|(&x, &y)| dot(x, y)).sum()
}
fn axpy(a: &mut [Point], s: f64, b: &[Point]) {
    for (x, &y) in a.iter_mut().zip(b) {
        *x = add(*x, mul(y, s));
    }
}

/// Projected limited-memory BFGS, with an energy-decreasing backtracking search.
/// This improves positions only. It cannot repair a poor quad layout.
pub fn optimize<F: ScalarField>(
    mesh: &mut QuadMesh,
    geometry: &Band<'_, F>,
    options: OptimizeOptions,
) -> Result<Optimization> {
    optimize_impl(mesh, geometry, options, false)
}

/// Optimize a smooth approximation to worst-corner shape times relative face area.
pub fn optimize_quality<F: ScalarField>(
    mesh: &mut QuadMesh,
    geometry: &Band<'_, F>,
    options: OptimizeOptions,
) -> Result<Optimization> {
    optimize_impl(mesh, geometry, options, true)
}

fn optimize_impl<F: ScalarField>(
    mesh: &mut QuadMesh,
    geometry: &Band<'_, F>,
    options: OptimizeOptions,
    direct_quality: bool,
) -> Result<Optimization> {
    let started = std::time::Instant::now();
    let mut projection_seconds = 0.0;
    if !options.size_weight.is_finite()
        || options.size_weight < 0.0
        || !options.minimum_shape.is_finite()
        || !(0.0..1.0).contains(&options.minimum_shape)
    {
        return Err(invalid("Invalid optimizer weights"));
    }
    let masks = geometry.constraints(mesh)?;
    let mut x = mesh
        .points
        .iter()
        .zip(&masks)
        .map(|(&p, &m)| geometry.project(p, m))
        .collect::<Result<Vec<_>>>()?;
    let normals = mesh
        .quads
        .iter()
        .zip(&mesh.labels)
        .map(|(q, &l)| geometry.normal(mul(q.iter().fold([0.0; 3], |s, &v| add(s, x[v])), 0.25), l))
        .collect::<Vec<_>>();
    let area = mesh
        .quads
        .iter()
        .map(|q| norm(cross(sub(x[q[1]], x[q[0]]), sub(x[q[3]], x[q[0]]))))
        .sum::<f64>()
        / mesh.quads.len() as f64;
    if !area.is_finite() || area <= 0.0 {
        return Err(failed("Degenerate input area"));
    }
    let objective_start = std::time::Instant::now();
    let (mut value, mut grad) = energy(&x, mesh, &normals, area, options, direct_quality);
    let mut objective_seconds = objective_start.elapsed().as_secs_f64();
    let mut objective_evaluations = 1;
    if !value.is_finite() || grad.iter().flatten().any(|v| !v.is_finite()) {
        return Err(failed("Non-finite optimization objective"));
    }
    let initial = value;
    let mut history: Vec<(Vec<Point>, Vec<Point>, f64)> = Vec::new();
    let mut stalled = false;
    let mut iterations = 0;
    for i in 0..options.iterations {
        for j in 0..grad.len() {
            grad[j] = geometry.tangent(x[j], grad[j], masks[j]);
        }
        let mut direction = grad.clone();
        let mut alphas = Vec::new();
        for (s, y, rho) in history.iter().rev() {
            let alpha = rho * vdot(s, &direction);
            axpy(&mut direction, -alpha, y);
            alphas.push(alpha);
        }
        let factor = history
            .last()
            .map_or(area * x.len() as f64 * 0.1, |(s, y, _)| {
                vdot(s, y) / vdot(y, y).max(1e-300)
            });
        for v in &mut direction {
            *v = mul(*v, factor);
        }
        for ((s, y, rho), alpha) in history.iter().zip(alphas.iter().rev()) {
            let beta = rho * vdot(y, &direction);
            axpy(&mut direction, alpha - beta, s);
        }
        for j in 0..direction.len() {
            direction[j] = mul(geometry.tangent(x[j], direction[j], masks[j]), -1.0);
        }
        if vdot(&grad, &direction) >= 0.0 {
            history.clear();
            direction = grad.iter().map(|&g| mul(g, -factor.abs())).collect();
        }
        let max_move = direction.iter().map(|&v| norm(v)).fold(0.0, f64::max);
        if max_move < geometry.length() * 1e-13 {
            break;
        }
        let mut step = (area.sqrt() * 0.5 / max_move).min(1.0);
        let mut accepted = None;
        for _ in 0..40 {
            let projection_start = std::time::Instant::now();
            let candidate = x
                .iter()
                .zip(&direction)
                .zip(&masks)
                .map(|((&p, &d), &m)| geometry.project(add(p, mul(d, step)), m))
                .collect::<Result<Vec<_>>>();
            projection_seconds += projection_start.elapsed().as_secs_f64();
            if let Ok(next) = candidate {
                let objective_start = std::time::Instant::now();
                let (v, g) = energy(&next, mesh, &normals, area, options, direct_quality);
                objective_seconds += objective_start.elapsed().as_secs_f64();
                objective_evaluations += 1;
                let displacement = next
                    .iter()
                    .zip(&x)
                    .map(|(&a, &b)| sub(a, b))
                    .collect::<Vec<_>>();
                if v.is_finite()
                    && v < value
                    && v <= value + 1e-4 * vdot(&grad, &displacement).min(0.0)
                {
                    accepted = Some((next, v, g, displacement));
                    break;
                }
            }
            step *= 0.5;
        }
        let Some((next, v, mut g, s)) = accepted else {
            stalled = true;
            break;
        };
        for j in 0..g.len() {
            g[j] = geometry.tangent(next[j], g[j], masks[j]);
        }
        let y = g
            .iter()
            .zip(&grad)
            .map(|(&a, &b)| sub(a, b))
            .collect::<Vec<_>>();
        let sy = vdot(&s, &y);
        if sy > 1e-12 * (vdot(&s, &s) * vdot(&y, &y)).sqrt() {
            if history.len() == 12 {
                history.remove(0);
            }
            history.push((s, y, 1.0 / sy));
        }
        x = next;
        value = v;
        grad = g;
        iterations = i + 1;
    }
    geometry
        .field
        .check()
        .map_err(MeshingError::GenerationFailed)?;
    mesh.points = x;
    Ok(Optimization {
        iterations,
        initial_energy: initial,
        final_energy: value,
        line_search_stalled: stalled,
        objective_evaluations,
        objective_seconds,
        projection_seconds,
        total_seconds: started.elapsed().as_secs_f64(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn quality_energy_derivatives_match_finite_differences() {
        let options = OptimizeOptions {
            minimum_shape: 0.3,
            ..Default::default()
        };
        for height in [0.08, 0.4, 1.1] {
            let mesh = QuadMesh {
                points: vec![
                    [0.0, 0.0, 0.0],
                    [1.2, 0.0, 0.02],
                    [0.9, height, 0.01],
                    [0.05, height, 0.0],
                ],
                quads: vec![[0, 1, 2, 3]],
                labels: vec![7],
            };
            let normals = [[0.0, 0.0, 1.0]];
            let (_, gradient) = quality_energy(&mesh.points, &mesh, &normals, 0.8, options);
            for (i, expected) in gradient.iter().enumerate() {
                for axis in 0..3 {
                    let mut plus = mesh.points.clone();
                    let mut minus = mesh.points.clone();
                    plus[i][axis] += 1e-7;
                    minus[i][axis] -= 1e-7;
                    let finite = (quality_energy(&plus, &mesh, &normals, 0.8, options).0
                        - quality_energy(&minus, &mesh, &normals, 0.8, options).0)
                        / 2e-7;
                    assert!(
                        (finite - expected[axis]).abs() < 1e-4 * expected[axis].abs().max(1.0),
                        "{finite} != {}",
                        expected[axis]
                    );
                }
            }
        }
    }
    #[test]
    fn corner_derivatives_match_finite_differences() {
        let options = OptimizeOptions {
            minimum_shape: 0.3,
            ..Default::default()
        };
        for e in [[0.7, 0.2, 0.03], [0.7, -0.02, 0.0], [-0.7, 0.2, 0.0]] {
            let b = [0.1, 0.9, 0.02];
            let n = [0.0, 0.0, 1.0];
            let (_, ge, gb) = corner(e, b, n, 0.4, options);
            for i in 0..6 {
                let mut ea = e;
                let mut eb = e;
                let mut ba = b;
                let mut bb = b;
                let h = 1e-7;
                if i < 3 {
                    ea[i] += h;
                    eb[i] -= h;
                } else {
                    ba[i - 3] += h;
                    bb[i - 3] -= h;
                }
                let finite = (corner(ea, ba, n, 0.4, options).0
                    - corner(eb, bb, n, 0.4, options).0)
                    / (2.0 * h);
                let analytic = if i < 3 { ge[i] } else { gb[i - 3] };
                assert!(
                    (finite - analytic).abs() < 1e-4 * analytic.abs().max(1.0),
                    "{finite} != {analytic}"
                );
            }
        }
    }
    #[test]
    fn projection_is_invariant_to_field_scale() {
        for scale in [1e-15, 1.0, 1e15] {
            let field = |p: Point| scale * dot(p, p);
            let band = Band {
                field: &field,
                bounds: [[-2.0; 3], [2.0; 3]],
                levels: [scale * 0.25, scale],
            };
            let p = band.project([0.7, 0.7, 0.0], 1).unwrap();
            assert!((dot(p, p) - 1.0).abs() < 1e-10);
        }
    }

    fn cube() -> QuadMesh {
        QuadMesh {
            points: vec![
                [-1.0, -1.0, -1.0],
                [1.0, -1.0, -1.0],
                [1.0, 1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
                [1.0, -1.0, 1.0],
                [1.0, 1.0, 1.0],
                [-1.0, 1.0, 1.0],
            ],
            quads: vec![
                [0, 3, 2, 1],
                [4, 5, 6, 7],
                [0, 1, 5, 4],
                [3, 7, 6, 2],
                [0, 4, 7, 3],
                [1, 2, 6, 5],
            ],
            labels: vec![6, 7, 4, 5, 2, 3],
        }
    }
    #[test]
    fn quality_rejects_reversed_and_duplicate_faces() {
        let f = |_: Point| 0.0;
        let geometry = Band {
            field: &f,
            bounds: [[-1.0; 3], [1.0; 3]],
            levels: [-1.0, 1.0],
        };
        let mut mesh = cube();
        let q = quality(&mesh, &geometry).unwrap();
        assert!(q.is_valid());
        assert!((q.shape_size_mean - 1.0).abs() < 1e-12);
        mesh.quads[0].reverse();
        let q = quality(&mesh, &geometry).unwrap();
        assert!(!q.is_valid());
        assert_eq!(q.invalid_cells, 1);
        assert_eq!(q.orientation_errors, 4);
        mesh = cube();
        mesh.quads.push(mesh.quads[0]);
        mesh.labels.push(mesh.labels[0]);
        let q = quality(&mesh, &geometry).unwrap();
        assert!(!q.is_valid());
        assert_eq!(q.duplicate_faces, 1);
        assert_eq!(q.nonmanifold_edges, 4);
    }
    #[test]
    fn refinement_preserves_shared_edges_and_constraints() {
        let f = |_: Point| 0.0;
        let geometry = Band {
            field: &f,
            bounds: [[-1.0; 3], [1.0; 3]],
            levels: [-1.0, 1.0],
        };
        let mesh = refine(&cube()).unwrap();
        assert_eq!(mesh.points.len(), 26);
        assert_eq!(mesh.quads.len(), 24);
        assert!(quality(&mesh, &geometry).unwrap().is_valid());
        let masks = geometry.constraints(&mesh).unwrap();
        for (&p, &m) in mesh.points.iter().zip(&masks) {
            assert_eq!(geometry.project(p, m).unwrap(), p);
        }
    }
    #[test]
    fn optimizer_improves_a_distorted_cube_without_moving_creases() {
        let f = |_: Point| 0.0;
        let geometry = Band {
            field: &f,
            bounds: [[-1.0; 3], [1.0; 3]],
            levels: [-1.0, 1.0],
        };
        let mut mesh = refine(&cube()).unwrap();
        mesh.points[12][0] += 0.2;
        mesh.points[12][1] += 0.1;
        let before = quality(&mesh, &geometry).unwrap().shape_size_mean;
        let result = optimize(&mut mesh, &geometry, OptimizeOptions::default()).unwrap();
        assert!(result.final_energy < result.initial_energy);
        let after = quality(&mesh, &geometry).unwrap();
        assert!(after.is_valid());
        assert!(after.shape_size_mean > before);
        for (actual, expected) in mesh.points[..8].iter().zip(cube().points) {
            assert_eq!(*actual, expected);
        }
    }
    #[test]
    fn crease_subdivision_and_quality_optimization_preserve_cube() {
        let field = |_: Point| 0.0;
        let geometry = Band {
            field: &field,
            bounds: [[-1.0; 3], [1.0; 3]],
            levels: [-1.0, 1.0],
        };
        let mut mesh = subdivide(&cube(), &geometry).unwrap();
        assert_eq!(mesh.points.len(), 26);
        assert_eq!(mesh.quads.len(), 24);
        assert_eq!(&mesh.points[..8], &cube().points);
        assert!(quality(&mesh, &geometry).unwrap().is_valid());
        mesh.points[12][0] += 0.2;
        mesh.points[12][1] += 0.1;
        let before = quality(&mesh, &geometry).unwrap().shape_size_mean;
        optimize_quality(&mut mesh, &geometry, OptimizeOptions::default()).unwrap();
        let after = quality(&mesh, &geometry).unwrap();
        assert!(after.is_valid() && after.shape_size_mean > before);
        assert_eq!(&mesh.points[..8], &cube().points);
        for (&point, mask) in mesh.points.iter().zip(geometry.constraints(&mesh).unwrap()) {
            assert_eq!(geometry.project(point, mask).unwrap(), point);
        }
    }
    #[test]
    fn native_sphere_extraction_is_closed_and_valid() {
        let f = |p: Point| dot(p, p);
        let geometry = Band {
            field: &f,
            bounds: [[-1.2; 3], [1.2; 3]],
            levels: [-1.0, 1.0],
        };
        let mut mesh = extract_candidate(&geometry, 16).unwrap();
        optimize(
            &mut mesh,
            &geometry,
            OptimizeOptions {
                iterations: 60,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(quality(&mesh, &geometry).unwrap().is_valid());
        coarsen(&mut mesh, &geometry, 4).unwrap();
        regularize(&mut mesh, &geometry, 4).unwrap();
        improve_layout(&mut mesh, &geometry, 4).unwrap();
        mesh = subdivide(&mesh, &geometry).unwrap();
        assert!(quality(&mesh, &geometry).unwrap().is_valid());
        for p in mesh.points {
            assert!((dot(p, p) - 1.0).abs() < 1e-9);
        }
    }
    #[test]
    fn invalid_input_is_rejected_without_mutation() {
        let f = |_: Point| 0.0;
        let geometry = Band {
            field: &f,
            bounds: [[-1.0; 3], [1.0; 3]],
            levels: [-1.0, 1.0],
        };
        let mut mesh = cube();
        mesh.quads[0][0] = mesh.points.len();
        let original = mesh.points.clone();
        assert!(optimize(&mut mesh, &geometry, OptimizeOptions::default()).is_err());
        assert_eq!(mesh.points, original);
        assert!(extract_candidate(&geometry, 0).is_err());
    }
    #[test]
    fn quality_detects_a_vertex_pinch_with_two_faces_per_edge() {
        let f = |_: Point| 0.0;
        let geometry = Band {
            field: &f,
            bounds: [[-1.0; 3], [1.0; 3]],
            levels: [-1.0, 1.0],
        };
        let mut mesh = cube();
        let other = cube();
        mesh.points.extend_from_slice(&other.points[1..]);
        for q in other.quads {
            mesh.quads.push(q.map(|v| if v == 0 { 0 } else { v + 7 }));
        }
        mesh.labels.extend(other.labels);
        let report = quality(&mesh, &geometry).unwrap();
        assert_eq!(report.nonmanifold_edges, 0);
        assert_eq!(report.bad_vertex_links, 1);
        assert!(!report.is_valid());
    }
    #[test]
    fn feature_rows_share_vertices_and_preserve_closed_topology() {
        let f = |_: Point| 0.0;
        let geometry = Band {
            field: &f,
            bounds: [[-1.0; 3], [1.0; 3]],
            levels: [-1.0, 1.0],
        };
        let mesh = boundary_rows(&cube(), &geometry).unwrap();
        assert_eq!(mesh.quads.len(), 30);
        assert_eq!(mesh.points.len(), 32);
        assert!(quality(&mesh, &geometry).unwrap().is_valid());
    }
}
