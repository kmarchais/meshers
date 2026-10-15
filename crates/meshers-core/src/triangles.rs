//! Experimental direct triangle surface extraction on a Freudenthal sampling grid.
//! Only crossing edges and boundary polygons are retained. No volume mesh is built.
use crate::{MeshingError, Point, cross, dot, implicit::ScalarField, quads::Band, sub};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct TriangleMesh {
    pub points: Vec<Point>,
    pub faces: Vec<[usize; 3]>,
    pub labels: Vec<u8>,
}

/// Check oriented manifold topology and triangle quality. Intersections are not certified.
pub fn validate(mesh: &TriangleMesh) -> Result<(f64, f64)> {
    if mesh.faces.len() != mesh.labels.len() || mesh.points.iter().flatten().any(|v| !v.is_finite())
    {
        return Err(fail("Invalid surface arrays"));
    }
    if mesh.faces.is_empty() {
        return Ok((0., 0.));
    }
    let mut minimum = 1f64;
    let mut volume = 0.;
    let origin = mesh.points[0];
    let n = mesh.points.len();
    let mut directed = Vec::with_capacity(3 * mesh.faces.len());
    for face in &mesh.faces {
        if face.iter().any(|&i| i >= n) {
            return Err(fail("Invalid triangle index"));
        }
        let [a, b, c] = face.map(|i| sub(mesh.points[i], origin));
        let normal = cross(sub(b, a), sub(c, a));
        let sum = dot(sub(b, a), sub(b, a)) + dot(sub(c, b), sub(c, b)) + dot(sub(a, c), sub(a, c));
        let q = 2. * 3f64.sqrt() * dot(normal, normal).sqrt() / sum;
        if !q.is_finite() || q <= 0. {
            return Err(fail("Degenerate surface triangle"));
        }
        minimum = minimum.min(q);
        volume += dot(a, cross(b, c)) / 6.;
        for i in 0..3 {
            directed.push([face[i], face[(i + 1) % 3]]);
        }
    }
    if !closed_oriented_manifold(&mesh.faces, directed, n) {
        return Err(fail("Surface is not a closed oriented two-manifold"));
    }
    if !volume.is_finite() || volume <= 0. {
        return Err(fail("Nonpositive surface volume"));
    }
    Ok((minimum, volume))
}
/// Every directed edge appears once and its reverse once, and the faces around each
/// vertex form one cycle. Sorting replaces hashing: this runs on every result.
fn closed_oriented_manifold(faces: &[[usize; 3]], mut directed: Vec<[usize; 2]>, n: usize) -> bool {
    directed.sort_unstable();
    if directed.windows(2).any(|w| w[0] == w[1])
        || directed
            .iter()
            .any(|&[a, b]| directed.binary_search(&[b, a]).is_err())
    {
        return false;
    }
    // Link of each vertex as directed edges (next, previous corner), grouped by vertex.
    let mut start = vec![0usize; n + 1];
    for face in faces {
        for &v in face {
            start[v + 1] += 1;
        }
    }
    for v in 0..n {
        start[v + 1] += start[v];
    }
    let mut fill = start.clone();
    let mut link = vec![[0usize; 2]; start[n]];
    for &[a, b, c] in faces {
        for (v, edge) in [(a, [b, c]), (b, [c, a]), (c, [a, b])] {
            link[fill[v]] = edge;
            fill[v] += 1;
        }
    }
    (0..n).all(|v| {
        let ring = &link[start[v]..start[v + 1]];
        if ring.is_empty() {
            return true;
        }
        // Walk the cycle from the first link edge; it must visit every edge.
        let mut current = ring[0][1];
        for steps in 1..=ring.len() {
            match ring.iter().find(|e| e[0] == current) {
                Some(e) if e == &ring[0] => return steps == ring.len(),
                Some(e) => current = e[1],
                None => return false,
            }
        }
        false
    })
}

type Result<T> = std::result::Result<T, MeshingError>;
/// Physical distance estimate at vertices, mid-edges and facet interiors.
pub fn sampled_error<F: ScalarField>(mesh: &TriangleMesh, band: &Band<'_, F>) -> Result<f64> {
    use rayon::prelude::*;
    band.field.check().map_err(MeshingError::GenerationFailed)?;
    // NaN marks a nonfinite distance; `f64::max` would drop it, so reduce explicitly.
    let worst = |a: f64, b: f64| {
        if a.is_nan() || b.is_nan() {
            f64::NAN
        } else {
            a.max(b)
        }
    };
    let maximum = mesh
        .faces
        .par_iter()
        .zip(&mesh.labels)
        .filter(|&(_, &label)| label < 2)
        .map(|(face, &label)| {
            let p = face.map(|i| mesh.points[i]);
            let level = band.levels[if label == 0 { 1 } else { 0 }];
            let mut maximum = 0f64;
            for sample in [
                p[0],
                p[1],
                p[2],
                std::array::from_fn(|a| (p[0][a] + p[1][a] + p[2][a]) / 3.),
                std::array::from_fn(|a| (p[0][a] + p[1][a]) * 0.5),
                std::array::from_fn(|a| (p[1][a] + p[2][a]) * 0.5),
                std::array::from_fn(|a| (p[2][a] + p[0][a]) * 0.5),
            ] {
                let residual = (band.field.value(sample) - level).abs();
                let gradient = crate::norm(band.gradient(sample));
                let distance = if residual == 0. {
                    0.
                } else {
                    residual / gradient
                };
                if !distance.is_finite() {
                    return f64::NAN;
                }
                maximum = maximum.max(distance);
            }
            maximum
        })
        .reduce(|| 0., worst);
    band.field.check().map_err(MeshingError::GenerationFailed)?;
    if maximum.is_nan() {
        return Err(fail("Nonfinite sampled surface distance"));
    }
    Ok(maximum)
}
fn fail(message: &str) -> MeshingError {
    MeshingError::GenerationFailed(message.into())
}
/// Freudenthal lattice edge directions: three axes, three face diagonals and the body
/// diagonal, each from the lower grid node.
const STEPS: [[usize; 3]; 7] = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 1, 0],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
];

enum Root {
    Node(usize),
    Point(Point),
}

/// Level crossings on the sampling lattice. Every crossing is solved once, in
/// parallel, before triangulation, so slabs of cubes can be triangulated
/// independently and the vertex numbering does not depend on scheduling.
struct Lattice<'a, 'b, F: ScalarField> {
    band: &'a Band<'b, F>,
    n: usize,
    offsets: [usize; 7],
    grid: Vec<Point>,
    values: Vec<f64>,
    /// Crossing lattice edges `(lower node, direction, label)`, sorted.
    keys: Vec<(usize, u8, u8)>,
    key_vertex: Vec<usize>,
    node_vertex: Vec<usize>,
    points: Vec<Point>,
}

type Faces = Vec<([usize; 3], u8)>;

impl<F: ScalarField> Lattice<'_, '_, F> {
    fn level(&self, label: u8) -> f64 {
        self.band.levels[if label == 0 { 1 } else { 0 }]
    }

    fn solve(&self, a: usize, b: usize, label: u8) -> Result<Root> {
        let level = self.level(label);
        let va = self.values[a] - level;
        let vb = self.values[b] - level;
        if va == 0.0 {
            return Ok(Root::Node(a));
        }
        if vb == 0.0 {
            return Ok(Root::Node(b));
        }
        if va.signum() == vb.signum() {
            return Err(fail("No bracket on surface edge"));
        }
        let mut lo = 0.;
        let mut hi = 1.;
        let mut flo = va;
        let mut fhi = vb;
        let value_scale = va
            .abs()
            .max(vb.abs())
            .max(level.abs())
            .max(f64::MIN_POSITIVE);
        let tolerance = 64. * f64::EPSILON * value_scale;
        let mut root = None;
        // Secant interpolation converges quickly on resolved smooth edges.
        // Periodic bisection prevents an endpoint from stalling. A bounded
        // bisection fallback retains the original bracket-resolution guarantee.
        for iteration in 0..52 {
            let secant = lo - flo * (hi - lo) / (fhi - flo);
            let t = if iteration >= 12
                || iteration % 4 == 3
                || !secant.is_finite()
                || secant <= lo
                || secant >= hi
            {
                (lo + hi) * 0.5
            } else {
                secant
            };
            let p = std::array::from_fn(|k| self.grid[a][k] * (1. - t) + self.grid[b][k] * t);
            let v = self.band.field.value(p) - level;
            if !v.is_finite() {
                return Err(fail("Non-finite field during edge solve"));
            }
            if v.abs() <= tolerance {
                root = Some(t);
                break;
            }
            if v.signum() == va.signum() {
                lo = t;
                flo = v;
            } else {
                hi = t;
                fhi = v;
            }
        }
        let t = root.unwrap_or((lo + hi) * 0.5);
        Ok(Root::Point(std::array::from_fn(|k| {
            self.grid[a][k] * (1. - t) + self.grid[b][k] * t
        })))
    }

    fn node(&self, id: usize) -> Result<usize> {
        match self.node_vertex[id] {
            usize::MAX => Err(fail("Unassigned lattice vertex")),
            v => Ok(v),
        }
    }

    fn vertex(&self, a: usize, b: usize, label: u8) -> Result<usize> {
        let (a, b) = if a < b { (a, b) } else { (b, a) };
        let direction = self
            .offsets
            .iter()
            .position(|&o| o == b - a)
            .ok_or_else(|| fail("Invalid lattice edge"))?;
        if let Ok(i) = self.keys.binary_search(&(a, direction as u8, label)) {
            return Ok(self.key_vertex[i]);
        }
        // Clipping treats a node exactly on the level as inside; it is its own root.
        let level = self.level(label);
        if self.values[a] == level {
            return self.node(a);
        }
        if self.values[b] == level {
            return self.node(b);
        }
        Err(fail("No bracket on surface edge"))
    }

    fn polygon(
        &self,
        mut ids: Vec<usize>,
        label: u8,
        normal: Option<Point>,
        out: &mut Faces,
    ) -> Result<()> {
        ids.sort_unstable();
        ids.dedup();
        if ids.len() < 3 {
            return Ok(());
        }
        let center = std::array::from_fn(|k| {
            ids.iter().map(|&i| self.points[i][k]).sum::<f64>() / ids.len() as f64
        });
        let normal = normal.unwrap_or_else(|| self.band.normal(center, label));
        let axis = sub(self.points[ids[0]], center);
        let perpendicular = cross(normal, axis);
        ids.sort_by(|&a, &b| {
            let a = sub(self.points[a], center);
            let b = sub(self.points[b], center);
            dot(a, perpendicular)
                .atan2(dot(a, axis))
                .total_cmp(&dot(b, perpendicular).atan2(dot(b, axis)))
        });
        // Choose the fan with the best worst triangle, avoiding a fixed diagonal bias.
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
                let [a, b, c] = face.map(|i| self.points[i]);
                // A polygon spans one lattice tetrahedron, so its centre normal serves
                // every candidate triangle.
                let e = sub(b, a);
                let d = sub(c, a);
                let edge = sub(c, b);
                let score = dot(cross(e, d), normal) / (dot(e, e) + dot(d, d) + dot(edge, edge));
                worst = worst.min(score);
                faces.push(face);
            }
            if best.as_ref().is_none_or(|(score, _)| worst > *score) {
                best = Some((worst, faces))
            }
        }
        let (score, faces) = best.unwrap();
        if score <= 0. {
            return Err(fail("Inverted triangle after analytic edge intersection"));
        }
        out.extend(faces.into_iter().map(|face| (face, label)));
        Ok(())
    }

    fn cap(&self, face: [usize; 3], label: u8, out: &mut Faces) -> Result<()> {
        // Each clipped point retains its original grid-edge support so both
        // implicit walls and clipping caps share exactly the same vertex IDs.
        let mut poly = face.map(|id| ([id, id], self.values[id], None)).to_vec();
        for wall in [1u8, 0] {
            let level = self.level(wall);
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
                    let id = self.vertex(support[0], support[1], wall)?;
                    next.push(([support[0], support[1]], level, Some(id)));
                }
            }
            poly = next;
        }
        let ids = poly
            .into_iter()
            .map(|(support, _, id)| id.map_or_else(|| self.node(support[0]), Ok))
            .collect::<Result<Vec<_>>>()?;
        self.polygon(ids, label, None, out)
    }

    /// Triangles of one slab of cubes.
    fn slab(&self, z: usize) -> Result<Faces> {
        let (n, cells) = (self.n, self.n - 1);
        let stride = [1, n, n * n];
        let h: Point = std::array::from_fn(|k| {
            (self.band.bounds[1][k] - self.band.bounds[0][k]) / cells as f64
        });
        let mut out = Vec::new();
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
                        let level = self.level(label);
                        if !level.is_finite() {
                            continue;
                        }
                        let mut ids = Vec::new();
                        for i in 0..4 {
                            for j in i + 1..4 {
                                if (self.values[tet[i]] < level) != (self.values[tet[j]] < level) {
                                    ids.push(self.vertex(tet[i], tet[j], label)?)
                                }
                            }
                        }
                        if ids.is_empty() {
                            continue;
                        }
                        // The tetrahedron's edges step along `perm` one axis at a time,
                        // so the linear interpolant's gradient needs no field evaluation.
                        let mut g = [0.; 3];
                        for (k, &axis) in perm.iter().enumerate() {
                            g[axis] = (self.values[tet[k + 1]] - self.values[tet[k]]) / h[axis];
                        }
                        let length = crate::norm(g);
                        let sign = if label == 0 { 1. } else { -1. };
                        let normal = (length > 0.).then(|| g.map(|v| sign * v / length));
                        self.polygon(ids, label, normal, &mut out)?;
                    }
                    // Only cubes touching the box can hold cap faces.
                    if ![x, y, z].iter().any(|&c| c == 0 || c == cells - 1) {
                        continue;
                    }
                    for omit in 0..4 {
                        let mut face = [0; 3];
                        for (slot, j) in (0..4).filter(|&j| j != omit).enumerate() {
                            face[slot] = tet[j];
                        }
                        for label in 2..8 {
                            let axis = (label - 2) / 2;
                            let boundary_index = if (label - 2) % 2 == 0 { 0 } else { cells };
                            if face
                                .iter()
                                .all(|&i| (i / stride[axis]) % n == boundary_index)
                            {
                                self.cap(face, label as u8, &mut out)?;
                            }
                        }
                    }
                }
            }
        }
        Ok(out)
    }
}

pub fn extract<F: ScalarField>(band: &Band<'_, F>, cells: usize) -> Result<TriangleMesh> {
    use rayon::prelude::*;
    band.validate()?;
    if !(4..=128).contains(&cells) {
        return Err(MeshingError::InvalidOptions(
            "Surface cells must be in 4..=128".into(),
        ));
    }
    let n = cells + 1;
    let index = |i: usize| [i % n, (i / n) % n, i / (n * n)];
    let grid: Vec<Point> = (0..n * n * n)
        .map(|i| {
            let c = index(i);
            std::array::from_fn(|k| {
                band.bounds[0][k]
                    + (band.bounds[1][k] - band.bounds[0][k]) * c[k] as f64 / cells as f64
            })
        })
        .collect();
    let mut values: Vec<f64> = if band.field.prefers_batches() {
        band.field.values(&grid)
    } else {
        grid.par_iter().map(|&p| band.field.value(p)).collect()
    };
    band.field.check().map_err(MeshingError::GenerationFailed)?;
    if values.iter().any(|v| !v.is_finite()) {
        return Err(fail("Non-finite sampled field"));
    }
    // A zero iso-value has no relative scale of its own. Use the sampled
    // field range so near-zero cancellation ties share one lattice vertex.
    // The tolerance scales with the field, preserving positive rescaling.
    let value_scale = values
        .iter()
        .copied()
        .map(f64::abs)
        .chain(
            band.levels
                .iter()
                .copied()
                .filter(|v| v.is_finite())
                .map(f64::abs),
        )
        .fold(f64::MIN_POSITIVE, f64::max);
    for value in &mut values {
        for &level in &band.levels {
            if level.is_finite() && (*value - level).abs() <= 64. * f64::EPSILON * value_scale {
                *value = level;
            }
        }
    }
    let offsets = STEPS.map(|s| s[0] + n * s[1] + n * n * s[2]);
    let mut lattice = Lattice {
        band,
        n,
        offsets,
        grid,
        values,
        keys: Vec::new(),
        key_vertex: Vec::new(),
        node_vertex: Vec::new(),
        points: Vec::new(),
    };
    let level = |label: u8| band.levels[if label == 0 { 1 } else { 0 }];
    lattice.keys = (0..n * n * n)
        .into_par_iter()
        .flat_map_iter(|a| {
            let c = index(a);
            let values = &lattice.values;
            STEPS
                .iter()
                .enumerate()
                .filter(move |(_, s)| (0..3).all(|k| c[k] + s[k] < n))
                .flat_map(move |(d, _)| {
                    let b = a + offsets[d];
                    (0..2u8).filter_map(move |label| {
                        let level = level(label);
                        (level.is_finite() && (values[a] < level) != (values[b] < level))
                            .then_some((a, d as u8, label))
                    })
                })
        })
        .collect();
    let roots = lattice
        .keys
        .par_iter()
        .map(|&(a, d, label)| lattice.solve(a, a + offsets[d as usize], label))
        .collect::<Result<Vec<_>>>()?;
    band.field.check().map_err(MeshingError::GenerationFailed)?;
    // Nodes become vertices when a cap keeps them or they lie exactly on a level.
    // Unused ones are dropped after triangulation.
    lattice.node_vertex = vec![usize::MAX; n * n * n];
    for id in 0..n * n * n {
        let v = lattice.values[id];
        let on_box = index(id).iter().any(|&c| c == 0 || c == cells);
        let inside = v >= band.levels[0] && v <= band.levels[1];
        if (on_box && inside) || band.levels.iter().any(|&l| l.is_finite() && v == l) {
            lattice.node_vertex[id] = lattice.points.len();
            lattice.points.push(lattice.grid[id]);
        }
    }
    lattice.key_vertex = roots
        .into_iter()
        .map(|root| match root {
            Root::Node(id) => lattice.node_vertex[id],
            Root::Point(p) => {
                lattice.points.push(p);
                lattice.points.len() - 1
            }
        })
        .collect();
    let slabs = (0..cells)
        .into_par_iter()
        .map(|z| lattice.slab(z))
        .collect::<Result<Vec<_>>>()?;
    // Compact to the referenced vertices, keeping their order.
    let mut map = vec![usize::MAX; lattice.points.len()];
    let mut mesh = TriangleMesh {
        points: Vec::new(),
        faces: Vec::new(),
        labels: Vec::new(),
    };
    for (face, label) in slabs.into_iter().flatten() {
        for &v in &face {
            if map[v] == usize::MAX {
                map[v] = 0;
            }
        }
        mesh.faces.push(face);
        mesh.labels.push(label);
    }
    for (v, m) in map.iter_mut().enumerate() {
        if *m == 0 {
            *m = mesh.points.len();
            mesh.points.push(lattice.points[v]);
        }
    }
    for face in &mut mesh.faces {
        *face = face.map(|v| map[v]);
    }
    // Avoid a silent success for unresolved or empty geometry in the benchmark.
    if mesh.faces.is_empty() {
        return Err(fail("No resolved triangle surface"));
    }
    band.field.check().map_err(MeshingError::GenerationFailed)?;
    Ok(mesh)
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
    let e = sub(b, a);
    let d = sub(c, a);
    let other = sub(c, b);
    let normal = band.normal(std::array::from_fn(|k| (a[k] + b[k] + c[k]) / 3.0), label);
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
    use std::collections::HashMap;

    #[test]
    fn manifold_check_rejects_open_and_pinched_surfaces() {
        let tetra = |o: usize| {
            [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]].map(|f: [usize; 3]| f.map(|v| v + o))
        };
        let check = |faces: &[[usize; 3]], n: usize| {
            let directed = faces
                .iter()
                .flat_map(|f| [[f[0], f[1]], [f[1], f[2]], [f[2], f[0]]])
                .collect();
            closed_oriented_manifold(faces, directed, n)
        };
        assert!(check(&tetra(0), 4));
        assert!(!check(&tetra(0)[..3], 4));
        // Two tetrahedra sharing vertex 3: every edge is fine, the vertex link is not.
        let mut pinched = tetra(0).to_vec();
        pinched.extend(tetra(3));
        assert!(!check(&pinched, 7));
    }
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
    fn zero_level_roundoff_ties_do_not_create_collapsed_sphere_edges() {
        let mut reference: Option<(usize, Vec<Point>)> = None;
        for scale in [1e-12, 1., 1e12] {
            let field = |p: Point| scale * (dot(p, p).sqrt() - 1.);
            let band = Band {
                field: &field,
                bounds: [[-1.2; 3], [1.2; 3]],
                levels: [f64::NEG_INFINITY, 0.],
            };
            let mut mesh = extract(&band, 12).unwrap();
            assert!(validate(&mesh).unwrap().0 > 1e-8);
            if let Some((faces, points)) = &reference {
                assert_eq!(mesh.faces.len(), *faces);
                assert_eq!(mesh.points.len(), points.len());
                for (a, b) in mesh.points.iter().zip(points) {
                    assert!(crate::norm(crate::sub(*a, *b)) < 1e-11);
                }
            } else {
                reference = Some((mesh.faces.len(), mesh.points.clone()));
            }
            smooth(&mut mesh, &band, 10).unwrap();
            improve(&mut mesh, &band, 12).unwrap();
            polish(&mut mesh, &band, 40).unwrap();
            assert!(validate(&mesh).unwrap().0 > 1e-8);
            let volume = validate(&mesh).unwrap().1;
            assert!((volume / (4. * std::f64::consts::PI / 3.) - 1.).abs() < 0.03);
            assert!(sampled_error(&mesh, &band).unwrap() < 0.02);
        }
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
