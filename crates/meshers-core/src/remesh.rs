//! Curvature-adaptive remeshing of implicit triangle surfaces.
//!
//! Target edge lengths follow the chord sag of the exact surface. Two surface points
//! whose normals differ by an angle `theta` over a chord of length `l` bound an arc of
//! curvature about `theta / l`; a triangle of side `l` on such a surface sags by about
//! `kappa l^2 / 6` at its centroid. Edges are split, collapsed, flipped and relaxed
//! towards that length, and every vertex is projected back onto its exact constraint
//! (implicit wall, box plane, or their intersection curve). A final pass splits edges
//! until sampled centroid and midpoint distances meet the tolerance.
use crate::{
    MeshingError, Point, cross, dot, implicit::ScalarField, norm, quads::Band, sub,
    triangles::TriangleMesh,
};

type Result<T> = std::result::Result<T, MeshingError>;

#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Maximum surface distance of the triangles, in physical units.
    pub tolerance: f64,
    pub min_length: f64,
    pub max_length: f64,
    /// Maximum growth of the target length per unit distance (0.3 is typical).
    pub gradation: f64,
    pub iterations: usize,
}

#[derive(Clone, Debug, Default)]
pub struct Stats {
    pub splits: usize,
    pub collapses: usize,
    pub flips: usize,
    pub moves: usize,
    pub error_splits: usize,
    pub maximum_error: f64,
    /// Seconds per phase: sizing, split, collapse, flip, relax, enforce+validate.
    pub seconds: [f64; 6],
}

const SPLIT: f64 = 4. / 3.;
const COLLAPSE: f64 = 4. / 5.;
/// Fraction of the tolerance used for sizing, leaving room for estimate error.
const SIZING_SAFETY: f64 = 0.8;

fn fail(message: &str) -> MeshingError {
    MeshingError::GenerationFailed(message.into())
}

fn unit(a: Point) -> Point {
    let n = norm(a);
    if n > 0. { a.map(|v| v / n) } else { a }
}

fn angle(a: Point, b: Point) -> f64 {
    norm(cross(a, b)).atan2(dot(a, b))
}

fn shape(p: [Point; 3], normal: Point) -> f64 {
    let e = sub(p[1], p[0]);
    let d = sub(p[2], p[0]);
    let f = sub(p[2], p[1]);
    2. * 3f64.sqrt() * dot(cross(e, d), normal) / (dot(e, e) + dot(d, d) + dot(f, f)).max(1e-300)
}

fn wall_label(mask: u8) -> Option<u8> {
    if mask & 1 != 0 {
        Some(0)
    } else if mask & 2 != 0 {
        Some(1)
    } else {
        None
    }
}

struct Surface<'a, 'b, F: ScalarField> {
    band: &'a Band<'b, F>,
    points: Vec<Point>,
    masks: Vec<u8>,
    /// Unit wall normal for vertices on an implicit wall, zero otherwise.
    normals: Vec<Point>,
    vf: Vec<Vec<usize>>,
    faces: Vec<[usize; 3]>,
    labels: Vec<u8>,
    alive: Vec<bool>,
}

impl<'a, 'b, F: ScalarField> Surface<'a, 'b, F> {
    fn new(mesh: &TriangleMesh, band: &'a Band<'b, F>) -> Result<Self> {
        let n = mesh.points.len();
        let mut masks = vec![0u8; n];
        let mut vf = vec![Vec::new(); n];
        for (f, (face, &label)) in mesh.faces.iter().zip(&mesh.labels).enumerate() {
            for &v in face {
                masks[v] |= 1 << label;
                vf[v].push(f);
            }
        }
        if masks.iter().any(|&m| m == 0 || (0..4).any(|j| (m >> (2 * j)) & 3 == 3)) {
            return Err(fail("Unresolved or incompatible surface constraints"));
        }
        let mut surface = Self {
            band,
            points: mesh.points.clone(),
            masks,
            normals: vec![[0.; 3]; n],
            vf,
            faces: mesh.faces.clone(),
            labels: mesh.labels.clone(),
            alive: vec![true; mesh.faces.len()],
        };
        for v in 0..n {
            surface.normals[v] = surface.wall_normal(surface.points[v], surface.masks[v]);
        }
        Ok(surface)
    }

    fn wall_normal(&self, p: Point, mask: u8) -> Point {
        wall_label(mask).map_or([0.; 3], |l| self.band.normal(p, l))
    }

    /// Shape against the reference normal: cached vertex wall normals, or the box axis.
    /// `moved` substitutes one vertex (or the new vertex `usize::MAX`) and its normal.
    fn face_shape(
        &self,
        face: [usize; 3],
        label: u8,
        moved: Option<(usize, Point, Point)>,
    ) -> f64 {
        let at = |i: usize| match moved {
            Some((j, q, n)) if i == j => (q, n),
            _ => (self.points[i], self.normals[i]),
        };
        let p = face.map(|i| at(i).0);
        let normal = if label >= 2 {
            self.band.normal(p[0], label)
        } else {
            unit(face.iter().fold([0.; 3], |m, &i| {
                let n = at(i).1;
                [m[0] + n[0], m[1] + n[1], m[2] + n[2]]
            }))
        };
        shape(p, normal)
    }

    fn edge_faces(&self, u: usize, v: usize) -> Vec<usize> {
        self.vf[u]
            .iter()
            .copied()
            .filter(|&f| self.faces[f].contains(&v))
            .collect()
    }

    fn neighbors(&self, u: usize) -> Vec<usize> {
        let mut out: Vec<usize> = self.vf[u]
            .iter()
            .flat_map(|&f| self.faces[f])
            .filter(|&w| w != u)
            .collect();
        out.sort_unstable();
        out.dedup();
        out
    }

    /// Undirected edges, each once: an oriented closed surface holds every edge in
    /// exactly one face in increasing order.
    fn edges(&self) -> Vec<[usize; 2]> {
        let mut out = Vec::new();
        for (f, face) in self.faces.iter().enumerate() {
            if !self.alive[f] {
                continue;
            }
            for i in 0..3 {
                let (a, b) = (face[i], face[(i + 1) % 3]);
                if a < b {
                    out.push([a, b]);
                }
            }
        }
        out
    }

    fn length(&self, u: usize, v: usize) -> f64 {
        norm(sub(self.points[u], self.points[v]))
    }

    /// Sag estimate of the chord between two vertices on the same wall.
    fn sag(&self, u: usize, v: usize, p: Option<(usize, Point, Point)>) -> f64 {
        let (pu, nu) = match p {
            Some((j, q, n)) if j == u => (q, n),
            _ => (self.points[u], self.normals[u]),
        };
        let (pv, nv) = match p {
            Some((j, q, n)) if j == v => (q, n),
            _ => (self.points[v], self.normals[v]),
        };
        if norm(nu) == 0. || norm(nv) == 0. {
            return 0.;
        }
        norm(sub(pu, pv)) * angle(nu, nv) / 8.
    }

    fn curve_tangent(&self, v: usize) -> Point {
        let cap = (2..8).find(|&l| self.masks[v] & (1 << l) != 0).unwrap_or(2);
        let mut n = [0.; 3];
        n[(cap as usize - 2) / 2] = 1.;
        unit(cross(self.normals[v], n))
    }

    fn sizes(&self, o: &Options) -> Vec<f64> {
        let tol = o.tolerance * SIZING_SAFETY;
        let mut size = vec![o.max_length; self.points.len()];
        let edges = self.edges();
        for &[u, v] in &edges {
            let l = self.length(u, v);
            if l <= 0. {
                continue;
            }
            let shared = self.masks[u] & self.masks[v];
            let target = if shared.count_ones() >= 2 && shared & 3 != 0 {
                // Feature curve: its chord sags with the curve's own curvature.
                let t = self.curve_tangent(u);
                let s = self.curve_tangent(v);
                let theta = norm(cross(t, s)).atan2(dot(t, s).abs());
                (8. * tol * l / theta.max(1e-12)).sqrt()
            } else if norm(self.normals[u]) > 0. && norm(self.normals[v]) > 0. && shared & 3 != 0
            {
                let theta = angle(self.normals[u], self.normals[v]);
                (6. * tol * l / theta.max(1e-12)).sqrt()
            } else {
                o.max_length
            };
            let target = target.clamp(o.min_length, o.max_length);
            size[u] = size[u].min(target);
            size[v] = size[v].min(target);
        }
        // Lipschitz gradation, so neighbouring triangles change size gradually.
        for _ in 0..64 {
            let mut changed = false;
            for &[u, v] in &edges {
                let l = self.length(u, v);
                for (a, b) in [(u, v), (v, u)] {
                    let limit = size[a] + o.gradation * l;
                    if size[b] > limit * (1. + 1e-12) {
                        size[b] = limit;
                        changed = true;
                    }
                }
            }
            if !changed {
                break;
            }
        }
        size
    }

    fn target(size: &[f64], u: usize, v: usize) -> f64 {
        0.5 * (size[u] + size[v])
    }

    fn edge_mask(&self, faces: &[usize]) -> u8 {
        faces.iter().fold(0, |m, &f| m | (1 << self.labels[f]))
    }

    fn add_vertex(&mut self, p: Point, mask: u8) -> usize {
        self.points.push(p);
        self.masks.push(mask);
        self.normals.push(self.wall_normal(p, mask));
        self.vf.push(Vec::new());
        self.points.len() - 1
    }

    fn add_face(&mut self, face: [usize; 3], label: u8) -> usize {
        let f = self.faces.len();
        self.faces.push(face);
        self.labels.push(label);
        self.alive.push(true);
        for v in face {
            self.vf[v].push(f);
        }
        f
    }

    fn kill_face(&mut self, f: usize) {
        self.alive[f] = false;
        for v in self.faces[f] {
            self.vf[v].retain(|&g| g != f);
        }
    }

    /// Rotate a face so that the edge `{u, v}` comes first; returns `[x, y, opposite]`.
    fn rotate(face: [usize; 3], u: usize, v: usize) -> [usize; 3] {
        for i in 0..3 {
            let (a, b) = (face[i], face[(i + 1) % 3]);
            if (a == u && b == v) || (a == v && b == u) {
                return [a, b, face[(i + 2) % 3]];
            }
        }
        unreachable!("edge not in face")
    }

    fn split(&mut self, u: usize, v: usize) -> Option<usize> {
        let faces = self.edge_faces(u, v);
        if faces.len() != 2 {
            return None;
        }
        let mask = self.edge_mask(&faces);
        if mask & !(self.masks[u] & self.masks[v]) != 0 {
            return None;
        }
        let mid = std::array::from_fn(|k| 0.5 * (self.points[u][k] + self.points[v][k]));
        let p = self.band.project(mid, mask).ok()?;
        let mut planned = Vec::new();
        for &f in &faces {
            let [x, y, a] = Self::rotate(self.faces[f], u, v);
            planned.push((f, [x, usize::MAX, a], [usize::MAX, y, a]));
        }
        let n = self.wall_normal(p, mask);
        let ok = planned.iter().all(|(f, first, second)| {
            [first, second]
                .iter()
                .all(|face| self.face_shape(**face, self.labels[*f], Some((usize::MAX, p, n))) > 1e-6)
        });
        if !ok {
            return None;
        }
        let m = self.add_vertex(p, mask);
        for (f, first, second) in planned {
            let label = self.labels[f];
            self.kill_face(f);
            self.add_face(first.map(|i| if i == usize::MAX { m } else { i }), label);
            self.add_face(second.map(|i| if i == usize::MAX { m } else { i }), label);
        }
        Some(m)
    }

    /// Remove `u` by merging it into `v`.
    fn collapse(&mut self, u: usize, v: usize, size: &[f64], o: &Options) -> bool {
        let mu = self.masks[u];
        if mu.count_ones() >= 3 || mu & !self.masks[v] != 0 {
            return false;
        }
        let faces = self.edge_faces(u, v);
        if faces.len() != 2 || self.edge_mask(&faces) != mu {
            return false;
        }
        let opposite: Vec<usize> = faces
            .iter()
            .map(|&f| Self::rotate(self.faces[f], u, v)[2])
            .collect();
        let nu = self.neighbors(u);
        let nv = self.neighbors(v);
        let common: Vec<usize> = nu.iter().copied().filter(|w| nv.contains(w)).collect();
        if common.len() != 2 || !opposite.iter().all(|w| common.contains(w)) {
            return false;
        }
        // On a closed manifold a vertex has as many neighbours as incident faces.
        if opposite.iter().any(|&w| self.vf[w].len() <= 3) {
            return false;
        }
        for &w in &nu {
            if w == v {
                continue;
            }
            if self.length(v, w) > SPLIT * Self::target(size, v, w) {
                return false;
            }
            if self.sag(v, w, None) > o.tolerance * SIZING_SAFETY {
                return false;
            }
        }
        let before = self.vf[u]
            .iter()
            .map(|&f| self.face_shape(self.faces[f], self.labels[f], None))
            .fold(1f64, f64::min);
        for &f in &self.vf[u] {
            if faces.contains(&f) {
                continue;
            }
            let next = self.faces[f].map(|i| if i == u { v } else { i });
            if self.face_shape(next, self.labels[f], None) < before.min(0.3) {
                return false;
            }
        }
        for f in faces {
            self.kill_face(f);
        }
        for f in std::mem::take(&mut self.vf[u]) {
            self.faces[f] = self.faces[f].map(|i| if i == u { v } else { i });
            self.vf[v].push(f);
        }
        true
    }

    fn valence_target(&self, v: usize) -> i64 {
        match self.masks[v].count_ones() {
            1 => 6,
            2 => 4,
            _ => 3,
        }
    }

    fn flip(&mut self, u: usize, v: usize, o: &Options) -> bool {
        let faces = self.edge_faces(u, v);
        if faces.len() != 2 || self.labels[faces[0]] != self.labels[faces[1]] {
            return false;
        }
        let label = self.labels[faces[0]];
        let f0 = Self::rotate(self.faces[faces[0]], u, v);
        let f1 = Self::rotate(self.faces[faces[1]], u, v);
        let (x, y, a, b) = (f0[0], f0[1], f0[2], f1[2]);
        if a == b || self.vf[a].iter().any(|&f| self.faces[f].contains(&b)) {
            return false;
        }
        let next = [[x, b, a], [b, y, a]];
        let valence = |w: usize| self.vf[w].len() as i64;
        let deviation = |w: usize, delta: i64| (valence(w) + delta - self.valence_target(w)).pow(2);
        let gain = [x, y, a, b].iter().map(|&w| deviation(w, 0)).sum::<i64>()
            - (deviation(x, -1) + deviation(y, -1) + deviation(a, 1) + deviation(b, 1));
        let before = [f0, f1]
            .iter()
            .map(|&f| self.face_shape(f, label, None))
            .fold(1f64, f64::min);
        let after = next
            .iter()
            .map(|&f| self.face_shape(f, label, None))
            .fold(1f64, f64::min);
        if after <= 1e-6 || (label < 2 && self.sag(a, b, None) > o.tolerance * SIZING_SAFETY) {
            return false;
        }
        let accept = (gain > 0 && after >= before.min(0.4)) || after > before + 0.05;
        if !accept {
            return false;
        }
        self.kill_face(faces[0]);
        self.kill_face(faces[1]);
        self.add_face(next[0], label);
        self.add_face(next[1], label);
        true
    }

    fn relax(&mut self, v: usize) -> bool {
        let mask = self.masks[v];
        let bits = mask.count_ones();
        if bits >= 3 || self.vf[v].is_empty() {
            return false;
        }
        let p = self.points[v];
        let target = if bits == 1 {
            // Area-weighted barycentre of the one-ring, restricted to the tangent plane.
            let mut sum = [0.; 3];
            let mut weight = 0.;
            for &f in &self.vf[v] {
                let q = self.faces[f].map(|i| self.points[i]);
                let area = 0.5 * norm(cross(sub(q[1], q[0]), sub(q[2], q[0])));
                for k in 0..3 {
                    sum[k] += area * (q[0][k] + q[1][k] + q[2][k]) / 3.;
                }
                weight += area;
            }
            if weight <= 0. {
                return false;
            }
            let c: Point = sum.map(|s| s / weight);
            let n = if mask & 3 != 0 {
                self.normals[v]
            } else {
                let cap = (2..8).find(|&l| mask & (1 << l) != 0).unwrap();
                let mut n = [0.; 3];
                n[(cap as usize - 2) / 2] = 1.;
                n
            };
            let d = sub(c, p);
            let d: Point = std::array::from_fn(|k| d[k] - dot(d, n) * n[k]);
            std::array::from_fn(|k| p[k] + 0.5 * d[k])
        } else {
            // Curve vertex: move towards the midpoint of its two curve neighbours.
            let curve: Vec<usize> = self
                .neighbors(v)
                .into_iter()
                .filter(|&w| self.masks[w] & mask == mask && {
                    let faces = self.edge_faces(v, w);
                    faces.len() == 2 && self.edge_mask(&faces) == mask
                })
                .collect();
            if curve.len() != 2 {
                return false;
            }
            let (a, b) = (self.points[curve[0]], self.points[curve[1]]);
            std::array::from_fn(|k| p[k] + 0.5 * (0.5 * (a[k] + b[k]) - p[k]))
        };
        let Ok(q) = self.band.project(target, mask) else {
            return false;
        };
        let n = self.wall_normal(q, mask);
        let before = self.vf[v]
            .iter()
            .map(|&f| self.face_shape(self.faces[f], self.labels[f], None))
            .fold(1f64, f64::min);
        let after = self.vf[v]
            .iter()
            .map(|&f| self.face_shape(self.faces[f], self.labels[f], Some((v, q, n))))
            .fold(1f64, f64::min);
        if after < before.min(0.4) || after <= 1e-6 {
            return false;
        }
        self.points[v] = q;
        self.normals[v] = n;
        true
    }

    fn distance(&self, p: Point, label: u8) -> f64 {
        let level = self.band.levels[if label == 0 { 1 } else { 0 }];
        let residual = (self.band.field.value(p) - level).abs();
        if residual == 0. {
            return 0.;
        }
        residual / norm(self.band.gradient(p))
    }

    /// Split until sampled midpoint and centroid distances meet the tolerance.
    fn enforce(&mut self, o: &Options, stats: &mut Stats) -> Result<()> {
        for _ in 0..16 {
            self.band.field.check().map_err(MeshingError::GenerationFailed)?;
            let mut flagged = Vec::new();
            let mut worst = 0f64;
            for f in 0..self.faces.len() {
                if !self.alive[f] || self.labels[f] >= 2 {
                    continue;
                }
                let face = self.faces[f];
                let p = face.map(|i| self.points[i]);
                let centroid = std::array::from_fn(|k| (p[0][k] + p[1][k] + p[2][k]) / 3.);
                let dc = self.distance(centroid, self.labels[f]);
                worst = worst.max(dc);
                let mut longest = (0., [face[0], face[1]]);
                let mut any = false;
                for i in 0..3 {
                    let (a, b) = (face[i], face[(i + 1) % 3]);
                    if a > b {
                        // Each wall edge is sampled once, from its increasing face.
                        let l = self.length(a, b);
                        if l > longest.0 {
                            longest = (l, [a, b]);
                        }
                        continue;
                    }
                    let mid = std::array::from_fn(|k| 0.5 * (p[i][k] + p[(i + 1) % 3][k]));
                    let d = self.distance(mid, self.labels[f]);
                    worst = worst.max(d);
                    let l = self.length(a, b);
                    if l > longest.0 {
                        longest = (l, [a, b]);
                    }
                    if d > o.tolerance {
                        flagged.push((l, [a, b]));
                        any = true;
                    }
                }
                if dc > o.tolerance && !any {
                    flagged.push(longest);
                }
            }
            stats.maximum_error = worst;
            if flagged.is_empty() {
                return Ok(());
            }
            flagged.sort_by(|a, b| b.0.total_cmp(&a.0));
            let mut progress = false;
            for (_, [a, b]) in flagged {
                if self.length(a, b) < o.min_length {
                    continue;
                }
                if self.split(a, b).is_some() {
                    stats.error_splits += 1;
                    progress = true;
                }
            }
            if !progress {
                return Ok(());
            }
        }
        Ok(())
    }

    fn compact(&self) -> TriangleMesh {
        let mut map = vec![usize::MAX; self.points.len()];
        let mut mesh = TriangleMesh {
            points: Vec::new(),
            faces: Vec::new(),
            labels: Vec::new(),
        };
        for (f, face) in self.faces.iter().enumerate() {
            if !self.alive[f] {
                continue;
            }
            let face = face.map(|v| {
                if map[v] == usize::MAX {
                    map[v] = mesh.points.len();
                    mesh.points.push(self.points[v]);
                }
                map[v]
            });
            mesh.faces.push(face);
            mesh.labels.push(self.labels[f]);
        }
        mesh
    }
}

/// Remesh a closed labelled surface towards curvature-adaptive edge lengths.
/// The input must come from [`crate::triangles::extract`] or an equivalent labelled
/// surface whose vertices satisfy their constraints.
pub fn adapt<F: ScalarField>(
    mesh: &mut TriangleMesh,
    band: &Band<'_, F>,
    o: &Options,
) -> Result<Stats> {
    band.validate()?;
    if !(o.tolerance > 0. && o.min_length > 0. && o.max_length >= o.min_length)
        || !(o.gradation > 0.)
        || o.iterations > 100
    {
        return Err(MeshingError::InvalidOptions(
            "positive tolerance, 0 < min_length <= max_length, positive gradation, iterations <= 100"
                .into(),
        ));
    }
    let mut s = Surface::new(mesh, band)?;
    let mut stats = Stats::default();
    let mut clock = std::time::Instant::now();
    let mut lap = |stats: &mut Stats, phase: usize| {
        stats.seconds[phase] += clock.elapsed().as_secs_f64();
        clock = std::time::Instant::now();
    };
    for _ in 0..o.iterations {
        band.field.check().map_err(MeshingError::GenerationFailed)?;
        let size = s.sizes(o);
        lap(&mut stats, 0);
        let mut long: Vec<(f64, [usize; 2])> = s
            .edges()
            .into_iter()
            .filter_map(|[u, v]| {
                let ratio = s.length(u, v) / Surface::<F>::target(&size, u, v);
                (ratio > SPLIT).then_some((ratio, [u, v]))
            })
            .collect();
        long.sort_by(|a, b| b.0.total_cmp(&a.0));
        for (_, [u, v]) in long {
            if s.split(u, v).is_some() {
                stats.splits += 1;
            }
        }
        lap(&mut stats, 1);
        let size = s.sizes(o);
        lap(&mut stats, 0);
        let mut short: Vec<(f64, [usize; 2])> = s
            .edges()
            .into_iter()
            .filter_map(|[u, v]| {
                let ratio = s.length(u, v) / Surface::<F>::target(&size, u, v);
                (ratio < COLLAPSE).then_some((ratio, [u, v]))
            })
            .collect();
        short.sort_by(|a, b| a.0.total_cmp(&b.0));
        for (_, [u, v]) in short {
            if s.vf[u].is_empty() || s.vf[v].is_empty() || s.edge_faces(u, v).len() != 2 {
                continue;
            }
            // Keep the more constrained vertex.
            let (remove, keep) = if s.masks[u].count_ones() <= s.masks[v].count_ones() {
                (u, v)
            } else {
                (v, u)
            };
            if s.collapse(remove, keep, &size, o) || s.collapse(keep, remove, &size, o) {
                stats.collapses += 1;
            }
        }
        lap(&mut stats, 2);
        for [u, v] in s.edges() {
            if s.alive_edge(u, v) && s.flip(u, v, o) {
                stats.flips += 1;
            }
        }
        lap(&mut stats, 3);
        for v in 0..s.points.len() {
            if s.relax(v) {
                stats.moves += 1;
            }
        }
        lap(&mut stats, 4);
    }
    s.enforce(o, &mut stats)?;
    *mesh = s.compact();
    crate::triangles::validate(mesh)?;
    lap(&mut stats, 5);
    Ok(stats)
}

impl<F: ScalarField> Surface<'_, '_, F> {
    fn alive_edge(&self, u: usize, v: usize) -> bool {
        self.edge_faces(u, v).len() == 2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::triangles::{extract, sampled_error, validate};

    fn options(tolerance: f64) -> Options {
        Options {
            tolerance,
            min_length: 1e-3,
            max_length: 0.5,
            gradation: 0.3,
            iterations: 8,
        }
    }

    #[test]
    fn sphere_meets_tolerance_with_fewer_triangles() {
        let field = |p: Point| dot(p, p).sqrt() - 0.35;
        let band = Band {
            field: &field,
            bounds: [[-0.5; 3], [0.5; 3]],
            levels: [f64::NEG_INFINITY, 0.],
        };
        let mut mesh = extract(&band, 32).unwrap();
        let before = mesh.faces.len();
        let stats = adapt(&mut mesh, &band, &options(1e-3)).unwrap();
        let (minimum, volume) = validate(&mesh).unwrap();
        assert!(sampled_error(&mesh, &band).unwrap() <= 1e-3, "{stats:?}");
        assert!(mesh.faces.len() < before, "{} vs {before}", mesh.faces.len());
        assert!(minimum > 0.3, "minimum shape {minimum}");
        let exact = 4. / 3. * std::f64::consts::PI * 0.35f64.powi(3);
        assert!((volume / exact - 1.).abs() < 0.01);
    }

    #[test]
    fn clipped_band_keeps_caps_and_curves() {
        let field = |p: Point| {
            let t = std::f64::consts::TAU;
            ((t * p[0]).sin() * (t * p[1]).cos()
                + (t * p[1]).sin() * (t * p[2]).cos()
                + (t * p[2]).sin() * (t * p[0]).cos())
                / 0.55
        };
        let band = Band {
            field: &field,
            bounds: [[-0.5; 3], [0.5; 3]],
            levels: [-1., 1.],
        };
        let mut mesh = extract(&band, 24).unwrap();
        adapt(&mut mesh, &band, &options(2e-3)).unwrap();
        validate(&mesh).unwrap();
        assert!(sampled_error(&mesh, &band).unwrap() <= 2e-3);
        for (face, &label) in mesh.faces.iter().zip(&mesh.labels) {
            if label >= 2 {
                let axis = (label as usize - 2) / 2;
                let side = (label as usize - 2) % 2;
                for &v in face {
                    assert_eq!(mesh.points[v][axis], band.bounds[side][axis]);
                }
            }
        }
    }
}
