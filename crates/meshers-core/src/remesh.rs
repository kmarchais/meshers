//! Curvature-adaptive remeshing of implicit triangle surfaces.
//!
//! Target edge lengths follow the chord sag of the exact surface. Two surface points
//! whose normals differ by an angle `theta` over a chord of length `l` bound an arc of
//! curvature about `theta / l`; a triangle of side `l` on such a surface sags by about
//! `kappa l^2 / 6` at its centroid. Edges are split, collapsed, flipped and relaxed
//! towards that length, and every vertex is projected back onto its exact constraint
//! (implicit wall, box plane, or their intersection curve). A final pass splits edges
//! until sampled centroid and midpoint distances meet the tolerance.
use rayon::prelude::*;

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
    /// Offset wall vertices outward by half the local chord sag, so triangles deviate
    /// to both sides of the surface. Halves the deviation at equal triangle count and
    /// removes the volume bias; vertices then no longer lie exactly on the surface.
    pub balanced: bool,
    /// Coarsen by collapses only, then split where the sampled distance exceeds the
    /// tolerance. No flips or smoothing: fewer triangles than a uniform extraction at
    /// a fraction of the full remeshing cost, with lower but valid triangle shapes.
    pub fast: bool,
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

/// An edge with its sort key (a length or a length ratio).
type Edge = (f64, [usize; 2]);

const SPLIT: f64 = 4. / 3.;
const COLLAPSE: f64 = 4. / 5.;
/// Fraction of the tolerance used for sizing, leaving room for estimate error.
const SIZING_SAFETY: f64 = 0.8;
/// Collapses must keep every changed triangle above this shape, or not worsen it.
const SHAPE_GUARD: f64 = 0.3;
const FAST_SHAPE_GUARD: f64 = 0.05;

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

/// Largest magnitude of the quadratic through corner values and the mid-edge values
/// of edges (0,1), (1,2), (2,0). The deviation of a flat triangle from a smooth
/// surface is quadratic to leading order, and its peak is generally at none of the
/// six samples, so this catches what they miss without more field evaluations.
fn quadratic_max(corners: [f64; 3], mids: [f64; 3]) -> f64 {
    // q(s, t) = a + b s + c t + d s² + e s t + f t², corner 0 at the origin, corner 1
    // at s = 1, corner 2 at t = 1.
    let [c0, c1, c2] = corners;
    let [m01, m12, m20] = mids;
    let (a, b, c) = (c0, 4. * m01 - 3. * c0 - c1, 4. * m20 - 3. * c0 - c2);
    let (d, f) = (2. * c1 - 4. * m01 + 2. * c0, 2. * c2 - 4. * m20 + 2. * c0);
    let e = 4. * (m12 - a - 0.5 * b - 0.5 * c - 0.25 * d - 0.25 * f);
    let q = |s: f64, t: f64| a + b * s + c * t + d * s * s + e * s * t + f * t * t;
    let mut worst = c0.abs().max(c1.abs()).max(c2.abs());
    let mut candidate = |s: f64, t: f64| {
        if s >= 0. && t >= 0. && s + t <= 1. + 1e-12 {
            worst = worst.max(q(s, t).abs());
        }
    };
    let det = 4. * d * f - e * e;
    if det != 0. {
        candidate((-2. * f * b + e * c) / det, (-2. * d * c + e * b) / det);
    }
    if d != 0. {
        candidate(-b / (2. * d), 0.);
    }
    if f != 0. {
        candidate(0., -c / (2. * f));
    }
    let curvature = d - e + f;
    if curvature != 0. {
        let s = -(b - c + e - 2. * f) / (2. * curvature);
        candidate(s, 1. - s);
    }
    worst
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
        if masks
            .iter()
            .any(|&m| m == 0 || (0..4).any(|j| (m >> (2 * j)) & 3 == 3))
        {
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
    fn face_shape(&self, face: [usize; 3], label: u8, moved: Option<(usize, Point, Point)>) -> f64 {
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
            } else if norm(self.normals[u]) > 0. && norm(self.normals[v]) > 0. && shared & 3 != 0 {
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
            [first, second].iter().all(|face| {
                self.face_shape(**face, self.labels[*f], Some((usize::MAX, p, n))) > 1e-6
            })
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
            let guard = if o.fast {
                FAST_SHAPE_GUARD
            } else {
                SHAPE_GUARD
            };
            if self.face_shape(next, self.labels[f], None) < before.min(guard) {
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
                .filter(|&w| {
                    self.masks[w] & mask == mask && {
                        let faces = self.edge_faces(v, w);
                        faces.len() == 2 && self.edge_mask(&faces) == mask
                    }
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

    /// First-order wall distance; `mask` box bits hold those axes, so a chord on a
    /// wall/plane curve is measured against the curve inside its plane.
    /// Signed outward version of the same estimate.
    fn signed_in(&self, p: Point, label: u8, mask: u8) -> f64 {
        let f = self.band.field.value(p);
        let residual = if label == 0 {
            f - self.band.levels[1]
        } else {
            self.band.levels[0] - f
        };
        if residual == 0. {
            return 0.;
        }
        let mut g = self.band.gradient(p);
        for (a, v) in g.iter_mut().enumerate() {
            if (mask >> (2 + 2 * a)) & 3 != 0 {
                *v = 0.;
            }
        }
        residual / norm(g)
    }

    /// Evaluate `eval` in parallel on alive wall faces that are new since the last
    /// call or `dirty`; other entries keep their cached value, dead faces get none.
    fn refresh<T: Send>(
        &self,
        cache: &mut Vec<Option<T>>,
        dirty: impl Fn(usize) -> bool + Sync,
        eval: impl Fn(usize) -> T + Sync,
    ) {
        cache.resize_with(self.faces.len(), || None);
        cache.par_iter_mut().enumerate().for_each(|(f, slot)| {
            if !self.alive[f] || self.labels[f] >= 2 {
                *slot = None;
            } else if slot.is_none() || dirty(f) {
                *slot = Some(eval(f));
            }
        });
    }

    /// Distance of a wall face with vertices on the surface, and the edges to split
    /// when it exceeds `tolerance`.
    fn wall_error(&self, f: usize, tolerance: f64) -> (f64, Vec<Edge>) {
        let face = self.faces[f];
        let p = face.map(|i| self.points[i]);
        let centroid = std::array::from_fn(|k| (p[0][k] + p[1][k] + p[2][k]) / 3.);
        let dc = self.signed_in(centroid, self.labels[f], 0).abs();
        let mut flagged = Vec::new();
        let mut longest = (0., [face[0], face[1]]);
        let mut mids = [0.; 3];
        for i in 0..3 {
            // Every edge of a wall face is sampled: an edge shared with a cap
            // face has no second wall face to sample it from.
            let (a, b) = (face[i], face[(i + 1) % 3]);
            let mid = std::array::from_fn(|k| 0.5 * (p[i][k] + p[(i + 1) % 3][k]));
            mids[i] = self.signed_in(mid, self.labels[f], self.masks[a] & self.masks[b]);
            let l = self.length(a, b);
            if l > longest.0 {
                longest = (l, [a, b]);
            }
            if mids[i].abs() > tolerance {
                flagged.push((l, [a, b]));
            }
        }
        // Vertices lie on the surface here.
        let worst = dc.max(quadratic_max([0.; 3], mids));
        if worst > tolerance && flagged.is_empty() {
            flagged.push(longest);
        }
        (worst, flagged)
    }

    /// Split until sampled midpoint and centroid distances meet the tolerance.
    fn enforce(&mut self, o: &Options, stats: &mut Stats) -> Result<()> {
        // A split changes only the faces it creates, so only those are sampled again.
        let mut cache = Vec::new();
        for _ in 0..16 {
            self.band
                .field
                .check()
                .map_err(MeshingError::GenerationFailed)?;
            self.refresh(&mut cache, |_| false, |f| self.wall_error(f, o.tolerance));
            let worst = cache.iter().flatten().map(|s| s.0).fold(0f64, f64::max);
            let flagged: Vec<_> = cache.iter().flatten().flat_map(|s| s.1.clone()).collect();
            stats.maximum_error = worst;
            if flagged.is_empty() {
                return Ok(());
            }
            let splits = self.split_all(flagged, o);
            stats.error_splits += splits;
            if splits == 0 {
                return Ok(());
            }
        }
        Ok(())
    }

    /// Signed outward distance from the wall that `label` belongs to.
    fn signed(&self, p: Point, label: u8) -> f64 {
        self.signed_in(p, label, 0)
    }

    /// Largest deviation at wall vertices, edge midpoints, centroids and the maximum
    /// of their quadratic interpolant, and the longest edge of every face whose
    /// deviation exceeds `tolerance`.
    /// Faces are re-sampled only when new or when one of their vertices moved since
    /// the positions in `previous`.
    fn deviation(
        &self,
        tolerance: f64,
        cache: &mut Vec<Option<(f64, Option<Edge>)>>,
        previous: &[Point],
    ) -> (f64, Vec<Edge>) {
        let moved = |f: usize| {
            self.faces[f]
                .iter()
                .any(|&v| previous.get(v) != Some(&self.points[v]))
        };
        self.refresh(cache, moved, |f| {
            let label = self.labels[f];
            let face = self.faces[f];
            let p = face.map(|i| self.points[i]);
            let centroid = std::array::from_fn(|k| (p[0][k] + p[1][k] + p[2][k]) / 3.);
            let mut longest = (0., [face[0], face[1]]);
            let mut corners = [0.; 3];
            let mut mids = [0.; 3];
            for i in 0..3 {
                let mid = std::array::from_fn(|k| 0.5 * (p[i][k] + p[(i + 1) % 3][k]));
                corners[i] = self.signed(p[i], label);
                mids[i] = self.signed(mid, label);
                let edge = [face[i], face[(i + 1) % 3]];
                let l = self.length(edge[0], edge[1]);
                if l > longest.0 {
                    longest = (l, edge);
                }
            }
            let worst = self
                .signed(centroid, label)
                .abs()
                .max(quadratic_max(corners, mids));
            (worst, (worst > tolerance).then_some(longest))
        });
        let worst = cache.iter().flatten().map(|s| s.0).fold(0., f64::max);
        (worst, cache.iter().flatten().filter_map(|s| s.1).collect())
    }

    /// Split edges longest first; returns how many splits succeeded.
    fn split_all(&mut self, mut edges: Vec<Edge>, o: &Options) -> usize {
        edges.sort_by(|a, b| b.0.total_cmp(&a.0));
        edges.dedup_by_key(|e| {
            let [a, b] = e.1;
            [a.min(b), a.max(b)]
        });
        let mut count = 0;
        for (_, [a, b]) in edges {
            if self.edge_faces(a, b).len() != 2 || self.length(a, b) < o.min_length {
                continue;
            }
            if self.split(a, b).is_some() {
                count += 1;
            }
        }
        count
    }

    /// Move wall vertices outward by half the mean centroid sag of their faces, holding
    /// box coordinates. Returns the on-surface positions for undoing.
    /// Centroid sags are cached in `sags`: on-surface positions of surviving faces do
    /// not change between calls.
    fn balance(&mut self, sags: &mut Vec<Option<f64>>) -> Vec<Point> {
        let base = self.points.clone();
        let n = self.points.len();
        let mut sum = vec![0.; n];
        let mut count = vec![0usize; n];
        self.refresh(
            sags,
            |_| false,
            |f| {
                let p = self.faces[f].map(|i| self.points[i]);
                let c = std::array::from_fn(|k| (p[0][k] + p[1][k] + p[2][k]) / 3.);
                self.signed(c, self.labels[f])
            },
        );
        for (f, sag) in sags.iter().enumerate() {
            let Some(sag) = *sag else { continue };
            for &v in &self.faces[f] {
                sum[v] += sag;
                count[v] += 1;
            }
        }
        for v in 0..n {
            let mask = self.masks[v];
            if count[v] == 0 || mask & 3 == 0 || mask.count_ones() >= 3 {
                continue;
            }
            let mut direction = self.normals[v];
            for (a, d) in direction.iter_mut().enumerate() {
                if (mask >> (2 + 2 * a)) & 3 != 0 {
                    *d = 0.;
                }
            }
            let direction = unit(direction);
            let offset = -0.5 * sum[v] / count[v] as f64;
            self.points[v] = std::array::from_fn(|k| base[v][k] + offset * direction[k]);
        }
        base
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
        || o.gradation.is_nan()
        || o.gradation <= 0.
        || o.iterations > 100
    {
        return Err(MeshingError::InvalidOptions(
            "positive tolerance, 0 < min_length <= max_length, positive gradation, iterations <= 100"
                .into(),
        ));
    }
    let mut s = Surface::new(mesh, band)?;
    let mut stats = Stats::default();
    let sizing = Options {
        tolerance: if o.balanced {
            1.8 * o.tolerance
        } else {
            o.tolerance
        },
        ..*o
    };
    let o_final = o;
    let o = &sizing;
    let mut clock = std::time::Instant::now();
    let mut lap = |stats: &mut Stats, phase: usize| {
        stats.seconds[phase] += clock.elapsed().as_secs_f64();
        clock = std::time::Instant::now();
    };
    for _ in 0..o.iterations {
        band.field.check().map_err(MeshingError::GenerationFailed)?;
        if o.fast {
            let size = s.sizes(o);
            lap(&mut stats, 0);
            stats.collapses += s.coarsen(&size, o);
            lap(&mut stats, 2);
            continue;
        }
        let size = s.sizes(o);
        lap(&mut stats, 0);
        let mut long: Vec<Edge> = s
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
        stats.collapses += s.coarsen(&size, o);
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
    if o_final.balanced {
        // Chords may sag up to about twice the tolerance before balancing halves it.
        let relaxed = Options {
            tolerance: 1.8 * o_final.tolerance,
            ..*o_final
        };
        s.enforce(&relaxed, &mut stats)?;
        let mut done = false;
        let (mut sags, mut deviations, mut previous) = (Vec::new(), Vec::new(), Vec::new());
        for _ in 0..16 {
            band.field.check().map_err(MeshingError::GenerationFailed)?;
            let base = s.balance(&mut sags);
            let (deviation, flagged) = s.deviation(o_final.tolerance, &mut deviations, &previous);
            previous.clone_from(&s.points);
            if flagged.is_empty() {
                stats.maximum_error = deviation;
                done = true;
                break;
            }
            // Refine only where the balanced surface still deviates, on the exact
            // positions, then balance again.
            s.points = base;
            let splits = s.split_all(flagged, o_final);
            stats.error_splits += splits;
            if splits == 0 {
                break;
            }
        }
        if !done {
            s.enforce(o_final, &mut stats)?;
        }
    } else {
        s.enforce(o_final, &mut stats)?;
    }
    *mesh = s.compact();
    crate::triangles::validate(mesh)?;
    lap(&mut stats, 5);
    Ok(stats)
}

impl<F: ScalarField> Surface<'_, '_, F> {
    /// Collapse edges shorter than their target, shortest first. Returns the count.
    fn coarsen(&mut self, size: &[f64], o: &Options) -> usize {
        let mut short: Vec<Edge> = self
            .edges()
            .into_iter()
            .filter_map(|[u, v]| {
                let ratio = self.length(u, v) / Self::target(size, u, v);
                (ratio < COLLAPSE).then_some((ratio, [u, v]))
            })
            .collect();
        short.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut count = 0;
        for (_, [u, v]) in short {
            if self.vf[u].is_empty() || self.vf[v].is_empty() || self.edge_faces(u, v).len() != 2 {
                continue;
            }
            // Keep the more constrained vertex.
            let (remove, keep) = if self.masks[u].count_ones() <= self.masks[v].count_ones() {
                (u, v)
            } else {
                (v, u)
            };
            if self.collapse(remove, keep, size, o) || self.collapse(keep, remove, size, o) {
                count += 1;
            }
        }
        count
    }

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
            balanced: false,
            fast: false,
        }
    }

    #[test]
    fn quadratic_max_matches_dense_sampling() {
        let cases: [fn(f64, f64) -> f64; 3] = [
            |s, t| 0.1 + s - 2. * t * t + 3. * s * t,
            |s, t| -(s - 0.3).powi(2) - (t - 0.2).powi(2) + 0.01,
            |s, t| 2. * s * s - s + 0.5 * t,
        ];
        for q in cases {
            let corners = [q(0., 0.), q(1., 0.), q(0., 1.)];
            let mids = [q(0.5, 0.), q(0.5, 0.5), q(0., 0.5)];
            let n = 400;
            let dense = (0..=n)
                .flat_map(|i| (0..=n - i).map(move |j| (i as f64 / n as f64, j as f64 / n as f64)))
                .map(|(s, t)| q(s, t).abs())
                .fold(0., f64::max);
            let exact = quadratic_max(corners, mids);
            assert!(
                exact >= dense - 1e-12 && exact <= dense + 1e-5,
                "{exact} vs {dense}"
            );
        }
    }

    #[test]
    fn fast_mode_meets_tolerance_with_fewer_triangles() {
        let field = |p: Point| dot(p, p).sqrt() - 0.35;
        let band = Band {
            field: &field,
            bounds: [[-0.5; 3], [0.5; 3]],
            levels: [f64::NEG_INFINITY, 0.],
        };
        let mesh = extract(&band, 24).unwrap();
        let extracted = mesh.faces.len();
        for balanced in [false, true] {
            let mut out = TriangleMesh {
                points: mesh.points.clone(),
                faces: mesh.faces.clone(),
                labels: mesh.labels.clone(),
            };
            let o = Options {
                iterations: 2,
                balanced,
                fast: true,
                ..options(1e-3)
            };
            adapt(&mut out, &band, &o).unwrap();
            validate(&out).unwrap();
            assert!(
                out.faces.len() * 2 < extracted,
                "{} of {extracted}",
                out.faces.len()
            );
            if !balanced {
                assert!(sampled_error(&out, &band).unwrap() <= 1e-3);
            }
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
        assert!(
            mesh.faces.len() < before,
            "{} vs {before}",
            mesh.faces.len()
        );
        assert!(minimum > 0.3, "minimum shape {minimum}");
        let exact = 4. / 3. * std::f64::consts::PI * 0.35f64.powi(3);
        assert!((volume / exact - 1.).abs() < 0.01);
    }

    #[test]
    fn balanced_sphere_halves_deviation_and_volume_bias() {
        let field = |p: Point| dot(p, p).sqrt() - 0.35;
        let band = Band {
            field: &field,
            bounds: [[-0.5; 3], [0.5; 3]],
            levels: [f64::NEG_INFINITY, 0.],
        };
        let exact = 4. / 3. * std::f64::consts::PI * 0.35f64.powi(3);
        let mut plain = extract(&band, 24).unwrap();
        adapt(&mut plain, &band, &options(1e-3)).unwrap();
        let mut balanced = extract(&band, 24).unwrap();
        let stats = adapt(
            &mut balanced,
            &band,
            &Options {
                balanced: true,
                ..options(1e-3)
            },
        )
        .unwrap();
        assert!(stats.maximum_error <= 1e-3, "{stats:?}");
        // Same tolerance with fewer triangles. Half the centroid sag minimizes the
        // largest deviation; the volume bias shrinks but does not vanish.
        assert!(balanced.faces.len() * 10 < plain.faces.len() * 7);
        let bias = |m: &TriangleMesh| (validate(m).unwrap().1 / exact - 1.).abs();
        assert!(
            bias(&balanced) < 0.75 * bias(&plain),
            "{} {}",
            bias(&balanced),
            bias(&plain)
        );
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
