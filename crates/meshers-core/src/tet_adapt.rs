//! Size-driven tetrahedral adaptation with exact implicit boundary projection.
//!
//! The operations follow MMG's local remeshing loop (split long edges, collapse short
//! ones, swap for quality, relocate vertices), but boundary vertices are projected onto
//! the implicit field instead of a reconstructed surface. Boundary sizes come from the
//! chord sag of the exact surface: a triangle of side `l` whose vertex normals differ by
//! `theta` sags by about `l * theta / 6`. Interior sizes grade away from the boundary.
//!
//! Periodic axes are adapted on the torus: copies of a vertex on opposite box faces are
//! merged, geometry uses minimal-image displacements, and the periodic box planes become
//! internal constrained faces. Unfolding the result recreates identical opposite faces.
use crate::{Mesh, Point, cross, determinant, dot, geometry::Geometry, norm, quality, sub};
use std::collections::HashMap;

#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Maximum sampled distance between boundary triangles and the implicit walls.
    pub tolerance: f64,
    pub min_length: f64,
    pub max_length: f64,
    /// Maximum growth of the target length per unit distance (MMG's hgrad - 1).
    pub gradation: f64,
    /// Sheet bands: at least this many element lengths across the local thickness.
    /// Zero disables the thickness limit.
    pub layers: f64,
    pub iterations: usize,
    /// Smallest MMG quality an operation may create, unless the worst affected
    /// element was already worse.
    pub quality_guard: f64,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            tolerance: 0.01,
            min_length: 1e-3,
            max_length: 0.25,
            gradation: 0.3,
            layers: 0.,
            iterations: 4,
            quality_guard: 0.3,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct Stats {
    pub splits: usize,
    pub collapses: usize,
    pub swaps: usize,
    pub moves: usize,
    pub error_splits: usize,
    pub maximum_error: f64,
    /// Wall-clock seconds per phase: sizing, split, collapse, edge removal,
    /// face swap, relocation, error enforcement and unfolding.
    pub seconds: [f64; 7],
}

const SPLIT: f64 = 4. / 3.;
const COLLAPSE: f64 = 4. / 5.;
const SIZING_SAFETY: f64 = 0.8;
const WALL: u8 = 1;

type Face = [usize; 3];

/// Total order for finite sizes in a heap.
#[derive(Clone, Copy)]
struct Ordered(f64);
impl PartialEq for Ordered {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}
impl Eq for Ordered {}
impl PartialOrd for Ordered {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Ordered {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}

fn sorted(mut f: Face) -> Face {
    f.sort_unstable();
    f
}

fn unit(a: Point) -> Point {
    let n = norm(a);
    if n > 0. { a.map(|v| v / n) } else { a }
}

fn angle(a: Point, b: Point) -> f64 {
    norm(cross(a, b)).atan2(dot(a, b))
}

/// MMG-normalized quality: 1 for the regular tetrahedron, 0 when flat or inverted.
fn mmg(p: [Point; 4]) -> f64 {
    if determinant(p) <= 0. {
        return 0.;
    }
    quality(p).powf(1.5)
}

struct Tets<'g, G: Geometry> {
    g: &'g G,
    threshold: f64,
    points: Vec<Point>,
    /// Bit 0: implicit wall. Bits 2 + 2a + side: box plane.
    masks: Vec<u8>,
    features: Vec<u64>,
    /// Outward wall normal for wall vertices, zero otherwise.
    normals: Vec<Point>,
    frozen: Vec<bool>,
    tets: Vec<[usize; 4]>,
    alive: Vec<bool>,
    live_tetrahedra: usize,
    vt: Vec<Vec<usize>>,
    /// Boundary faces (sorted) and their label: 0 wall, 2.. box plane.
    faces: HashMap<Face, u8>,
    constraint_faces: Option<HashMap<Face, u64>>,
    vb: Vec<Vec<Face>>,
    /// Box length along periodic axes.
    period: [Option<f64>; 3],
    /// Wall vertices inserted short of the surface, still to be projected.
    pending: Vec<usize>,
}

impl<'g, G: Geometry> Tets<'g, G> {
    fn new(mesh: &Mesh, g: &'g G, threshold: f64, periodic: [bool; 3]) -> Self {
        Self::with_constraints(mesh, g, threshold, periodic, None)
    }

    fn with_constraints(
        mesh: &Mesh,
        g: &'g G,
        threshold: f64,
        periodic: [bool; 3],
        constraints: Option<&[crate::implicit::intersection::ConstraintMask]>,
    ) -> Self {
        let [lo, hi] = g.bounds();
        let period = std::array::from_fn(|a| periodic[a].then(|| hi[a] - lo[a]));
        // Coordinate buckets identify periodic copies, not distinct vertices of a
        // nonperiodic seed. Even an ulp-wide valid cell must retain its topology.
        let (points, map): (Vec<Point>, Vec<usize>) = if periodic.iter().any(|&axis| axis) {
            let mut index = HashMap::new();
            let mut points = Vec::new();
            let map: Vec<usize> = mesh
                .points
                .iter()
                .map(|&p| {
                    let c: Point = std::array::from_fn(|a| {
                        if periodic[a] && p[a] == hi[a] {
                            lo[a]
                        } else {
                            p[a]
                        }
                    });
                    let key: [i64; 3] = std::array::from_fn(|a| {
                        ((c[a] - lo[a]) / (hi[a] - lo[a]) * 1e11).round() as i64
                    });
                    *index.entry(key).or_insert_with(|| {
                        points.push(c);
                        points.len() - 1
                    })
                })
                .collect();
            (points, map)
        } else {
            (mesh.points.clone(), (0..mesh.points.len()).collect())
        };
        let n = points.len();
        let tets: Vec<[usize; 4]> = mesh.tets.iter().map(|t| t.map(|i| map[i])).collect();
        let mut s = Self {
            g,
            threshold,
            points,
            masks: vec![0; n],
            features: if let Some(masks) = constraints {
                masks.iter().map(|&mask| u64::from(mask)).collect()
            } else if g.halfspaces().is_some() {
                vec![0; n]
            } else {
                Vec::new()
            },
            normals: vec![[0.; 3]; n],
            frozen: vec![false; n],
            alive: vec![true; tets.len()],
            live_tetrahedra: tets.len(),
            tets,
            vt: vec![Vec::new(); n],
            faces: HashMap::new(),
            constraint_faces: constraints.map(|_| HashMap::new()),
            vb: vec![Vec::new(); n],
            period,
            pending: Vec::new(),
        };
        for (t, tet) in s.tets.iter().enumerate() {
            for &v in tet {
                s.vt[v].push(t);
            }
        }
        for &f in &mesh.surface {
            let p = f.map(|i| mesh.points[i]);
            let mut label = 0;
            for a in 0..3 {
                if p.iter().all(|q| q[a] == lo[a]) {
                    label = 2 + 2 * a as u8;
                } else if p.iter().all(|q| q[a] == hi[a]) {
                    // Opposite periodic faces are one internal plane on the torus.
                    label = if periodic[a] { 2 } else { 3 } + 2 * a as u8;
                }
            }
            let f = f.map(|i| map[i]);
            if !s.faces.contains_key(&sorted(f)) {
                let features = f.iter().fold(u64::MAX, |mask, &v| mask & s.feature_mask(v));
                s.add_feature_face(f, label, features);
            }
        }
        for v in 0..n {
            s.refresh(v);
        }
        if constraints.is_some() {
            s.pending = (0..n).filter(|&v| s.masks[v] & WALL != 0).collect();
        }
        s
    }

    /// Minimal-image displacement from `a` to `b`.
    fn delta(&self, a: Point, b: Point) -> Point {
        std::array::from_fn(|k| {
            let d = b[k] - a[k];
            self.period[k].map_or(d, |l| d - l * (d / l).round())
        })
    }

    /// Bring a point into the box along periodic axes.
    fn wrap(&self, mut p: Point) -> Point {
        let [lo, _] = self.g.bounds();
        for a in 0..3 {
            if let Some(l) = self.period[a] {
                p[a] = lo[a] + (p[a] - lo[a]).rem_euclid(l);
            }
        }
        p
    }

    /// Unwrapped coordinates of a simplex, relative to its first vertex.
    fn local<const N: usize>(&self, raw: [Point; N]) -> [Point; N] {
        let base = raw[0];
        raw.map(|p| {
            let d = self.delta(base, p);
            std::array::from_fn(|k| base[k] + d[k])
        })
    }

    fn add_feature_face(&mut self, f: Face, label: u8, features: u64) {
        let key = sorted(f);
        if let Some(faces) = &mut self.constraint_faces {
            faces.insert(key, features);
        }
        self.faces.insert(key, label);
        for v in key {
            self.vb[v].push(key);
            self.masks[v] |= if label == 0 { WALL } else { 1 << label };
        }
    }

    fn remove_face(&mut self, key: Face) -> Option<u8> {
        let label = self.faces.remove(&key)?;
        if let Some(faces) = &mut self.constraint_faces {
            faces.remove(&key);
        }
        for v in key {
            self.vb[v].retain(|&f| f != key);
        }
        Some(label)
    }

    /// Recompute a vertex's constraint mask, normal and frozen state.
    fn refresh(&mut self, v: usize) {
        let mut mask = 0;
        for f in &self.vb[v] {
            let label = self.faces[f];
            mask |= if label == 0 { WALL } else { 1 << label };
        }
        self.masks[v] = mask;
        if let Some(faces) = &self.constraint_faces {
            self.features[v] = self.vb[v].iter().fold(0, |mask, f| mask | faces[f]);
        }
        if let Some(planes) = self.g.halfspaces() {
            let tolerance = self.feature_tolerance();
            self.features[v] = if self.g.halfspace_regions().is_some() {
                self.vb[v]
                    .iter()
                    .filter(|face| self.faces[*face] == 0)
                    .fold(0, |mask, face| {
                        mask | face.iter().fold(u64::MAX, |common, &w| {
                            common & crate::polyhedral::active(planes, self.points[w], tolerance)
                        })
                    })
            } else {
                crate::polyhedral::active(planes, self.points[v], tolerance)
            };
        }
        self.normals[v] = if mask & WALL != 0 {
            self.outward(self.points[v])
        } else {
            [0.; 3]
        };
        // Points fixed by three constraints (box corners, a wall crossing a box edge)
        // never move; curves and box edges slide along their line.
        self.frozen[v] = if let Some(planes) = self.g.halfspaces()
            && self.feature_mask(v) != 0
        {
            let locked = std::array::from_fn(|a| (mask >> (2 + 2 * a)) & 3 != 0);
            crate::polyhedral::constraint_rank(planes, self.feature_mask(v), locked) == 3
        } else {
            mask.count_ones() >= 3 || self.feature_mask(v).count_ones() >= 2
        };
    }

    fn feature_mask(&self, v: usize) -> u64 {
        self.features.get(v).copied().unwrap_or(0)
    }

    fn feature_tolerance(&self) -> f64 {
        let [lo, hi] = self.g.bounds();
        (0..3).map(|a| hi[a] - lo[a]).fold(0., f64::max) * 1e-11
    }

    fn outward(&self, p: Point) -> Point {
        let v = self.g.value(p);
        let sign = if self.g.band() && v < 0. { -1. } else { 1. };
        unit(self.g.gradient(p).map(|x| sign * x))
    }

    fn level(&self, p: Point) -> f64 {
        self.g.level(self.g.value(p), self.threshold)
    }

    /// Newton projection onto the wall through `p`, holding box coordinates in `mask`.
    fn project(&self, mut p: Point, mask: u8, level: f64, features: u64) -> Option<Point> {
        let [lo, hi] = self.g.bounds();
        let locked: [bool; 3] = std::array::from_fn(|a| (mask >> (2 + 2 * a)) & 3 != 0);
        for a in 0..3 {
            if mask & (1 << (2 + 2 * a)) != 0 {
                p[a] = lo[a];
            }
            if mask & (1 << (3 + 2 * a)) != 0 {
                p[a] = hi[a];
            }
        }
        if let Some(fields) = self.g.constraints() {
            let p =
                crate::projection::onto_constraints(fields, features, p, locked, self.g.bounds())?;
            return (self.g.contains(p)
                && fields.iter().all(|field| {
                    let residual = field.value(p);
                    residual <= 0.
                        || residual
                            <= crate::norm(self.constraint_gradient(*field, p))
                                * self.feature_tolerance()
                }))
            .then_some(p);
        }
        if let Some(planes) = self.g.halfspaces() {
            let p = crate::polyhedral::project(
                planes,
                p,
                features,
                locked,
                self.feature_tolerance(),
                self.g.halfspace_regions(),
            )?;
            return self.g.contains(p).then_some(p);
        }
        if mask & WALL == 0 {
            return self.g.contains(p).then_some(p);
        }
        let scale = (0..3).map(|a| hi[a] - lo[a]).fold(0., f64::max);
        let p = crate::projection::onto_level(self.g, p, level, locked, 0.05 * scale)?;
        self.g.contains(p).then_some(p)
    }

    fn tet_points(&self, t: [usize; 4], moved: Option<(usize, Point)>) -> [Point; 4] {
        self.local(t.map(|i| match moved {
            Some((j, q)) if i == j => q,
            _ => self.points[i],
        }))
    }

    fn shell(&self, u: usize, v: usize) -> Vec<usize> {
        self.vt[u]
            .iter()
            .copied()
            .filter(|&t| self.tets[t].contains(&v))
            .collect()
    }

    fn adjacent(&self, u: usize) -> Vec<usize> {
        let mut out: Vec<usize> = self.vt[u]
            .iter()
            .flat_map(|&t| self.tets[t])
            .filter(|&w| w != u)
            .collect();
        out.sort_unstable();
        out.dedup();
        out
    }

    fn edges(&self) -> Vec<[usize; 2]> {
        let mut out = Vec::with_capacity(self.tets.len() * 2);
        for (t, tet) in self.tets.iter().enumerate() {
            if !self.alive[t] {
                continue;
            }
            for i in 0..4 {
                for j in i + 1..4 {
                    let (a, b) = (tet[i].min(tet[j]), tet[i].max(tet[j]));
                    out.push([a, b]);
                }
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    fn length(&self, u: usize, v: usize) -> f64 {
        norm(self.delta(self.points[u], self.points[v]))
    }

    /// Boundary faces containing edge `uv`.
    fn edge_faces(&self, u: usize, v: usize) -> Vec<Face> {
        self.vb[u]
            .iter()
            .copied()
            .filter(|f| f.contains(&v))
            .collect()
    }

    fn edge_mask(&self, faces: &[Face]) -> u8 {
        faces.iter().fold(0, |m, f| {
            let label = self.faces[f];
            m | if label == 0 { WALL } else { 1 << label }
        })
    }

    /// A canonical cut-plane vertex can represent the opposite endpoint image.
    /// Use that side's wall normal for the edge, without changing its node identity.
    fn edge_image_normal(&self, p: Point, toward: Point, normal: Point) -> Point {
        let [lo, hi] = self.g.bounds();
        let mut image = p;
        for a in 0..3 {
            if let Some(period) = self.period[a]
                && p[a] == lo[a]
            {
                let difference = toward[a] - p[a];
                if difference - period * (difference / period).round() < 0. {
                    image[a] = hi[a];
                }
            }
        }
        if image == p || normal == [0.; 3] {
            normal
        } else {
            self.outward(image)
        }
    }

    fn edge_turn(&self, p: Point, n: Point, q: Point, m: Point) -> f64 {
        angle(
            self.edge_image_normal(p, q, n),
            self.edge_image_normal(q, p, m),
        )
    }

    fn sag(&self, pu: Point, nu: Point, pv: Point, nv: Point) -> f64 {
        if norm(nu) == 0. || norm(nv) == 0. {
            return 0.;
        }
        norm(self.delta(pu, pv)) * self.edge_turn(pu, nu, pv, nv) / 6.
    }

    fn constraint_gradient(&self, field: &dyn crate::implicit::ScalarField, p: Point) -> Point {
        let [lo, hi] = self.g.bounds();
        let step = (0..3).map(|a| hi[a] - lo[a]).fold(0., f64::max) * 1e-5;
        field.gradient(p).unwrap_or_else(|| {
            std::array::from_fn(|a| {
                let mut lo = p;
                let mut hi = p;
                lo[a] -= step;
                hi[a] += step;
                (field.value(hi) - field.value(lo)) / (2. * step)
            })
        })
    }

    fn sizes(&self, o: &Options, edges: &[[usize; 2]]) -> Result<Vec<f64>, String> {
        let tol = o.tolerance * SIZING_SAFETY;
        let mut size = vec![o.max_length; self.points.len()];
        if let Some(targets) = self.g.target_lengths(&self.points) {
            if targets.len() != size.len()
                || targets
                    .iter()
                    .any(|value| !value.is_finite() || *value <= 0.)
            {
                return Err("size_field must return finite positive lengths".into());
            }
            for (length, target) in size.iter_mut().zip(targets) {
                *length = length.min(target);
            }
        }
        for &[u, v] in edges {
            // Flat patches and straight feature lines have zero curvature.
            // Endpoint subgradients across their junction are not a curvature metric.
            if self.g.halfspaces().is_some() {
                continue;
            }
            if self.masks[u] & self.masks[v] & WALL == 0 || self.edge_faces(u, v).is_empty() {
                continue;
            }
            let l = self.length(u, v);
            let theta = if let Some(fields) = self.g.constraints() {
                let common = self.feature_mask(u) & self.feature_mask(v);
                fields
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| common & (1u64 << j) != 0)
                    .map(|(_, field)| {
                        angle(
                            unit(self.constraint_gradient(*field, self.points[u])),
                            unit(self.constraint_gradient(*field, self.points[v])),
                        )
                    })
                    .fold(0., f64::max)
            } else {
                self.edge_turn(
                    self.points[u],
                    self.normals[u],
                    self.points[v],
                    self.normals[v],
                )
            };
            let mut target = (6. * tol * l / theta.max(1e-12)).sqrt();
            let locked = self.locked(u, v);
            if locked.iter().filter(|&&x| x).count() == 1 {
                // Curve chord: it sags with the curve's in-plane turning, which can far
                // exceed the wall's where the wall meets the plane at a shallow angle.
                let (t, s) = (self.curve_tangent(u, locked), self.curve_tangent(v, locked));
                let turn = norm(cross(t, s)).atan2(dot(t, s).abs());
                target = target.min((8. * tol * l / turn.max(1e-12)).sqrt());
            }
            for w in [u, v] {
                size[w] = size[w].min(target);
            }
        }
        if o.layers > 0. && self.g.band() {
            for ((s, &p), &mask) in size.iter_mut().zip(&self.points).zip(&self.masks) {
                if mask & WALL != 0 {
                    // Linear field across the sheet: thickness ~ band width / |grad|.
                    let g = norm(self.g.gradient(p));
                    if g > 0. {
                        *s = s.min(2. * self.threshold / g / o.layers);
                    }
                }
            }
        }
        for s in &mut size {
            *s = s.clamp(o.min_length, o.max_length);
        }
        // Lipschitz gradation in one Dijkstra pass from the smallest sizes.
        let mut adjacency = vec![Vec::new(); self.points.len()];
        for &[u, v] in edges {
            let l = self.length(u, v);
            adjacency[u].push((v, l));
            adjacency[v].push((u, l));
        }
        let mut heap: std::collections::BinaryHeap<_> = size
            .iter()
            .enumerate()
            .map(|(v, &h)| std::cmp::Reverse((Ordered(h), v)))
            .collect();
        while let Some(std::cmp::Reverse((Ordered(h), v))) = heap.pop() {
            if h > size[v] {
                continue;
            }
            for &(w, l) in &adjacency[v] {
                let limit = h + o.gradation * l;
                if limit < size[w] {
                    size[w] = limit;
                    heap.push(std::cmp::Reverse((Ordered(limit), w)));
                }
            }
        }
        Ok(size)
    }

    fn worst(&self, tets: &[usize]) -> f64 {
        tets.iter()
            .map(|&t| mmg(self.tet_points(self.tets[t], None)))
            .fold(1., f64::min)
    }

    fn add_tet(&mut self, t: [usize; 4]) -> usize {
        self.live_tetrahedra += 1;
        let i = self.tets.len();
        self.tets.push(t);
        self.alive.push(true);
        for v in t {
            self.vt[v].push(i);
        }
        i
    }

    fn kill_tet(&mut self, t: usize) {
        self.live_tetrahedra -= 1;
        self.alive[t] = false;
        for v in self.tets[t] {
            self.vt[v].retain(|&s| s != t);
        }
    }

    fn add_vertex(&mut self, p: Point) -> usize {
        self.points.push(p);
        self.masks.push(0);
        if !self.features.is_empty() {
            self.features.push(0);
        }
        self.normals.push([0.; 3]);
        self.frozen.push(false);
        self.vt.push(Vec::new());
        self.vb.push(Vec::new());
        self.points.len() - 1
    }

    fn split(&mut self, u: usize, v: usize, guard: f64, during_sizing: bool) -> Option<usize> {
        let shell = self.shell(u, v);
        if shell.is_empty() {
            return None;
        }
        let bfaces = self.edge_faces(u, v);
        let mask = self.edge_mask(&bfaces);
        let d = self.delta(self.points[u], self.points[v]);
        // Midpoint children are naturally worse than their parent; only reject slivers.
        let floor = (0.5 * self.worst(&shell).min(guard))
            .max(self.g.minimum_quality())
            .max(1e-9);
        let level = if mask & WALL != 0 {
            self.level(self.points[u])
        } else {
            0.
        };
        let mut chosen = None;
        // Early sizing is followed by a full quality stage. Late error repair
        // only uses off-center insertion on the constrained wall/box curves.
        let positions: &[f64] = if during_sizing || (mask & WALL != 0 && mask & !WALL != 0) {
            &[0.5, 0.25, 0.75]
        } else {
            &[0.5]
        };
        'positions: for &along in positions {
            let mid = self.wrap(std::array::from_fn(|a| self.points[u][a] + along * d[a]));
            let Some(target) = self.project(
                mid,
                mask,
                level,
                self.feature_mask(u) & self.feature_mask(v),
            ) else {
                continue;
            };
            // A deep concave sag can invert neighbours at full projection;
            // insert part of the way and finish projection after local repair.
            let offset = self.delta(mid, target);
            for fraction in [1., 0.75, 0.5, 0.25] {
                let p = self.wrap(std::array::from_fn(|a| mid[a] + fraction * offset[a]));
                let mut planned = Vec::new();
                let valid = shell.iter().all(|&t| {
                    let tet = self.tets[t];
                    let a = tet.map(|i| if i == v { usize::MAX } else { i });
                    let b = tet.map(|i| if i == u { usize::MAX } else { i });
                    planned.push((t, a, b));
                    [a, b].iter().all(|c| {
                        mmg(self.local(c.map(|i| if i == usize::MAX { p } else { self.points[i] })))
                            >= floor
                    })
                });
                if valid {
                    chosen = Some((p, fraction, planned));
                    break 'positions;
                }
            }
        }
        let (p, fraction, planned) = chosen?;
        let m = self.add_vertex(p);
        if !self.features.is_empty() {
            self.features[m] = self.feature_mask(u) & self.feature_mask(v);
        }
        if fraction < 1. {
            self.pending.push(m);
        }
        for (t, a, b) in planned {
            self.kill_tet(t);
            self.add_tet(a.map(|i| if i == usize::MAX { m } else { i }));
            self.add_tet(b.map(|i| if i == usize::MAX { m } else { i }));
        }
        for f in bfaces {
            let features = self.constraint_faces.as_ref().map_or(0, |faces| faces[&f]);
            let label = self.remove_face(f).unwrap();
            let c = f.into_iter().find(|&w| w != u && w != v).unwrap();
            self.add_feature_face([u, m, c], label, features);
            self.add_feature_face([m, v, c], label, features);
        }
        for w in [u, v, m] {
            self.refresh(w);
        }
        Some(m)
    }

    /// Apply the existing guarded short-edge collapses, shortest first.
    fn coarsen(&mut self, edges: &[[usize; 2]], size: &[f64], o: &Options) -> usize {
        let mut count = 0;
        let mut short: Vec<(f64, [usize; 2])> = edges
            .iter()
            .filter_map(|&[u, v]| {
                let r = self.length(u, v) / (0.5 * (size[u] + size[v]));
                (r < COLLAPSE).then_some((r, [u, v]))
            })
            .collect();
        short.sort_by(|a, b| a.0.total_cmp(&b.0));
        for (_, [u, v]) in short {
            if self.vt[u].is_empty() || self.vt[v].is_empty() || self.shell(u, v).is_empty() {
                continue;
            }
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

    /// Merge `u` into `v`.
    fn collapse(&mut self, u: usize, v: usize, size: &[f64], o: &Options) -> bool {
        let (mu, mv) = (self.masks[u], self.masks[v]);
        if self.feature_mask(u) & !self.feature_mask(v) != 0 {
            return false;
        }
        // Coincident cuts may leave separate copies of the same corner a few
        // ulps apart. Collapse those through the normal link and quality checks;
        // freezing both would preserve a roundoff-width sliver forever.
        let [lo, hi] = self.g.bounds();
        let roundoff = 64. * f64::EPSILON * (0..3).map(|a| hi[a] - lo[a]).fold(0., f64::max);
        let same_corner = (self.g.halfspaces().is_some() || self.g.constraints().is_some())
            && self.feature_mask(u) == self.feature_mask(v)
            && self.feature_mask(u).count_ones() >= 2
            && self.length(u, v) <= roundoff;
        if (self.frozen[u] && !same_corner) || mu & !mv != 0 || mu.count_ones() >= 3 {
            return false;
        }
        let bfaces = self.edge_faces(u, v);
        if mu != 0 && self.edge_mask(&bfaces) != mu {
            // A boundary vertex may only slide along its own patch or curve.
            return false;
        }
        if mu == 0 && !bfaces.is_empty() {
            return false;
        }
        let shell = self.shell(u, v);
        if shell.is_empty() {
            return false;
        }
        // Link condition, vertex level: common neighbours are exactly the shell ring.
        let au = self.adjacent(u);
        let av = self.adjacent(v);
        let mut ring: Vec<usize> = shell
            .iter()
            .flat_map(|&t| self.tets[t])
            .filter(|&w| w != u && w != v)
            .collect();
        ring.sort_unstable();
        ring.dedup();
        let common: Vec<usize> = au
            .iter()
            .copied()
            .filter(|w| av.binary_search(w).is_ok())
            .collect();
        if common != ring {
            return false;
        }
        if mu != 0 {
            // Surface link condition on the boundary triangulation.
            let bu: Vec<usize> = self.vb[u]
                .iter()
                .flatten()
                .copied()
                .filter(|&w| w != u)
                .collect();
            let bv: Vec<usize> = self.vb[v]
                .iter()
                .flatten()
                .copied()
                .filter(|&w| w != v)
                .collect();
            let mut common: Vec<usize> = bu.iter().copied().filter(|w| bv.contains(w)).collect();
            common.sort_unstable();
            common.dedup();
            let mut opposite: Vec<usize> = bfaces
                .iter()
                .flatten()
                .copied()
                .filter(|&w| w != u && w != v)
                .collect();
            opposite.sort_unstable();
            if common != opposite {
                return false;
            }
        }
        for &w in &au {
            if w != v && self.length(v, w) > SPLIT * 0.5 * (size[v] + size[w]) {
                return false;
            }
        }
        let others: Vec<usize> = self.vt[u]
            .iter()
            .copied()
            .filter(|t| !shell.contains(t))
            .collect();
        let before = self.worst(&self.vt[u].clone()).min(o.quality_guard);
        for &t in &others {
            let q = self.tet_points(self.tets[t], Some((u, self.points[v])));
            if mmg(q) < before {
                return false;
            }
        }
        if mu & WALL != 0 && self.g.halfspaces().is_none() {
            for f in &self.vb[u] {
                if f.contains(&v) || self.faces[f] != 0 {
                    continue;
                }
                let [a, b] = [f[0], f[1], f[2]]
                    .into_iter()
                    .filter(|&w| w != u)
                    .collect::<Vec<_>>()[..]
                else {
                    unreachable!()
                };
                for w in [a, b] {
                    let old_error = self.edge_error(
                        self.points[u],
                        self.normals[u],
                        self.points[w],
                        self.normals[w],
                        self.feature_mask(u) & self.feature_mask(w),
                    );
                    if self.masks[w] & WALL != 0
                        && self.edge_error(
                            self.points[v],
                            self.normals[v],
                            self.points[w],
                            self.normals[w],
                            self.feature_mask(v) & self.feature_mask(w),
                        ) > (o.tolerance * SIZING_SAFETY).max(old_error)
                    {
                        return false;
                    }
                }
            }
        }
        for t in shell {
            self.kill_tet(t);
        }
        for t in std::mem::take(&mut self.vt[u]) {
            self.tets[t] = self.tets[t].map(|i| if i == u { v } else { i });
            self.vt[v].push(t);
        }
        for f in bfaces {
            self.remove_face(f);
        }
        for f in std::mem::take(&mut self.vb[u]) {
            let features = self.constraint_faces.as_ref().map_or(0, |faces| faces[&f]);
            let label = self.remove_face(f).unwrap_or(0);
            self.add_feature_face(f.map(|i| if i == u { v } else { i }), label, features);
        }
        self.masks[u] = 0;
        for w in au {
            self.refresh(w);
        }
        self.refresh(v);
        true
    }

    /// Ordered ring of an interior edge's shell; None for boundary edges or odd shells.
    fn ring(&self, u: usize, v: usize, shell: &[usize]) -> Option<Vec<usize>> {
        let mut pairs = Vec::new();
        for &t in shell {
            let rest: Vec<usize> = self.tets[t]
                .into_iter()
                .filter(|&w| w != u && w != v)
                .collect();
            pairs.push([rest[0], rest[1]]);
        }
        let mut ring = vec![pairs[0][0], pairs[0][1]];
        let mut used = vec![false; pairs.len()];
        used[0] = true;
        while ring.len() < pairs.len() {
            let last = *ring.last().unwrap();
            let next = (0..pairs.len()).find(|&i| !used[i] && pairs[i].contains(&last))?;
            used[next] = true;
            let w = if pairs[next][0] == last {
                pairs[next][1]
            } else {
                pairs[next][0]
            };
            ring.push(w);
        }
        let closing = (0..pairs.len()).find(|&i| !used[i])?;
        if !(pairs[closing].contains(&ring[0]) && pairs[closing].contains(ring.last().unwrap())) {
            return None;
        }
        Some(ring)
    }

    /// Replace an interior edge shell by the best triangulation of its ring.
    fn remove_edge(&mut self, u: usize, v: usize) -> bool {
        if !self.edge_faces(u, v).is_empty() {
            return false;
        }
        let shell = self.shell(u, v);
        let n = shell.len();
        if !(3..=7).contains(&n) {
            return false;
        }
        let Some(ring) = self.ring(u, v, &shell) else {
            return false;
        };
        let before = self.worst(&shell);
        // Each ring triangle joins u and v, which must lie on opposite sides of it.
        let oriented = |a: usize, b: usize, c: usize| -> Option<[[usize; 4]; 2]> {
            let du = determinant(self.tet_points([a, b, c, u], None));
            let dv = determinant(self.tet_points([a, b, c, v], None));
            if du * dv >= 0. {
                return None;
            }
            let tu = if du > 0. { [a, b, c, u] } else { [a, c, b, u] };
            let tv = if dv > 0. { [a, b, c, v] } else { [a, c, b, v] };
            Some([tu, tv])
        };
        let q = |i: usize, k: usize, j: usize| -> f64 {
            oriented(ring[i], ring[k], ring[j]).map_or(f64::NEG_INFINITY, |ts| {
                ts.iter()
                    .map(|&t| mmg(self.tet_points(t, None)))
                    .fold(1., f64::min)
            })
        };
        // Dynamic programming over polygon triangulations, maximizing the worst quality.
        let mut best = vec![vec![f64::NEG_INFINITY; n]; n];
        let mut split = vec![vec![usize::MAX; n]; n];
        for i in 0..n - 1 {
            best[i][i + 1] = f64::INFINITY;
        }
        for gap in 2..n {
            for i in 0..n - gap {
                let j = i + gap;
                for k in i + 1..j {
                    let value = best[i][k].min(best[k][j]).min(q(i, k, j));
                    if value > best[i][j] {
                        best[i][j] = value;
                        split[i][j] = k;
                    }
                }
            }
        }
        if best[0][n - 1] <= before + 1e-6 {
            return false;
        }
        let mut new = Vec::new();
        let mut stack = vec![(0, n - 1)];
        while let Some((i, j)) = stack.pop() {
            if j - i < 2 {
                continue;
            }
            let k = split[i][j];
            let Some(ts) = oriented(ring[i], ring[k], ring[j]) else {
                return false;
            };
            new.extend(ts);
            stack.push((i, k));
            stack.push((k, j));
        }
        let volume = |ts: &mut dyn Iterator<Item = [usize; 4]>| -> f64 {
            ts.map(|t| determinant(self.tet_points(t, None)).abs())
                .sum()
        };
        let old = volume(&mut shell.iter().map(|&t| self.tets[t]));
        if (volume(&mut new.iter().copied()) - old).abs() > 1e-9 * old {
            return false;
        }
        for t in shell {
            self.kill_tet(t);
        }
        for t in new {
            self.add_tet(t);
        }
        true
    }

    /// 2-3 swap across an interior face shared by two tetrahedra.
    fn swap_face(&mut self, t: usize, i: usize) -> bool {
        let tet = self.tets[t];
        let d = tet[i];
        let face: Vec<usize> = tet.into_iter().filter(|&w| w != d).collect();
        if self
            .faces
            .contains_key(&sorted([face[0], face[1], face[2]]))
        {
            return false;
        }
        let Some(&s) = self.vt[face[0]]
            .iter()
            .find(|&&s| s != t && face.iter().all(|w| self.tets[s].contains(w)))
        else {
            return false;
        };
        let Some(e) = self.tets[s].into_iter().find(|w| !face.contains(w)) else {
            return false;
        };
        if self.adjacent(d).contains(&e) {
            return false;
        }
        let before = self.worst(&[t, s]);
        let [a, b, c] = [face[0], face[1], face[2]];
        let mut new = Vec::new();
        for [x, y] in [[a, b], [b, c], [c, a]] {
            let mut tet = [x, y, d, e];
            if determinant(self.tet_points(tet, None)) < 0. {
                tet.swap(0, 1);
            }
            new.push(tet);
        }
        let after = new
            .iter()
            .map(|&t| mmg(self.tet_points(t, None)))
            .fold(1., f64::min);
        // The three new tetrahedra must exactly fill the old pair.
        let volume = |ts: &[[usize; 4]]| -> f64 {
            ts.iter()
                .map(|&t| determinant(self.tet_points(t, None)).abs())
                .sum()
        };
        let old = volume(&[self.tets[t], self.tets[s]]);
        if after <= before + 1e-6 || (volume(&new) - old).abs() > 1e-9 * old {
            return false;
        }
        self.kill_tet(t);
        self.kill_tet(s);
        for tet in new {
            self.add_tet(tet);
        }
        true
    }

    fn relocate(&mut self, v: usize, o: &Options) -> bool {
        let mask = self.masks[v];
        if self.frozen[v] || self.vt[v].is_empty() {
            return false;
        }
        let p = self.points[v];
        let adj = self.adjacent(v);
        let mean: Point = adj.iter().fold([0.; 3], |m, &w| {
            let d = self.delta(p, self.points[w]);
            std::array::from_fn(|a| m[a] + d[a] / adj.len() as f64)
        });
        let centroid: Point = std::array::from_fn(|a| p[a] + mean[a]);
        let mut target = std::array::from_fn(|a| p[a] + 0.5 * (centroid[a] - p[a]));
        if mask & WALL != 0 {
            let n = self.normals[v];
            let d = sub(target, p);
            target = std::array::from_fn(|a| p[a] + d[a] - dot(d, n) * n[a]);
        }
        let level = if mask & WALL != 0 { self.level(p) } else { 0. };
        let Some(q) = self.project(target, mask, level, self.feature_mask(v)) else {
            return false;
        };
        let incident = self.vt[v].clone();
        let before = self.worst(&incident);
        let after = incident
            .iter()
            .map(|&t| mmg(self.tet_points(self.tets[t], Some((v, q)))))
            .fold(1., f64::min);
        if after <= before.min(o.quality_guard) || after <= 0. {
            return false;
        }
        if mask & WALL != 0 && self.g.halfspaces().is_none() {
            let n = self.outward(q);
            for f in &self.vb[v] {
                if self.faces[f] != 0 {
                    continue;
                }
                for &w in f {
                    if w != v
                        && self.masks[w] & WALL != 0
                        && self.edge_error(
                            q,
                            n,
                            self.points[w],
                            self.normals[w],
                            self.feature_mask(v) & self.feature_mask(w),
                        ) > (o.tolerance * SIZING_SAFETY).max(self.edge_error(
                            p,
                            self.normals[v],
                            self.points[w],
                            self.normals[w],
                            self.feature_mask(v) & self.feature_mask(w),
                        ))
                    {
                        return false;
                    }
                }
            }
        }
        self.points[v] = q;
        self.refresh(v);
        true
    }

    fn edge_error(&self, p: Point, n: Point, q: Point, m: Point, features: u64) -> f64 {
        let [lo, hi] = self.g.bounds();
        let locked = std::array::from_fn(|a| {
            let p_cap = p[a] == lo[a] || p[a] == hi[a];
            let q_cap = q[a] == lo[a] || q[a] == hi[a];
            p_cap && q_cap && (p[a] == q[a] || self.period[a].is_some())
        });
        let delta = self.delta(p, q);
        let midpoint = std::array::from_fn(|a| p[a] + 0.5 * delta[a]);
        if self.g.constraints().is_none() {
            let sag = self.sag(p, n, q, m);
            // Near a tangent cap, a small 3D wall error can hide a much larger
            // error of its trace. Collapse and relocation must preserve that trace.
            return if locked.iter().any(|&axis| axis) {
                sag.max(self.distance(midpoint, locked))
            } else {
                sag
            };
        }
        self.constraint_distance(midpoint, locked, features)
    }

    fn constraint_distance(&self, p: Point, locked: [bool; 3], features: u64) -> f64 {
        let distance = self.distance(p, locked);
        self.g.constraints().map_or(distance, |fields| {
            crate::projection::onto_constraints(fields, features, p, locked, self.g.bounds())
                .map_or(f64::INFINITY, |q| distance.max(norm(sub(q, p))))
        })
    }

    /// First-order wall distance moving only along unlocked axes: a chord in a box
    /// plane must follow the wall/plane curve inside that plane.
    fn distance(&self, mut p: Point, locked: [bool; 3]) -> f64 {
        // Local simplices can extend outside the fundamental cell. Evaluate
        // their actual periodic image, preserving a join's one-sided normal.
        let [lo, hi] = self.g.bounds();
        for a in 0..3 {
            if let Some(period) = self.period[a]
                && (p[a] < lo[a] || p[a] > hi[a])
            {
                p[a] = lo[a] + (p[a] - lo[a]).rem_euclid(period);
            }
        }
        let r = self.g.residual(p, self.threshold).abs();
        if r == 0. {
            return 0.;
        }
        let mut g = self.g.gradient(p);
        for a in 0..3 {
            if locked[a] {
                g[a] = 0.;
            }
        }
        r / norm(g)
    }

    /// Box axes held by both ends of an edge.
    fn locked(&self, u: usize, v: usize) -> [bool; 3] {
        let shared = self.masks[u] & self.masks[v];
        std::array::from_fn(|a| (shared >> (2 + 2 * a)) & 3 != 0)
    }

    /// Unit tangent of the wall/plane curve through a curve vertex.
    fn curve_tangent(&self, v: usize, locked: [bool; 3]) -> Point {
        let a = (0..3).find(|&a| locked[a]).unwrap_or(0);
        let mut n = [0.; 3];
        n[a] = 1.;
        unit(cross(self.normals[v], n))
    }

    fn check_budget(&self, maximum: usize) -> Result<(), String> {
        if self.live_tetrahedra > maximum {
            return Err(format!(
                "adaptation exceeded max_tetrahedra: {} > {maximum}",
                self.live_tetrahedra
            ));
        }
        Ok(())
    }

    /// Move partially inserted wall vertices onto the surface where now valid.
    fn settle(&mut self, guard: f64, maximum: usize) -> Result<bool, String> {
        let mut changed = false;
        let pending = std::mem::take(&mut self.pending);
        for v in pending {
            if self.vt[v].is_empty() {
                continue;
            }
            let mask = self.masks[v];
            let level = self.level(self.points[v]);
            // Same floor as a split: geometry may not be bought with slivers.
            let floor = (0.5 * self.worst(&self.vt[v].clone()).min(guard))
                .max(self.g.minimum_quality())
                .max(1e-6);
            let Some(q) = self.project(self.points[v], mask, level, self.feature_mask(v)) else {
                self.pending.push(v);
                continue;
            };
            let blocking: Vec<usize> = self.vt[v]
                .iter()
                .copied()
                .filter(|&t| mmg(self.tet_points(self.tets[t], Some((v, q)))) < floor)
                .collect();
            if blocking.is_empty() {
                changed |= self.points[v] != q;
                self.points[v] = q;
                self.refresh(v);
                continue;
            }
            let origin = self.points[v];
            let residual = (self.g.value(origin) - level).abs();
            let partial = (1..=8).find_map(|halvings| {
                let fraction = 0.5f64.powi(halvings);
                let candidate =
                    std::array::from_fn(|axis| origin[axis] + fraction * (q[axis] - origin[axis]));
                if !self.g.contains(candidate) {
                    return None;
                }
                let value = self.g.value(candidate);
                (self.g.level(value, self.threshold) == level
                    && (value - level).abs() < 0.999 * residual
                    && self.vt[v].iter().all(|&tetrahedron| {
                        mmg(self.tet_points(self.tets[tetrahedron], Some((v, candidate)))) >= floor
                    }))
                .then_some(candidate)
            });
            if let Some(candidate) = partial {
                changed = true;
                self.points[v] = candidate;
                self.refresh(v);
                self.pending.push(v);
                continue;
            }
            // Make room: split the longest edge facing the vertex in each element the
            // move would flatten, then try again in the next round.
            let mut facing = Vec::new();
            for t in blocking {
                let rest: Vec<usize> = self.tets[t].into_iter().filter(|&w| w != v).collect();
                let mut best = (0., [rest[0], rest[1]]);
                for (i, j) in [(0, 1), (1, 2), (0, 2)] {
                    let l = self.length(rest[i], rest[j]);
                    if l > best.0 {
                        best = (l, [rest[i], rest[j]]);
                    }
                }
                facing.push(best.1);
            }
            for [a, b] in facing {
                if !self.shell(a, b).is_empty() {
                    changed |= self.split(a, b, guard, false).is_some();
                    self.check_budget(maximum)?;
                }
            }
            self.pending.push(v);
        }
        Ok(changed)
    }

    fn enforce(&mut self, o: &Options, stats: &mut Stats, maximum: usize) -> Result<(), String> {
        for _ in 0..16 {
            self.g.check()?;
            let settled = self.settle(o.quality_guard, maximum)?;
            let mut flagged = Vec::new();
            let mut worst = 0f64;
            let wall: Vec<Face> = self
                .faces
                .iter()
                .filter(|(_, l)| **l == 0)
                .map(|(f, _)| *f)
                .collect();
            for f in wall {
                let p = self.local(f.map(|i| self.points[i]));
                let c = std::array::from_fn(|a| (p[0][a] + p[1][a] + p[2][a]) / 3.);
                let features = self.constraint_faces.as_ref().map_or(0, |faces| faces[&f]);
                let dc = self.constraint_distance(c, [false; 3], features);
                worst = worst.max(dc);
                let mut any = false;
                let mut longest = (0., [f[0], f[1]]);
                for (i, j) in [(0, 1), (1, 2), (0, 2)] {
                    let m = std::array::from_fn(|a| 0.5 * (p[i][a] + p[j][a]));
                    let d = self.constraint_distance(
                        m,
                        self.locked(f[i], f[j]),
                        self.feature_mask(f[i]) & self.feature_mask(f[j]),
                    );
                    worst = worst.max(d);
                    let l = self.length(f[i], f[j]);
                    let edge = [f[i].min(f[j]), f[i].max(f[j])];
                    if l > longest.0 {
                        longest = (l, edge);
                    }
                    if d > o.tolerance {
                        flagged.push((d, l, edge));
                        any = true;
                    }
                }
                if dc > o.tolerance && !any {
                    flagged.push((dc, longest.0, longest.1));
                }
            }
            stats.maximum_error = worst;
            if flagged.is_empty() {
                return Ok(());
            }
            // Preserve the worst report for a shared edge, then repair the largest
            // defect before another split spends its cavity's quality margin.
            flagged.sort_by(|a, b| a.2.cmp(&b.2).then_with(|| b.0.total_cmp(&a.0)));
            flagged.dedup_by_key(|e| e.2);
            flagged.sort_by(|a, b| {
                b.0.total_cmp(&a.0)
                    .then_with(|| b.1.total_cmp(&a.1))
                    .then_with(|| a.2.cmp(&b.2))
            });
            let mut progress = settled;
            let mut blocked_vertices = Vec::new();
            for (_, _, [a, b]) in flagged {
                let shell = self.shell(a, b);
                if shell.is_empty() {
                    continue;
                }
                if self.length(a, b) >= o.min_length
                    && self.split(a, b, o.quality_guard, false).is_some()
                {
                    self.check_budget(maximum)?;
                    stats.error_splits += 1;
                    progress = true;
                } else {
                    for t in shell {
                        blocked_vertices.extend(self.tets[t]);
                    }
                }
            }
            if !blocked_vertices.is_empty() {
                // Error insertion may need room that the preceding quality pass
                // could not see before its new wall nodes existed.
                blocked_vertices.sort_unstable();
                blocked_vertices.dedup();
                for v in blocked_vertices {
                    if self.relocate(v, o) {
                        stats.moves += 1;
                        progress = true;
                    }
                }
            }
            if !progress {
                return Ok(());
            }
        }
        Ok(())
    }

    /// Unfold the torus into the box: each tetrahedron moves to the cell holding its
    /// centroid, and plane vertices get one copy per side they are used on.
    fn into_mesh(self) -> Result<Mesh, String> {
        let [lo, hi] = self.g.bounds();
        let mut map = HashMap::<(usize, [u8; 3]), usize>::new();
        let mut points = Vec::new();
        let mut tets = Vec::new();
        for (t, tet) in self.tets.iter().enumerate() {
            if !self.alive[t] {
                continue;
            }
            let mut p = self.tet_points(*tet, None);
            for a in 0..3 {
                if let Some(l) = self.period[a] {
                    let c = p.iter().map(|q| q[a]).sum::<f64>() / 4.;
                    let shift = ((c - lo[a]) / l).floor() * l;
                    for q in &mut p {
                        q[a] -= shift;
                    }
                }
            }
            let mut out = [0; 4];
            for (k, &v) in tet.iter().enumerate() {
                let mut image = [0u8; 3];
                let mut q = self.points[v];
                for a in 0..3 {
                    if let Some(l) = self.period[a] {
                        let m = ((p[k][a] - q[a]) / l).round();
                        if m == 1. && q[a] == lo[a] {
                            image[a] = 1;
                            q[a] = hi[a];
                        } else if m != 0. {
                            return Err(format!(
                                "a tetrahedron crosses a periodic plane: tet={t}, vertex={v}, axis={a}, image={m}, canonical={q:?}, local={:?}, mask={}, tet_points={:?}",
                                p[k],
                                self.masks[v],
                                tet.map(|vertex| self.points[vertex])
                            ));
                        }
                    }
                }
                out[k] = *map.entry((v, image)).or_insert_with(|| {
                    points.push(q);
                    points.len() - 1
                });
            }
            tets.push(out);
        }
        let surface = crate::boundary(&tets)?;
        Ok(Mesh {
            points,
            tets,
            surface,
        })
    }
}

/// Adapt a valid boundary-conforming mesh towards curvature and thickness sizes.
pub(crate) fn adapt(
    mesh: &mut Mesh,
    g: &impl Geometry,
    threshold: f64,
    periodic: [bool; 3],
    o: &Options,
    maximum: usize,
) -> Result<Stats, String> {
    adapt_with_polish(mesh, g, threshold, periodic, o, maximum, 0, None)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn adapt_with_polish(
    mesh: &mut Mesh,
    g: &impl Geometry,
    threshold: f64,
    periodic: [bool; 3],
    o: &Options,
    maximum: usize,
    polish: usize,
    constraints: Option<&[crate::implicit::intersection::ConstraintMask]>,
) -> Result<Stats, String> {
    if !(o.tolerance > 0. && o.min_length > 0. && o.max_length >= o.min_length)
        || o.gradation.is_nan()
        || o.gradation <= 0.
        || !(0.0..1.).contains(&o.quality_guard)
        || o.layers.is_nan()
        || o.layers < 0.
        || o.iterations > 100
    {
        return Err("invalid adaptation options".into());
    }
    let [lo, hi] = g.bounds();
    if (0..3).any(|a| periodic[a] && o.max_length > 0.25 * (hi[a] - lo[a])) {
        return Err("periodic adaptation needs max_length <= a quarter of the period".into());
    }
    if constraints
        .is_some_and(|masks| masks.len() != mesh.points.len() || periodic.iter().any(|&a| a))
    {
        return Err("retained constraints require one mask per nonperiodic vertex".into());
    }
    let mut s = if constraints.is_some() {
        Tets::with_constraints(mesh, g, threshold, periodic, constraints)
    } else {
        Tets::new(mesh, g, threshold, periodic)
    };
    s.check_budget(maximum)?;
    let mut stats = Stats::default();
    let mut clock = std::time::Instant::now();
    let mut lap = |stats: &mut Stats, phase: usize| {
        stats.seconds[phase] += clock.elapsed().as_secs_f64();
        clock = std::time::Instant::now();
    };
    for _ in 0..o.iterations {
        g.check()?;
        let edges = s.edges();
        let size = s.sizes(o, &edges)?;
        lap(&mut stats, 0);
        let mut long: Vec<(f64, [usize; 2])> = edges
            .iter()
            .filter_map(|&[u, v]| {
                let r = s.length(u, v) / (0.5 * (size[u] + size[v]));
                (r > SPLIT).then_some((r, [u, v]))
            })
            .collect();
        long.sort_by(|a, b| b.0.total_cmp(&a.0));
        let mut reclaimed_budget = false;
        for (_, [u, v]) in long {
            if !reclaimed_budget && s.live_tetrahedra.saturating_add(s.shell(u, v).len()) > maximum
            {
                // A split adds one cell per incident tet. Reclaim overresolved
                // cells before failing the cap, using the normal size, link,
                // feature, geometry and quality guards. Try once per pass.
                lap(&mut stats, 1);
                let edges = s.edges();
                let size = s.sizes(o, &edges)?;
                lap(&mut stats, 0);
                stats.collapses += s.coarsen(&edges, &size, o);
                lap(&mut stats, 2);
                reclaimed_budget = true;
            }
            if s.split(u, v, o.quality_guard, true).is_some() {
                s.check_budget(maximum)?;
                stats.splits += 1;
            }
        }
        lap(&mut stats, 1);
        let edges = s.edges();
        let size = s.sizes(o, &edges)?;
        lap(&mut stats, 0);
        stats.collapses += s.coarsen(&edges, &size, o);
        lap(&mut stats, 2);
        // Edge removal only pays off around poor elements.
        let mut candidates = Vec::new();
        for t in 0..s.tets.len() {
            if s.alive[t] && mmg(s.tet_points(s.tets[t], None)) < 0.4 {
                let tet = s.tets[t];
                for i in 0..4 {
                    for j in i + 1..4 {
                        candidates.push([tet[i].min(tet[j]), tet[i].max(tet[j])]);
                    }
                }
            }
        }
        candidates.sort_unstable();
        candidates.dedup();
        for [u, v] in candidates {
            if s.remove_edge(u, v) {
                s.check_budget(maximum)?;
                stats.swaps += 1;
            }
        }
        lap(&mut stats, 3);
        for t in 0..s.tets.len() {
            if !s.alive[t] || mmg(s.tet_points(s.tets[t], None)) > 0.5 {
                continue;
            }
            for i in 0..4 {
                if s.alive[t] && s.swap_face(t, i) {
                    s.check_budget(maximum)?;
                    stats.swaps += 1;
                    break;
                }
            }
        }
        lap(&mut stats, 4);
        s.settle(o.quality_guard, maximum)?;
        for v in 0..s.points.len() {
            if s.relocate(v, o) {
                stats.moves += 1;
            }
        }
        lap(&mut stats, 5);
    }
    s.enforce(o, &mut stats, maximum)?;
    s.settle(o.quality_guard, maximum)?;
    lap(&mut stats, 6);
    for _ in 0..polish {
        g.check()?;
        for v in 0..s.points.len() {
            if s.relocate(v, o) {
                stats.moves += 1;
            }
        }
    }
    lap(&mut stats, 5);
    *mesh = s.into_mesh()?;
    lap(&mut stats, 6);
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use crate::implicit::{Options, Region, generate};
    use crate::{Point, dot};

    #[derive(serde::Deserialize)]
    struct CapStar {
        points: Vec<Point>,
        tets: Vec<[usize; 4]>,
        surface: Vec<[usize; 3]>,
        bad_coordinates: [Point; 2],
        coefficients: Vec<[f64; 5]>,
        threshold: f64,
    }

    impl crate::geometry::Geometry for CapStar {
        fn value(&self, p: Point) -> f64 {
            self.coefficients
                .iter()
                .fold(self.threshold, |value, &[a, b, c, amplitude, phase]| {
                    value
                        - amplitude
                            * (std::f64::consts::TAU
                                * (a * (p[0] + 0.5) + b * (p[1] + 0.5) + c * (p[2] + 0.5))
                                + phase)
                                .cos()
                })
        }
        fn gradient(&self, p: Point) -> Point {
            self.coefficients
                .iter()
                .fold([0.; 3], |mut gradient, &[a, b, c, amplitude, phase]| {
                    let factor = amplitude
                        * std::f64::consts::TAU
                        * (std::f64::consts::TAU
                            * (a * (p[0] + 0.5) + b * (p[1] + 0.5) + c * (p[2] + 0.5))
                            + phase)
                            .sin();
                    for (value, mode) in gradient.iter_mut().zip([a, b, c]) {
                        *value += factor * mode;
                    }
                    gradient
                })
        }
        fn bounds(&self) -> [Point; 2] {
            [[-0.5; 3], [0.5; 3]]
        }
        fn rank(&self, p: Point) -> [i64; 3] {
            p.map(|x| (x * 1e11).round() as i64)
        }
        fn band(&self) -> bool {
            false
        }
        fn minimum_quality(&self) -> f64 {
            0.1
        }
        fn geometry_tolerance(&self) -> Option<f64> {
            Some(0.002)
        }
    }

    #[test]
    fn periodic_cap_curve_split_and_enforcement() {
        use crate::geometry::Geometry;
        let fixture: CapStar =
            serde_json::from_str(include_str!("testdata/spinodoid-cap-star.json")).unwrap();
        let mesh = crate::Mesh {
            points: fixture.points.clone(),
            tets: fixture.tets.clone(),
            surface: fixture.surface.clone(),
        };
        let mut state = super::Tets::new(&mesh, &fixture, 0., [true; 3]);
        let [u, v] = fixture
            .bad_coordinates
            .map(|point| state.points.iter().position(|p| *p == point).unwrap());
        let mid = std::array::from_fn(|a| 0.5 * (state.points[u][a] + state.points[v][a]));
        let faces = state.edge_faces(u, v);
        let mask = state.edge_mask(&faces);
        let target = state.project(mid, mask, 0., 0);
        assert_eq!(state.shell(u, v).len(), 3);
        assert!(state.distance(mid, [true, false, false]) > 0.002);
        assert!(target.is_some(), "the exact cap trace has a nearby root");
        assert!(fixture.value(target.unwrap()).abs() < 1e-10);
        let inserted = state.split(u, v, 0.3, false);
        assert!(
            inserted.is_some(),
            "a quarter projection keeps every child above quality 0.1"
        );
        let options = super::Options {
            tolerance: 0.001,
            min_length: 0.002,
            max_length: 0.1,
            ..Default::default()
        };
        let mut stats = super::Stats::default();
        state.enforce(&options, &mut stats, 5000).unwrap();
        assert!(
            stats.maximum_error <= fixture.geometry_tolerance().unwrap(),
            "periodic cap curve must pass the fixed geometry gate with its internal aim unchanged"
        );
        assert!(state.pending.is_empty());
        assert!(
            state
                .tets
                .iter()
                .enumerate()
                .filter(|(t, _)| state.alive[*t])
                .all(|(_, &tet)| super::mmg(state.tet_points(tet, None)) >= 0.1)
        );
        let mut rebuilt = super::Tets::new(&mesh, &fixture, 0., [true; 3]);
        let mut rebuilt_stats = super::Stats::default();
        rebuilt.enforce(&options, &mut rebuilt_stats, 5000).unwrap();
        assert!(
            rebuilt_stats.maximum_error <= fixture.geometry_tolerance().unwrap(),
            "the post-quality star must repair without a manually inserted midpoint"
        );
        assert!(rebuilt.pending.is_empty());
        assert!(
            rebuilt
                .tets
                .iter()
                .enumerate()
                .filter(|(t, _)| rebuilt.alive[*t])
                .all(|(_, &tet)| super::mmg(rebuilt.tet_points(tet, None)) >= 0.1)
        );
    }

    #[test]
    fn periodic_spinodoid_closes_geometry_after_quality() {
        struct Field<'a>(&'a CapStar);
        impl crate::implicit::ScalarField for Field<'_> {
            fn value(&self, p: Point) -> f64 {
                crate::geometry::Geometry::value(self.0, p)
            }
            fn gradient(&self, p: Point) -> Option<Point> {
                Some(crate::geometry::Geometry::gradient(self.0, p))
            }
        }
        let fixture: CapStar =
            serde_json::from_str(include_str!("testdata/spinodoid-cap-star.json")).unwrap();
        let output = generate(
            &Field(&fixture),
            Options {
                bounds: [[-0.5; 3], [0.5; 3]],
                cells: [24; 3],
                periodic: [true; 3],
                geometry_tolerance: 0.002,
                minimum_quality: 0.1,
                max_tetrahedra: 250_000,
                optimize_passes: 1,
                adapt: Some(super::Options {
                    tolerance: 0.001,
                    min_length: 0.002,
                    max_length: 0.1,
                    ..Default::default()
                }),
                ..Default::default()
            },
        )
        .unwrap();
        assert!(output.diagnostics.maximum_sampled_surface_error <= 0.002);
        assert!(output.diagnostics.minimum_mmg_quality >= 0.1);
        assert!(output.mesh.points.len() <= 250_000 && output.mesh.tets.len() <= 250_000);
        assert!(output.periodic_pairs.iter().all(|pairs| !pairs.is_empty()));
    }

    struct PlaneWave;

    #[test]
    fn cap_curve_edge_error_respects_the_held_plane() {
        struct ShallowWall;
        impl crate::geometry::Geometry for ShallowWall {
            fn value(&self, p: Point) -> f64 {
                (std::f64::consts::TAU * (p[0] + 0.5)).sin() + 0.01 * (p[2] - 10. * p[1] * p[1])
            }
            fn gradient(&self, p: Point) -> Point {
                [
                    std::f64::consts::TAU * (std::f64::consts::TAU * (p[0] + 0.5)).cos(),
                    -0.2 * p[1],
                    0.01,
                ]
            }
            fn bounds(&self) -> [Point; 2] {
                [[-0.5; 3], [0.5; 3]]
            }
            fn rank(&self, _: Point) -> [i64; 3] {
                [0; 3]
            }
            fn band(&self) -> bool {
                false
            }
        }
        let mesh = crate::Mesh {
            points: Vec::new(),
            tets: Vec::new(),
            surface: Vec::new(),
        };
        let state = super::Tets::new(&mesh, &ShallowWall, 0., [true, false, false]);
        for x in [-0.5, 0.5] {
            let p = [x, -0.025, 0.00625];
            let q = [x, 0.025, 0.00625];
            let (n, m) = (state.outward(p), state.outward(q));
            assert!(
                state.sag(p, n, q, m) < 0.0001,
                "3D normal turning understates the held-plane chord"
            );
            assert!(
                state.edge_error(p, n, q, m, 0) > 0.006,
                "cap-curve geometry guards must measure error inside the cap plane"
            );
        }
    }

    impl crate::geometry::Geometry for PlaneWave {
        fn value(&self, point: Point) -> f64 {
            (std::f64::consts::TAU * (point[1] + 0.5)).sin() + 0.001
        }
        fn gradient(&self, point: Point) -> Point {
            [
                0.,
                std::f64::consts::TAU * (std::f64::consts::TAU * (point[1] + 0.5)).cos(),
                0.,
            ]
        }
        fn bounds(&self) -> [Point; 2] {
            [[-0.5; 3], [0.5; 3]]
        }
        fn rank(&self, _point: Point) -> [i64; 3] {
            [0; 3]
        }
        fn band(&self) -> bool {
            false
        }
    }

    #[test]
    fn affine_curve_mobility_counts_independent_constraints() {
        struct Column;
        impl crate::geometry::Geometry for Column {
            fn value(&self, p: Point) -> f64 {
                p[0].abs().max(p[1].abs()) - 0.2
            }
            fn gradient(&self, _: Point) -> Point {
                [1., 0., 0.]
            }
            fn bounds(&self) -> [Point; 2] {
                [[-0.5; 3], [0.5; 3]]
            }
            fn rank(&self, _: Point) -> [i64; 3] {
                [0; 3]
            }
            fn band(&self) -> bool {
                false
            }
            fn halfspaces(&self) -> Option<&[[f64; 4]]> {
                Some(&[
                    [1., 0., 0., -0.2],
                    [-1., 0., 0., -0.2],
                    [0., 1., 0., -0.2],
                    [0., -1., 0., -0.2],
                ])
            }
        }
        let mesh = crate::Mesh {
            points: vec![
                [0.2, 0.2, 0.],
                [0.2, 0.2, 0.5],
                [0.2, 0., 0.5],
                [0., 0.2, 0.5],
            ],
            tets: Vec::new(),
            surface: vec![[0, 1, 2], [0, 1, 3], [1, 2, 3]],
        };
        let tets = super::Tets::new(&mesh, &Column, 0., [false; 3]);
        assert!(!tets.frozen[0], "two planes leave a movable line");
        assert!(tets.frozen[1], "the cap fixes the line endpoint");
        assert!(!tets.frozen[2] && !tets.frozen[3]);
        let moved = tets
            .project([0.21, 0.21, 0.1], tets.masks[0], 0., tets.features[0])
            .unwrap();
        assert_eq!(moved, [0.2, 0.2, 0.1]);
        let held = tets
            .project([0.21, 0.21, 0.1], tets.masks[1], 0., tets.features[1])
            .unwrap();
        assert_eq!(held, [0.2, 0.2, 0.5]);
    }

    #[test]
    fn periodic_distance_uses_the_correct_image_of_a_normal_jump() {
        struct C0Wall;
        impl crate::geometry::Geometry for C0Wall {
            fn value(&self, p: Point) -> f64 {
                p[0] + 0.3 * p[2].abs() - 0.35
            }
            fn gradient(&self, p: Point) -> Point {
                [1., 0., if p[2] < 0. { -0.3 } else { 0.3 }]
            }
            fn bounds(&self) -> [Point; 2] {
                [[-0.5; 3], [0.5; 3]]
            }
            fn rank(&self, _: Point) -> [i64; 3] {
                [0; 3]
            }
            fn band(&self) -> bool {
                false
            }
        }
        let mesh = crate::Mesh {
            points: Vec::new(),
            tets: Vec::new(),
            surface: Vec::new(),
        };
        let field = C0Wall;
        let periodic = super::Tets::new(&mesh, &field, 0., [false, false, true]);
        let nonperiodic = super::Tets::new(&mesh, &field, 0., [false; 3]);
        let x = 0.35 - 0.3 * 0.49;
        for z in [-0.51, 0.51] {
            assert!(periodic.distance([x, 0., z], [false; 3]) < 1e-12);
            assert!(periodic.distance([x, 0., z], [false, false, true]) < 1e-12);
            assert!(nonperiodic.distance([x, 0., z], [false; 3]) > 1e-3);
        }
        assert!(periodic.distance([x, 0., 0.49], [false; 3]) < 1e-12);
        // The canonical lower-cap node represents an upper-cap endpoint in
        // this edge's physical image. Its straight wall has zero turning.
        let join = [0.2, 0., -0.5];
        let upper_neighbor = [x, 0., 0.49];
        assert!(
            periodic.sag(
                join,
                periodic.outward(join),
                upper_neighbor,
                periodic.outward(upper_neighbor),
            ) < 1e-12
        );
        assert!(
            periodic.sag(
                upper_neighbor,
                periodic.outward(upper_neighbor),
                join,
                periodic.outward(join),
            ) < 1e-12
        );
        let lower = [0.347, 0., -0.01];
        let upper = [0.347, 0., 0.01];
        assert!(
            periodic.sag(
                lower,
                periodic.outward(lower),
                upper,
                periodic.outward(upper)
            ) > 1e-4
        );
        // A single genuine wall patch is enough to check the sizing path.
        let patch = crate::Mesh {
            points: vec![
                [0.2, 0., 0.5],
                upper_neighbor,
                [0.2, 0.01, 0.5],
                [0.19, 0., 0.49],
            ],
            tets: vec![[0, 2, 1, 3]],
            surface: vec![[0, 1, 2]],
        };
        let adapter = super::Tets::new(&patch, &field, 0., [false, false, true]);
        let sizes = adapter
            .sizes(
                &super::Options {
                    tolerance: 1e-5,
                    min_length: 1e-5,
                    max_length: 0.1,
                    ..Default::default()
                },
                &[[0, 1]],
            )
            .unwrap();
        assert_eq!(sizes[0], 0.1);
        assert_eq!(sizes[1], 0.1);
    }

    #[test]
    fn nonperiodic_seed_keeps_distinct_nearby_vertices() {
        let mesh = crate::Mesh {
            points: vec![[0., 0., 0.], [1e-13, 0., 0.], [0., 0.1, 0.], [0., 0., 0.1]],
            tets: vec![[0, 1, 2, 3]],
            surface: Vec::new(),
        };
        let tets = super::Tets::new(&mesh, &PlaneWave, 0., [false; 3]);
        assert_eq!(tets.points, mesh.points);
        assert_eq!(tets.tets, mesh.tets);
        assert!(crate::determinant(tets.tet_points(tets.tets[0], None)) > 0.);
    }

    #[test]
    fn projection_cannot_move_a_wall_vertex_across_a_periodic_plane() {
        let mesh = crate::Mesh {
            points: Vec::new(),
            tets: Vec::new(),
            surface: Vec::new(),
        };
        let tets = super::Tets::new(&mesh, &PlaneWave, 0., [true; 3]);
        assert!(tets.project([-0.4, -0.5001, 0.23], 0, 0., 0).is_none());
        assert!(
            tets.project([-0.4, -0.4999, 0.23], super::WALL, 0., 0)
                .is_none()
        );
        let point = tets
            .project([-0.4, 0.4999, 0.23], super::WALL, 0., 0)
            .unwrap();
        assert!(point[1] > 0.49 && point[1] < 0.5);
    }

    #[test]
    fn pending_wall_projection_makes_safe_partial_progress() {
        struct PendingPlane;

        impl crate::geometry::Geometry for PendingPlane {
            fn value(&self, point: Point) -> f64 {
                -point[2]
            }
            fn gradient(&self, _point: Point) -> Point {
                [0., 0., -1.]
            }
            fn bounds(&self) -> [Point; 2] {
                [[-1.; 3], [1.; 3]]
            }
            fn rank(&self, _point: Point) -> [i64; 3] {
                [0; 3]
            }
            fn band(&self) -> bool {
                false
            }
            fn minimum_quality(&self) -> f64 {
                0.1
            }
        }

        let points = vec![
            [0., 0., 0.2],
            [0.2, 0., 0.05],
            [-0.1, 0.03f64.sqrt(), 0.05],
            [-0.1, -0.03f64.sqrt(), 0.05],
        ];
        let mut tetrahedron = [0, 1, 2, 3];
        if crate::determinant(tetrahedron.map(|node| points[node])) < 0. {
            tetrahedron.swap(1, 2);
        }
        let mesh = crate::Mesh {
            points,
            tets: vec![tetrahedron],
            surface: vec![[0, 1, 2]],
        };
        let mut state = super::Tets::new(&mesh, &PendingPlane, 0., [false; 3]);
        state.pending.push(0);
        state.settle(0.3, 1000).unwrap();
        assert!(
            state.points[0][2] < 0.2,
            "pending wall vertex did not make safe partial progress"
        );
        assert!(state.points[0][2] > 0.);
        assert!(state.pending.contains(&0));
        assert_eq!(state.live_tetrahedra, 1);
        for (tetrahedron, alive) in state.tets.iter().zip(&state.alive) {
            if *alive {
                assert!(super::mmg(state.tet_points(*tetrahedron, None)) >= 0.1);
            }
        }
    }

    #[test]
    fn graded_sizing_reclaims_short_edges_before_exhausting_budget() {
        use crate::geometry::Geometry;
        struct GradedBox {
            fine: f64,
        }
        impl Geometry for GradedBox {
            fn value(&self, _: Point) -> f64 {
                -1.
            }
            fn gradient(&self, _: Point) -> Point {
                [0.; 3]
            }
            fn bounds(&self) -> [Point; 2] {
                [[0.; 3], [1.; 3]]
            }
            fn rank(&self, p: Point) -> [i64; 3] {
                p.map(|x| (x * 1e11).round() as i64)
            }
            fn band(&self) -> bool {
                false
            }
            fn minimum_quality(&self) -> f64 {
                0.1
            }
            fn target_lengths(&self, p: &[Point]) -> Option<Vec<f64>> {
                Some(
                    p.iter()
                        .map(|p| self.fine + (0.6 - self.fine) * p[0])
                        .collect(),
                )
            }
        }
        let mut seed = generate(
            &|_: Point| -1.,
            Options {
                cells: [4; 3],
                optimize_passes: 0,
                minimum_quality: 0.1,
                ..Default::default()
            },
        )
        .unwrap()
        .mesh;
        let maximum = seed.tets.len();
        let o = super::Options {
            min_length: 0.1,
            max_length: 0.6,
            gradation: 0.5,
            iterations: 1,
            quality_guard: 0.1,
            ..Default::default()
        };
        let stats = super::adapt(
            &mut seed,
            &GradedBox { fine: 0.2 },
            0.,
            [false; 3],
            &o,
            maximum,
        )
        .expect("coarsenable bulk must not block refinement under the same cap");
        assert!(stats.splits > 0 && stats.collapses > 0);
        assert!(seed.tets.len() <= maximum);
        let mut volume = 0.;
        for tet in &seed.tets {
            let points = tet.map(|i| seed.points[i]);
            assert!(super::mmg(points) >= 0.1);
            volume += crate::determinant(points) / 6.;
        }
        assert!((volume - 1.).abs() < 1e-12);
        assert!(
            seed.points
                .iter()
                .flatten()
                .all(|&x| (0.0..=1.0).contains(&x))
        );
        println!(
            "cap={maximum}, admitted={}, splits={}, collapses={}",
            seed.tets.len(),
            stats.splits,
            stats.collapses
        );
        let mut insufficient = generate(
            &|_: Point| -1.,
            Options {
                cells: [4; 3],
                optimize_passes: 0,
                minimum_quality: 0.1,
                ..Default::default()
            },
        )
        .unwrap()
        .mesh;
        let error = super::adapt(
            &mut insufficient,
            &GradedBox { fine: 0.1 },
            0.,
            [false; 3],
            &o,
            maximum,
        )
        .expect_err("a genuinely insufficient cap must still fail");
        assert!(error.contains("exceeded max_tetrahedra"), "{error}");
    }

    fn adapt(tolerance: f64) -> Option<super::Options> {
        Some(super::Options {
            tolerance,
            min_length: tolerance,
            max_length: 0.1,
            ..Default::default()
        })
    }

    #[test]
    fn sphere_adapts_with_fewer_tetrahedra_and_exact_walls() {
        let field = |p: Point| dot(p, p).sqrt() - 0.35;
        let base = Options {
            bounds: [[-0.5; 3], [0.5; 3]],
            cells: [24; 3],
            geometry_tolerance: 1e-3,
            optimize_passes: 2,
            minimum_quality: 0.1,
            ..Default::default()
        };
        let lattice = generate(&field, base).unwrap();
        let adapted = generate(
            &field,
            Options {
                adapt: adapt(1e-3),
                ..base
            },
        )
        .unwrap();
        let exact = 4. / 3. * std::f64::consts::PI * 0.35f64.powi(3);
        let d = &adapted.diagnostics;
        assert!(d.maximum_sampled_surface_error <= 1e-3);
        assert!((d.volume / exact - 1.).abs() < 0.01, "{}", d.volume / exact);
        assert!(
            adapted.mesh.tets.len() < lattice.mesh.tets.len(),
            "{} vs {} ({:?})",
            adapted.mesh.tets.len(),
            lattice.mesh.tets.len(),
            d.adapt
        );
        assert!(d.minimum_mmg_quality >= 0.1, "{}", d.minimum_mmg_quality);
    }

    #[test]
    fn periodic_gyroid_band_stays_periodic() {
        let field = |p: Point| {
            let t = std::f64::consts::TAU;
            (t * p[0]).sin() * (t * p[1]).cos()
                + (t * p[1]).sin() * (t * p[2]).cos()
                + (t * p[2]).sin() * (t * p[0]).cos()
        };
        let out = generate(
            &field,
            Options {
                bounds: [[-0.5; 3], [0.5; 3]],
                cells: [24; 3],
                region: Region::Band {
                    lower: -0.55,
                    upper: 0.55,
                },
                periodic: [true; 3],
                geometry_tolerance: 2e-3,
                optimize_passes: 2,
                minimum_quality: 0.1,
                adapt: adapt(2e-3),
                ..Default::default()
            },
        )
        .unwrap();
        assert!(out.periodic_pairs.iter().all(|p| !p.is_empty()));
        assert!(out.diagnostics.maximum_sampled_surface_error <= 2e-3);
    }
}
