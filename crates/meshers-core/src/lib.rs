//! Experimental direct, lattice-based meshing of a smooth gyroid band.
//! This is not an implementation of the certified isosurface-stuffing stencils.
//! CPU-only dependencies; export and accelerator implementations live in adapters.
use serde::Serialize;
use std::collections::{BTreeMap, HashMap};
use std::f64::consts::TAU;

pub mod accelerator;
mod error;
mod geometry;
pub mod implicit;
pub use error::MeshingError;
use geometry::Geometry;
mod profile;

/// Available CPU parallelism, capped at the supported maximum of 256 workers.
/// Respects operating-system availability where Rust can determine it.
#[must_use]
pub fn available_threads() -> usize {
    std::thread::available_parallelism()
        .map_or(1, std::num::NonZeroUsize::get)
        .min(256)
}

pub type Point = [f64; 3];
pub type Tet = [usize; 4];
pub fn sub(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}
pub fn dot(a: Point, b: Point) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
pub fn cross(a: Point, b: Point) -> Point {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub fn norm(a: Point) -> f64 {
    dot(a, a).sqrt()
}
fn lerp(a: Point, b: Point, t: f64) -> Point {
    std::array::from_fn(|i| a[i] + t * (b[i] - a[i]))
}
fn canonical(p: Point) -> Point {
    p.map(|x| if x == 1. { 0. } else { x })
}
pub fn field(p: Point) -> f64 {
    let [x, y, z] = canonical(p).map(|x| TAU * x);
    x.sin() * y.cos() + y.sin() * z.cos() + z.sin() * x.cos()
}
fn gradient(p: Point) -> Point {
    let [x, y, z] = canonical(p).map(|x| TAU * x);
    [
        TAU * (x.cos() * y.cos() - z.sin() * x.sin()),
        TAU * (y.cos() * z.cos() - x.sin() * y.sin()),
        TAU * (z.cos() * x.cos() - y.sin() * z.sin()),
    ]
}
fn rank(p: Point) -> [i64; 3] {
    canonical(p).map(|x| (x * 1e11).round() as i64)
}
pub fn determinant(p: [Point; 4]) -> f64 {
    dot(sub(p[1], p[0]), cross(sub(p[2], p[0]), sub(p[3], p[0])))
}
pub fn quality(p: [Point; 4]) -> f64 {
    let l2: f64 = (0..4)
        .flat_map(|i| (i + 1..4).map(move |j| dot(sub(p[i], p[j]), sub(p[i], p[j]))))
        .sum();
    12. * (determinant(p).abs() / 2.).powf(2. / 3.) / l2
}
pub fn min_dihedral(p: [Point; 4]) -> f64 {
    let normals: Vec<Point> = (0..4)
        .map(|i| {
            let f: Vec<usize> = (0..4).filter(|&j| j != i).collect();
            let mut n = cross(sub(p[f[1]], p[f[0]]), sub(p[f[2]], p[f[0]]));
            if dot(n, sub(p[i], p[f[0]])) < 0. {
                n = n.map(|x| -x);
            }
            let len = norm(n);
            n.map(|x| x / len)
        })
        .collect();
    (0..4)
        .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
        .map(|(i, j)| {
            (-dot(normals[i], normals[j]))
                .clamp(-1., 1.)
                .acos()
                .to_degrees()
        })
        .fold(180., f64::min)
}

#[derive(Clone, Copy)]
pub struct GyroidOptions {
    pub optimize_passes: usize,
    pub post_optimize_passes: usize,
    pub threads: usize,
    pub surface_weight: f64,
    pub surface_tolerance: f64,
    pub max_tetrahedra: usize,
    pub align_cube_edges: bool,
    pub resolution: usize,
    pub threshold: f64,
    pub snap: f64,
    pub periodic: bool,
}
impl Default for GyroidOptions {
    fn default() -> Self {
        Self {
            optimize_passes: 0,
            post_optimize_passes: 4,
            threads: 1,
            surface_weight: 0.,
            surface_tolerance: 0.,
            max_tetrahedra: usize::MAX,
            align_cube_edges: false,
            resolution: 40,
            threshold: 0.5,
            snap: 0.2,
            periodic: true,
        }
    }
}
pub struct Mesh {
    pub points: Vec<Point>,
    pub tets: Vec<Tet>,
    pub surface: Vec<[usize; 3]>,
}
#[derive(Serialize)]
pub struct Metrics {
    pub points: usize,
    pub tetrahedra: usize,
    pub surface_triangles: usize,
    pub volume: f64,
    pub mean_ratio_min: f64,
    pub mean_ratio_p01: f64,
    pub mean_ratio_median: f64,
    pub mean_ratio_below_01: usize,
    pub mmg_quality_min: f64,
    pub mmg_quality_p01: f64,
    pub mmg_quality_median: f64,
    pub mmg_quality_below_01: usize,
    pub minimum_dihedral_degrees: f64,
    pub maximum_surface_node_residual: f64,
    pub periodic_nodes_and_triangles: [bool; 3],
}

struct Builder<'a, G: Geometry> {
    geometry: &'a G,
    points: Vec<Point>,
    values: Vec<f64>,
    edges: HashMap<(usize, usize, i8), usize>,
    threshold: f64,
}
impl<G: Geometry> Builder<'_, G> {
    fn add(&mut self, p: Point, value: f64) -> usize {
        let i = self.points.len();
        self.points.push(p);
        self.values.push(value);
        i
    }
    fn intersection(&mut self, a: usize, b: usize, side: i8) -> usize {
        let level = self.threshold * f64::from(side);
        if (self.values[a] - level).abs() < 1e-12 {
            return a;
        }
        if (self.values[b] - level).abs() < 1e-12 {
            return b;
        }
        let (a, b) = if a < b { (a, b) } else { (b, a) };
        let key = (a, b, side);
        if let Some(&i) = self.edges.get(&key) {
            return i;
        }
        let (pa, pb) = (self.points[a], self.points[b]);
        let (mut lo, mut hi) = (0., 1.);
        let (mut flo, mut fhi) = (self.values[a] - level, self.values[b] - level);
        let mut t = flo / (flo - fhi);
        for _ in 0..40 {
            let f = self.geometry.value(lerp(pa, pb, t)) - level;
            if f.abs() < 1e-13 {
                break;
            }
            if f.signum() == flo.signum() {
                lo = t;
                flo = f;
            } else {
                hi = t;
                fhi = f;
            }
            let secant = (lo * fhi - hi * flo) / (fhi - flo);
            t = if secant > lo + 0.05 * (hi - lo) && secant < hi - 0.05 * (hi - lo) {
                secant
            } else {
                (lo + hi) / 2.
            };
        }
        let i = self.add(lerp(pa, pb, t), level);
        self.edges.insert(key, i);
        i
    }
    fn clip(&mut self, faces: Vec<Vec<usize>>, side: i8) -> Vec<Vec<usize>> {
        let mut out = Vec::new();
        let mut cap_edges = Vec::new();
        for face in faces {
            let mut poly = Vec::new();
            let mut cuts = Vec::new();
            for k in 0..face.len() {
                let a = face[k];
                let b = face[(k + 1) % face.len()];
                let da = f64::from(side) * self.values[a] - self.threshold;
                let db = f64::from(side) * self.values[b] - self.threshold;
                let ia = da <= 1e-12;
                let ib = db <= 1e-12;
                if ia {
                    poly.push(a);
                }
                if ia != ib {
                    let p = self.intersection(a, b, side);
                    poly.push(p);
                    cuts.push(p);
                }
            }
            poly.dedup();
            if poly.len() > 1 && poly.first() == poly.last() {
                poly.pop();
            }
            if poly.len() >= 3 {
                out.push(poly);
            }
            cuts.sort_unstable();
            cuts.dedup();
            if cuts.len() == 2 {
                cap_edges.push((cuts[0], cuts[1]));
            }
        }
        // Recover the cap's cycle from face intersections, avoiding angle sorting.
        if cap_edges.len() >= 3 {
            let first = cap_edges[0].0;
            let mut cap = vec![first];
            let mut cur = first;
            let mut prev = usize::MAX;
            for _ in 0..cap_edges.len() {
                let next = cap_edges.iter().find_map(|&(a, b)| {
                    if a == cur && b != prev {
                        Some(b)
                    } else if b == cur && a != prev {
                        Some(a)
                    } else {
                        None
                    }
                });
                let Some(next) = next else { break };
                if next == first {
                    break;
                }
                cap.push(next);
                prev = cur;
                cur = next;
            }
            if cap.len() >= 3 {
                out.push(cap);
            }
        }
        out
    }
    fn triangulate(&self, face: &[usize]) -> Vec<[usize; 3]> {
        let pivot = (0..face.len())
            .min_by_key(|&i| self.geometry.rank(self.points[face[i]]))
            .unwrap();
        (1..face.len() - 1)
            .map(|k| {
                [
                    face[pivot],
                    face[(pivot + k) % face.len()],
                    face[(pivot + k + 1) % face.len()],
                ]
            })
            .collect()
    }
}

/// Generate a unit-box gyroid band with the CPU implementation.
///
/// This preserves the experimental geometry and option semantics. General fields
/// and physical domains are planned separately; see docs/geometry-contract.md.
pub fn generate_gyroid(o: GyroidOptions) -> Result<Mesh, MeshingError> {
    generate_gyroid_with_accelerator(o, None)
}

/// Experimental extension point used by the optional CUDA adapter.
/// The adapter is trusted to honor the accelerator contract.
pub fn generate_gyroid_with_accelerator(
    o: GyroidOptions,
    accelerator: Option<&dyn accelerator::Factory>,
) -> Result<Mesh, MeshingError> {
    validate_options(o).map_err(MeshingError::InvalidOptions)?;
    generate_impl(o, accelerator).map_err(MeshingError::GenerationFailed)
}

fn generate_impl(
    o: GyroidOptions,
    accelerator: Option<&dyn accelerator::Factory>,
) -> Result<Mesh, String> {
    let mut profile = profile::Profile::new("generate");
    // Resolve both sides of the band on each background edge.
    let n = o.resolution;
    let h = 1. / n as f64;
    if 2. * 3_f64.sqrt() * TAU * h > 2. * o.threshold {
        return Err("resolution too coarse for the band; increase --resolution".into());
    }
    let mut b = Builder {
        geometry: &geometry::Legacy,
        points: Vec::new(),
        values: Vec::new(),
        edges: HashMap::new(),
        threshold: o.threshold,
    };
    // Align coordinate planes with analytic gyroid intersections on cube edges.
    // A separable monotone map preserves translated opposite-face patterns.
    let a = o.threshold.asin() / TAU;
    let knots: Vec<_> = [0., a, 0.5 - a, 0.5 + a, 1. - a, 1.]
        .into_iter()
        .map(|x| ((x * n as f64).round() / n as f64, x))
        .collect();
    if o.align_cube_edges
        && knots
            .windows(2)
            .any(|w| w[1].0 <= w[0].0 || w[1].1 <= w[0].1)
    {
        return Err("cube-edge alignment requires distinct resolved edge intersections".into());
    }
    let warp = |x: f64| {
        if !o.align_cube_edges || x == 0. || x == 1. {
            return x;
        }
        let w = knots.windows(2).find(|w| x <= w[1].0).unwrap();
        w[0].1 + (x - w[0].0) * (w[1].1 - w[0].1) / (w[1].0 - w[0].0)
    };
    let id = |x: usize, y: usize, z: usize| x + (n + 1) * (y + (n + 1) * z);
    for z in 0..=n {
        for y in 0..=n {
            for x in 0..=n {
                let p = [x as f64 * h, y as f64 * h, z as f64 * h].map(warp);
                b.add(p, field(p));
            }
        }
    }
    let center_start = b.points.len();
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let p = [
                    (x as f64 + 0.5) * h,
                    (y as f64 + 0.5) * h,
                    (z as f64 + 0.5) * h,
                ]
                .map(warp);
                b.add(p, field(p));
            }
        }
    }
    profile.mark("lattice_and_field");
    // Constrained Newton snapping: cube coordinates are held fixed. Identical
    // periodic representatives use the same function, derivative and decisions.
    if o.snap > 0. {
        for i in 0..b.points.len() {
            let original = b.points[i];
            let value = b.values[i];
            let level = value.signum() * o.threshold;
            let mut g = gradient(original);
            for a in 0..3 {
                if original[a] == 0. || original[a] == 1. {
                    g[a] = 0.;
                }
            }
            if norm(g) < 1e-12 || (value - level).abs() / norm(g) > o.snap * h {
                continue;
            }
            let mut p = original;
            for _ in 0..8 {
                let f = field(p) - level;
                let mut g = gradient(p);
                for a in 0..3 {
                    if original[a] == 0. || original[a] == 1. {
                        g[a] = 0.;
                    }
                }
                let g2 = dot(g, g);
                if g2 < 1e-20 {
                    break;
                }
                p = std::array::from_fn(|a| p[a] - f * g[a] / g2);
            }
            if p.iter().all(|&v| (0.0..=1.0).contains(&v))
                && norm(sub(p, original)) <= o.snap * h
                && (field(p) - level).abs() < 1e-11
            {
                b.points[i] = p;
                b.values[i] = level;
            }
        }
    }
    profile.mark("snapping");
    let cube_faces = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [3, 7, 6, 2],
        [0, 4, 7, 3],
        [1, 2, 6, 5],
    ];
    let mut tets = Vec::new();
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let corners = [
                    id(x, y, z),
                    id(x + 1, y, z),
                    id(x + 1, y + 1, z),
                    id(x, y + 1, z),
                    id(x, y, z + 1),
                    id(x + 1, y, z + 1),
                    id(x + 1, y + 1, z + 1),
                    id(x, y + 1, z + 1),
                ];
                let c = center_start + x + n * (y + n * z);
                // True BCC interior: adjacent cell centers plus one shared
                // face edge form a disphenoid. Boundary pyramids end at the cube.
                for (fi, face) in cube_faces.iter().enumerate() {
                    let face = face.map(|i| corners[i]);
                    let coord = [z, z, y, y, x, x][fi];
                    let positive = fi % 2 == 1;
                    let boundary_face = if positive { coord == n - 1 } else { coord == 0 };
                    let mut background = Vec::with_capacity(4);
                    if boundary_face {
                        for tri in b.triangulate(&face) {
                            background.push([c, tri[0], tri[1], tri[2]]);
                        }
                    } else if positive {
                        let neighbor = c + [n * n, n * n, n, n, 1, 1][fi];
                        for j in 0..4 {
                            background.push([c, neighbor, face[j], face[(j + 1) % 4]]);
                        }
                    } else {
                        continue;
                    }
                    for tet in background {
                        if tet.iter().all(|&i| b.values[i] > o.threshold + 1e-12)
                            || tet.iter().all(|&i| b.values[i] < -o.threshold - 1e-12)
                        {
                            continue;
                        }
                        if determinant(tet.map(|i| b.points[i])).abs() < h * h * h * 0.01 {
                            return Err("snapping degenerated a background tetrahedron".into());
                        }
                        let faces = vec![
                            vec![tet[0], tet[1], tet[2]],
                            vec![tet[0], tet[3], tet[1]],
                            vec![tet[0], tet[2], tet[3]],
                            vec![tet[1], tet[3], tet[2]],
                        ];
                        let faces = b.clip(faces, 1);
                        let faces = b.clip(faces, -1);
                        let Some(anchor) = faces
                            .iter()
                            .flatten()
                            .copied()
                            .min_by_key(|&i| rank(b.points[i]))
                        else {
                            continue;
                        };
                        // Pulling triangulation uses one global ordering, including on
                        // shared and periodic faces. It needs no new interior vertices.
                        for f in faces {
                            if f.contains(&anchor) {
                                continue;
                            }
                            for tri in b.triangulate(&f) {
                                let mut t = [anchor, tri[0], tri[1], tri[2]];
                                let det = determinant(t.map(|i| b.points[i]));
                                if det.abs() < 1e-22 {
                                    return Err("degenerate generated tetrahedron; change resolution or snap".into());
                                }
                                if det < 0. {
                                    t.swap(2, 3);
                                }
                                tets.push(t);
                            }
                        }
                    }
                }
            }
        }
    }
    profile.mark("initial_tetrahedra");
    let mut map = vec![usize::MAX; b.points.len()];
    let mut points = Vec::new();
    for tet in &mut tets {
        for i in tet {
            if map[*i] == usize::MAX {
                map[*i] = points.len();
                points.push(b.points[*i]);
            }
            *i = map[*i];
        }
    }
    if tets.len() > o.max_tetrahedra {
        return Err("initial lattice exceeds the tetrahedron budget; reduce resolution".into());
    }
    profile.mark("point_compaction");
    let surface = boundary(&tets)?;
    profile.mark("initial_boundary");
    let mut mesh = Mesh {
        points,
        tets,
        surface,
    };
    optimize::optimize(
        &mut mesh,
        o.threshold,
        h,
        o.optimize_passes,
        o.surface_weight,
        o.threads,
        accelerator,
    )?;
    profile.mark("initial_optimization");
    adaptive::refine(
        &mut mesh,
        o.threshold,
        o.surface_tolerance,
        o.max_tetrahedra,
    )?;
    profile.mark("adaptive_refinement");
    if o.surface_tolerance > 0. {
        optimize::optimize(
            &mut mesh,
            o.threshold,
            h,
            o.post_optimize_passes,
            o.surface_weight,
            o.threads,
            accelerator,
        )?;
    }
    profile.mark("post_optimization");
    let periodic = mesh.periodicity();
    profile.mark("periodicity_validation");
    if o.periodic && !periodic.iter().all(|&x| x) {
        return Err(format!("periodic boundary check failed: {periodic:?}"));
    }
    Ok(mesh)
}

fn tet_faces(t: &Tet) -> [[usize; 3]; 4] {
    [
        [t[1], t[2], t[3]],
        [t[0], t[3], t[2]],
        [t[0], t[1], t[3]],
        [t[0], t[2], t[1]],
    ]
}
fn face_key(mut f: [usize; 3]) -> [usize; 3] {
    f.sort_unstable();
    f
}
fn face_parity(f: [usize; 3]) -> u8 {
    ((f[0] > f[1]) as u8 + (f[0] > f[2]) as u8 + (f[1] > f[2]) as u8) % 2
}
fn boundary(tets: &[Tet]) -> Result<Vec<[usize; 3]>, String> {
    let mut faces: Vec<_> = tets
        .iter()
        .flat_map(tet_faces)
        .map(|f| (face_key(f), f))
        .collect();
    faces.sort_unstable_by_key(|v| v.0);
    let mut surface = Vec::new();
    let mut i = 0;
    while i < faces.len() {
        let mut j = i + 1;
        while j < faces.len() && faces[j].0 == faces[i].0 {
            j += 1;
        }
        match j - i {
            1 => surface.push(faces[i].1),
            2 => {
                if face_parity(faces[i].1) == face_parity(faces[i + 1].1) {
                    return Err("inconsistent orientation across a shared face".into());
                }
            }
            _ => return Err("nonmanifold tetrahedral face".into()),
        }
        i = j;
    }
    surface.sort_unstable();
    Ok(surface)
}
impl Mesh {
    pub fn periodicity(&self) -> [bool; 3] {
        std::array::from_fn(|axis| {
            let sides: Vec<_> = [0., 1.]
                .iter()
                .map(|&value| {
                    let nodes: BTreeMap<_, _> = self
                        .points
                        .iter()
                        .enumerate()
                        .filter(|(_, p)| (p[axis] - value).abs() < 1e-10)
                        .map(|(i, &p)| (rank(p), i))
                        .collect();
                    let mut faces: Vec<_> = self
                        .surface
                        .iter()
                        .filter(|f| {
                            f.iter()
                                .all(|&i| (self.points[i][axis] - value).abs() < 1e-10)
                        })
                        .map(|f| {
                            let mut key = f.map(|i| rank(self.points[i]));
                            key.sort_unstable();
                            key
                        })
                        .collect();
                    faces.sort_unstable();
                    (nodes.keys().copied().collect::<Vec<_>>(), faces)
                })
                .collect();
            sides[0] == sides[1] && !sides[0].0.is_empty()
        })
    }
    pub fn metrics(&self, threshold: f64) -> Metrics {
        let mut q: Vec<_> = self
            .tets
            .iter()
            .map(|t| quality(t.map(|i| self.points[i])))
            .collect();
        q.sort_by(f64::total_cmp);
        let residual = self
            .surface
            .iter()
            .filter(|f| {
                !(0..3).any(|a| {
                    f.iter().all(|&i| self.points[i][a].abs() < 1e-10)
                        || f.iter().all(|&i| (self.points[i][a] - 1.).abs() < 1e-10)
                })
            })
            .flat_map(|f| f.iter())
            .map(|&i| (field(self.points[i]).abs() - threshold).abs())
            .fold(0., f64::max);
        Metrics {
            points: self.points.len(),
            tetrahedra: self.tets.len(),
            surface_triangles: self.surface.len(),
            volume: self
                .tets
                .iter()
                .map(|t| determinant(t.map(|i| self.points[i])) / 6.)
                .sum(),
            mean_ratio_min: q[0],
            mean_ratio_p01: q[q.len() / 100],
            mean_ratio_median: q[q.len() / 2],
            mean_ratio_below_01: q.iter().filter(|&&v| v < 0.1).count(),
            mmg_quality_min: q[0].powf(1.5),
            mmg_quality_p01: q[q.len() / 100].powf(1.5),
            mmg_quality_median: q[q.len() / 2].powf(1.5),
            mmg_quality_below_01: q.iter().filter(|&&v| v.powf(1.5) < 0.1).count(),
            minimum_dihedral_degrees: self
                .tets
                .iter()
                .map(|t| min_dihedral(t.map(|i| self.points[i])))
                .fold(180., f64::min),
            maximum_surface_node_residual: residual,
            periodic_nodes_and_triangles: self.periodicity(),
        }
    }
}

mod adaptive;
mod optimize;

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn aligned_optimization_preserves_geometry_and_periodicity() {
        for (resolution, threshold) in [(24, 0.5), (25, 0.6)] {
            let options = GyroidOptions {
                resolution,
                threshold,
                align_cube_edges: true,
                ..GyroidOptions::default()
            };
            let before = generate_gyroid(options).unwrap().metrics(threshold);
            let after = generate_gyroid(GyroidOptions {
                optimize_passes: 2,
                ..options
            })
            .unwrap()
            .metrics(threshold);
            assert_eq!(before.tetrahedra, after.tetrahedra);
            assert!(after.mean_ratio_min >= before.mean_ratio_min - 1e-10);
            assert!(after.maximum_surface_node_residual < 1e-9);
            assert!(after.periodic_nodes_and_triangles.iter().all(|&v| v));
        }
    }
    #[test]
    fn rejects_unresolved_edge_alignment() {
        assert!(
            generate_gyroid(GyroidOptions {
                threshold: 1.,
                align_cube_edges: true,
                ..GyroidOptions::default()
            })
            .is_err()
        );
    }
    #[test]
    fn adaptive_parallel_is_deterministic() {
        let o = GyroidOptions {
            resolution: 24,
            align_cube_edges: true,
            optimize_passes: 2,
            post_optimize_passes: 2,
            surface_weight: 3.,
            surface_tolerance: 0.0008,
            threads: 1,
            ..GyroidOptions::default()
        };
        let serial = generate_gyroid(o).unwrap();
        let parallel = generate_gyroid(GyroidOptions { threads: 4, ..o }).unwrap();
        assert_eq!(serial.points, parallel.points);
        assert_eq!(serial.tets, parallel.tets);
        assert!(parallel.periodicity().iter().all(|&v| v));
        assert!(
            parallel
                .tets
                .iter()
                .all(|t| determinant(t.map(|i| parallel.points[i])) > 0.)
        );
        assert!(parallel.metrics(o.threshold).maximum_surface_node_residual < 1e-9);
        let unrefined = generate_gyroid(GyroidOptions {
            surface_tolerance: 0.,
            ..o
        })
        .unwrap();
        assert!(parallel.tets.len() > unrefined.tets.len());
        assert!(
            generate_gyroid(GyroidOptions {
                max_tetrahedra: 1,
                ..o
            })
            .is_err()
        );
    }
    #[test]
    fn periodic_band() {
        let mesh = generate_gyroid(GyroidOptions {
            resolution: 24,
            ..GyroidOptions::default()
        })
        .unwrap();
        let m = mesh.metrics(0.5);
        assert!(m.periodic_nodes_and_triangles.iter().all(|&v| v));
        assert!(m.mean_ratio_min > 0.);
        assert!(m.maximum_surface_node_residual < 1e-9);
        assert!((m.volume - 0.32356).abs() < 0.01);
        assert!(
            mesh.points
                .iter()
                .flatten()
                .all(|&x| (0.0..=1.0).contains(&x))
        );
    }
    #[test]
    fn rejects_underresolved_sheet() {
        assert!(
            generate_gyroid(GyroidOptions {
                resolution: 8,
                ..GyroidOptions::default()
            })
            .is_err()
        );
    }
    #[test]
    fn regular_tet_quality() {
        let p = [
            [0., 0., 0.],
            [1., 0., 0.],
            [0.5, 3_f64.sqrt() / 2., 0.],
            [0.5, 3_f64.sqrt() / 6., (2_f64 / 3.).sqrt()],
        ];
        assert!((quality(p) - 1.).abs() < 1e-12);
        assert!((min_dihedral(p) - 70.528779).abs() < 1e-5);
    }
    #[test]
    fn odd_resolution_and_other_band() {
        for o in [
            GyroidOptions {
                resolution: 25,
                ..GyroidOptions::default()
            },
            GyroidOptions {
                resolution: 32,
                threshold: 0.7,
                periodic: false,
                ..GyroidOptions::default()
            },
        ] {
            let mesh = generate_gyroid(o).unwrap();
            let m = mesh.metrics(o.threshold);
            assert!(m.maximum_surface_node_residual < 1e-9);
            assert!(m.mean_ratio_min > 0.0);
        }
    }
}

#[cfg(test)]
mod boundary_tests {
    use super::*;
    #[test]
    fn flat_boundary_rejects_same_orientation_and_nonmanifold_faces() {
        let t = [0, 1, 2, 3];
        assert!(boundary(&[t, t]).unwrap_err().contains("orientation"));
        assert!(boundary(&[t, t, t]).unwrap_err().contains("nonmanifold"));
    }
}

fn validate_options(o: GyroidOptions) -> Result<(), String> {
    if !o.surface_weight.is_finite()
        || o.surface_weight < 0.
        || !o.surface_tolerance.is_finite()
        || o.surface_tolerance < 0.
        || o.threads > 256
    {
        return Err(
            "surface parameters must be finite and nonnegative; threads must be at most 256".into(),
        );
    }

    if !(8..=192).contains(&o.resolution)
        || !o.threshold.is_finite()
        || !(0.1..=1.).contains(&o.threshold)
        || !o.snap.is_finite()
        || !(0.0..=0.2).contains(&o.snap)
    {
        return Err("resolution must be 8..192, threshold 0.1..1, snap 0..0.2".into());
    }
    Ok(())
}
