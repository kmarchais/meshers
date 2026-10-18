//! Quadratic (10-node) tetrahedra whose wall edges follow the exact implicit surface.
//!
//! Mid-edge nodes of implicit-wall edges are Newton-projected onto the wall, holding
//! any box coordinate the edge lies in, so the geometric error of the boundary drops
//! from O(h^2) to O(h^3). Interior and box-plane mid-edge nodes stay at midpoints.
//! Each element's Jacobian determinant is sampled at its nodes, face and cell centres
//! and fifteen interior points; corrections of failing elements are repeatedly halved,
//! together with their periodic copies. Unresolved edges are straightened after eight
//! reductions. Sampling is a strong check but not a proof of positivity everywhere.
use crate::{Mesh, Point, geometry::Geometry, norm};
use std::collections::{HashMap, HashSet};

/// VTK / fedoo / Gmsh-compatible mid-edge order after the four vertices.
pub const EDGES: [[usize; 2]; 6] = [[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]];

pub struct Quadratic {
    /// Vertices of the linear mesh followed by mid-edge nodes.
    pub points: Vec<Point>,
    pub tets: Vec<[usize; 10]>,
    /// Boundary triangles: three vertices, then mid-edge nodes of (0,1), (1,2), (2,0).
    pub surface: Vec<[usize; 6]>,
    pub curved_edges: usize,
    pub straightened_edges: usize,
    pub partially_curved_edges: usize,
    /// Smallest ratio of minimum to maximum sampled Jacobian over all elements.
    pub minimum_jacobian_ratio: f64,
}

fn key(a: usize, b: usize) -> [usize; 2] {
    [a.min(b), a.max(b)]
}

/// Sample points in barycentric coordinates: nodes, face centres, centroid and a
/// symmetric interior set.
fn samples() -> Vec<[f64; 4]> {
    let mut out = Vec::new();
    for i in 0..4 {
        let mut l = [0.; 4];
        l[i] = 1.;
        out.push(l);
    }
    for [i, j] in EDGES {
        let mut l = [0.; 4];
        l[i] = 0.5;
        l[j] = 0.5;
        out.push(l);
    }
    for skip in 0..4 {
        let mut l = [1. / 3.; 4];
        l[skip] = 0.;
        out.push(l);
    }
    out.push([0.25; 4]);
    for (a, b) in [(0.7, 0.1), (0.1, 0.3)] {
        for i in 0..4 {
            let mut l = [b; 4];
            l[i] = a;
            out.push(l);
        }
    }
    out
}

/// Jacobian determinant of a quadratic tetrahedron at barycentric point `l`.
fn jacobian(x: &[Point; 10], l: [f64; 4]) -> f64 {
    // dN/dL for the ten shape functions; reference coordinates r,s,t = L1,L2,L3.
    let mut grad = [[0.; 3]; 10];
    let dl = |i: usize| -> [f64; 3] {
        match i {
            0 => [-1., -1., -1.],
            1 => [1., 0., 0.],
            2 => [0., 1., 0.],
            _ => [0., 0., 1.],
        }
    };
    for i in 0..4 {
        let d = dl(i);
        let f = 4. * l[i] - 1.;
        grad[i] = d.map(|v| f * v);
    }
    for (k, [i, j]) in EDGES.iter().enumerate() {
        let (di, dj) = (dl(*i), dl(*j));
        grad[4 + k] = std::array::from_fn(|c| 4. * (l[*j] * di[c] + l[*i] * dj[c]));
    }
    let mut jac = [[0.; 3]; 3];
    for n in 0..10 {
        for a in 0..3 {
            for c in 0..3 {
                jac[a][c] += x[n][a] * grad[n][c];
            }
        }
    }
    jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1])
        - jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0])
        + jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0])
}

fn project(
    g: &impl Geometry,
    p: Point,
    level: f64,
    locked: [bool; 3],
    reach: f64,
) -> Option<Point> {
    let p = crate::projection::onto_level(g, p, level, locked, reach)?;
    g.contains(p).then_some(p)
}

/// Elevate a validated linear mesh. `tags` marks implicit-wall boundary faces with 0.
pub(crate) fn elevate(
    mesh: &Mesh,
    tags: &[u8],
    g: &impl Geometry,
    threshold: f64,
    periodic: [bool; 3],
) -> Result<Quadratic, String> {
    let [lo, hi] = g.bounds();
    let mut index = HashMap::<[usize; 2], usize>::new();
    let mut ends = Vec::new();
    let mut tets = Vec::with_capacity(mesh.tets.len());
    for t in &mesh.tets {
        let mut q = [0; 10];
        q[..4].copy_from_slice(t);
        for (k, [i, j]) in EDGES.iter().enumerate() {
            let e = key(t[*i], t[*j]);
            q[4 + k] = mesh.points.len()
                + *index.entry(e).or_insert_with(|| {
                    ends.push(e);
                    ends.len() - 1
                });
        }
        tets.push(q);
    }
    let straight: Vec<Point> = ends
        .iter()
        .map(|&[a, b]| std::array::from_fn(|c| 0.5 * (mesh.points[a][c] + mesh.points[b][c])))
        .collect();
    let mut wall = HashSet::new();
    for (f, &tag) in mesh.surface.iter().zip(tags) {
        if tag == 0 {
            for k in 0..3 {
                wall.insert(key(f[k], f[(k + 1) % 3]));
            }
        }
    }
    let mut mid = straight.clone();
    let mut curved = vec![false; ends.len()];
    g.check()?;
    for (e, &[a, b]) in ends.iter().enumerate() {
        if !wall.contains(&[a, b]) {
            continue;
        }
        let (pa, pb) = (mesh.points[a], mesh.points[b]);
        let locked: [bool; 3] =
            std::array::from_fn(|c| pa[c] == pb[c] && (pa[c] == lo[c] || pa[c] == hi[c]));
        let level = g.level(g.value(pa), threshold);
        let reach = 0.5 * norm(crate::sub(pa, pb));
        if let Some(p) = project(g, straight[e], level, locked, reach) {
            mid[e] = p;
            curved[e] = true;
        }
    }
    // Periodic copies of a box-plane edge share one canonical straight midpoint.
    let canonical = |p: Point| -> [i64; 3] {
        std::array::from_fn(|c| {
            let x = if periodic[c] && p[c] == hi[c] {
                lo[c]
            } else {
                p[c]
            };
            ((x - lo[c]) / (hi[c] - lo[c]) * 1e11).round() as i64
        })
    };
    let points_of = |q: &[usize; 10], mid: &[Point]| -> [Point; 10] {
        std::array::from_fn(|n| {
            if n < 4 {
                mesh.points[q[n]]
            } else {
                mid[q[n] - mesh.points.len()]
            }
        })
    };
    let sample = samples();
    let mut straightened = 0;
    let mut reductions = vec![0u8; ends.len()];
    let mut minimum_ratio;
    loop {
        g.check()?;
        minimum_ratio = f64::INFINITY;
        let mut bad = HashSet::new();
        for q in &tets {
            if (4..10).all(|n| !curved[q[n] - mesh.points.len()]) {
                continue;
            }
            let x = points_of(q, &mid);
            let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
            for &l in &sample {
                let d = jacobian(&x, l);
                low = low.min(d);
                high = high.max(d);
            }
            let ratio = if high > 0. {
                low / high
            } else {
                f64::NEG_INFINITY
            };
            minimum_ratio = minimum_ratio.min(ratio);
            if ratio.is_nan() || ratio <= 1e-3 {
                for &node in &q[4..] {
                    let e = node - mesh.points.len();
                    if curved[e] {
                        bad.insert(canonical(straight[e]));
                    }
                }
            }
        }
        if bad.is_empty() {
            break;
        }
        for e in 0..ends.len() {
            if curved[e] && bad.contains(&canonical(straight[e])) {
                if reductions[e] < 8 {
                    mid[e] = std::array::from_fn(|axis| 0.5 * (mid[e][axis] + straight[e][axis]));
                    reductions[e] += 1;
                } else {
                    curved[e] = false;
                    mid[e] = straight[e];
                    straightened += 1;
                }
            }
        }
    }
    if minimum_ratio == f64::INFINITY {
        minimum_ratio = 1.;
    }
    let mut surface = Vec::with_capacity(mesh.surface.len());
    for f in &mesh.surface {
        let mut s = [0; 6];
        s[..3].copy_from_slice(f);
        for k in 0..3 {
            s[3 + k] = mesh.points.len() + index[&key(f[k], f[(k + 1) % 3])];
        }
        surface.push(s);
    }
    let mut points = mesh.points.clone();
    points.extend(mid);
    if points.iter().flatten().any(|v| !v.is_finite()) {
        return Err("nonfinite quadratic node".into());
    }
    Ok(Quadratic {
        points,
        tets,
        surface,
        curved_edges: curved.iter().filter(|&&c| c).count(),
        straightened_edges: straightened,
        partially_curved_edges: (0..ends.len())
            .filter(|&edge| curved[edge] && reductions[edge] > 0)
            .count(),
        minimum_jacobian_ratio: minimum_ratio,
    })
}

/// Opposite node pairs over all quadratic nodes, per periodic axis.
pub(crate) fn node_pairs(
    points: &[Point],
    bounds: [Point; 2],
    periodic: [bool; 3],
) -> Result<[Vec<[usize; 2]>; 3], String> {
    let [lo, hi] = bounds;
    let mut pairs: [Vec<[usize; 2]>; 3] = std::array::from_fn(|_| Vec::new());
    for (a, pair) in pairs.iter_mut().enumerate() {
        if !periodic[a] {
            continue;
        }
        let key = |p: Point| -> [i64; 3] {
            std::array::from_fn(|c| {
                if c == a {
                    0
                } else {
                    ((p[c] - lo[c]) / (hi[c] - lo[c]) * 1e11).round() as i64
                }
            })
        };
        let mut maps = [
            std::collections::BTreeMap::new(),
            std::collections::BTreeMap::new(),
        ];
        for (i, &p) in points.iter().enumerate() {
            for (side, value) in [lo[a], hi[a]].iter().enumerate() {
                if p[a] == *value && maps[side].insert(key(p), i).is_some() {
                    return Err("ambiguous periodic quadratic node match".into());
                }
            }
        }
        if maps[0].keys().ne(maps[1].keys()) {
            return Err(format!("periodic quadratic node mismatch on axis {a}"));
        }
        for (k, &i) in &maps[0] {
            pair.push([i, maps[1][k]]);
        }
    }
    Ok(pairs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dot;

    struct ConcaveSphere {
        bounds: [Point; 2],
    }

    impl Geometry for ConcaveSphere {
        fn value(&self, point: Point) -> f64 {
            1. - norm(point)
        }
        fn gradient(&self, point: Point) -> Point {
            point.map(|coordinate| -coordinate / norm(point))
        }
        fn bounds(&self) -> [Point; 2] {
            self.bounds
        }
        fn rank(&self, _point: Point) -> [i64; 3] {
            [0; 3]
        }
        fn band(&self) -> bool {
            false
        }
    }

    #[test]
    fn validity_fallback_retains_safe_partial_curvature() {
        let height = 0.91f64.sqrt();
        let tetrahedra = vec![[0, 1, 2, 3]];
        let surface = crate::boundary(&tetrahedra).unwrap();
        let tags: Vec<_> = surface
            .iter()
            .map(|face| u8::from(face.contains(&3)))
            .collect();
        let mesh = Mesh {
            points: vec![
                [-0.3, 0., height],
                [0.3, 0., height],
                [0., 0.3, height],
                [0., 0.1, 1.01],
            ],
            tets: tetrahedra,
            surface,
        };
        let geometry = ConcaveSphere {
            bounds: [[-2.; 3], [2.; 3]],
        };
        let quadratic = elevate(&mesh, &tags, &geometry, 0., [false; 3]).unwrap();
        assert!(
            quadratic.curved_edges > 0,
            "all {} wall edges were straightened",
            quadratic.straightened_edges
        );
        assert!(quadratic.minimum_jacobian_ratio > 1e-3);
        assert!(quadratic.partially_curved_edges > 0);
        let straight: Point =
            std::array::from_fn(|axis| (mesh.points[0][axis] + mesh.points[1][axis]) / 2.);
        let midpoint = quadratic.points[quadratic.tets[0][4]];
        assert!((norm(midpoint) - 1.).abs() < (norm(straight) - 1.).abs());
        let coordinates = quadratic.tets[0].map(|node| quadratic.points[node]);
        assert!(
            samples()
                .into_iter()
                .all(|weights| jacobian(&coordinates, weights) > 0.)
        );
    }

    #[test]
    fn partial_curving_keeps_periodic_copies_identical() {
        let height = 0.66f64.sqrt();
        let mut points = vec![
            [-0.5, -0.3, height],
            [-0.5, 0.3, height],
            [-0.4, 0., 0.84f64.sqrt()],
            [-0.4, 0., 0.95],
        ];
        let upper: Vec<_> = points
            .iter()
            .map(|point| [-point[0], point[1], point[2]])
            .collect();
        points.extend(upper);
        points[7][2] = 1.1;
        let mut tetrahedra = vec![[0, 1, 2, 3], [4, 5, 6, 7]];
        for tetrahedron in &mut tetrahedra {
            if crate::determinant((*tetrahedron).map(|node| points[node])) < 0. {
                tetrahedron.swap(0, 1);
            }
        }
        let surface = crate::boundary(&tetrahedra).unwrap();
        let tags: Vec<_> = surface
            .iter()
            .map(|face| u8::from(face.contains(&3) || face.contains(&7)))
            .collect();
        let mesh = Mesh {
            points,
            tets: tetrahedra,
            surface,
        };
        let geometry = ConcaveSphere {
            bounds: [[-0.5, -2., -2.], [0.5, 2., 2.]],
        };
        let quadratic = elevate(&mesh, &tags, &geometry, 0., [true, false, false]).unwrap();
        assert!(quadratic.partially_curved_edges > 0);
        let pairs = node_pairs(&quadratic.points, geometry.bounds(), [true, false, false]).unwrap();
        assert_eq!(pairs[0].len(), 3);
        for &[lower, upper] in &pairs[0] {
            assert_eq!(
                crate::sub(quadratic.points[upper], quadratic.points[lower]),
                [1., 0., 0.]
            );
        }
        let midpoint = pairs[0]
            .iter()
            .map(|pair| quadratic.points[pair[0]])
            .find(|point| point[1] == 0.)
            .unwrap();
        assert!(midpoint[2] > height + 1e-6);
    }

    #[test]
    fn wall_mid_edge_nodes_lie_on_the_sphere() {
        use crate::implicit::{Options, generate};
        let field = |p: Point| dot(p, p).sqrt() - 0.35;
        let out = generate(
            &field,
            Options {
                bounds: [[-0.5; 3], [0.5; 3]],
                cells: [12; 3],
                geometry_tolerance: 0.1,
                optimize_passes: 1,
                order: 2,
                ..Default::default()
            },
        )
        .unwrap();
        let q = out.quadratic.unwrap().elements;
        assert!(q.curved_edges > 0 && q.minimum_jacobian_ratio > 0.);
        let mut worst_straight = 0f64;
        for f in &q.surface {
            for k in 0..3 {
                let m = q.points[f[3 + k]];
                assert!(field(m).abs() < 1e-12, "{}", field(m));
                let (a, b) = (q.points[f[k]], q.points[f[(k + 1) % 3]]);
                let straight = std::array::from_fn(|c| 0.5 * (a[c] + b[c]));
                worst_straight = worst_straight.max(field(straight).abs());
            }
        }
        assert!(
            worst_straight > 1e-4,
            "linear chords should sag: {worst_straight}"
        );
    }

    #[test]
    fn periodic_quadratic_nodes_pair_on_every_axis() {
        use crate::implicit::{Options, Region, generate};
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
                cells: [10; 3],
                region: Region::Band {
                    lower: -0.55,
                    upper: 0.55,
                },
                periodic: [true; 3],
                geometry_tolerance: 0.1,
                optimize_passes: 1,
                order: 2,
                ..Default::default()
            },
        )
        .unwrap();
        let q = out.quadratic.unwrap();
        for (axis, pairs) in q.periodic_pairs.iter().enumerate() {
            // More pairs than the linear mesh: mid-edge nodes on the faces pair too.
            assert!(pairs.len() > out.periodic_pairs[axis].len());
            for &[a, b] in pairs {
                let (pa, pb) = (q.elements.points[a], q.elements.points[b]);
                assert_eq!(pb[axis] - pa[axis], 1.);
                for c in (0..3).filter(|&c| c != axis) {
                    assert!((pa[c] - pb[c]).abs() < 1e-12);
                }
            }
        }
    }

    #[test]
    fn straight_jacobian_is_six_times_volume_everywhere() {
        let v = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
        let x: [Point; 10] = std::array::from_fn(|n| {
            if n < 4 {
                v[n]
            } else {
                let [i, j] = EDGES[n - 4];
                std::array::from_fn(|c| 0.5 * (v[i][c] + v[j][c]))
            }
        });
        for l in samples() {
            assert!((jacobian(&x, l) - 1.).abs() < 1e-12);
        }
    }
}
