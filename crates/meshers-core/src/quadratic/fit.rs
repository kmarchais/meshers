//! Bounded fitting of complete quadratic wall patches, with fixed corners.
use super::{Quadratic, jacobian, node_pairs, samples};
use crate::{Point, dot, geometry::Geometry, norm, sub};

#[derive(Clone, Debug, Default)]
pub struct FitStats {
    pub selected_groups: usize,
    pub affected_faces: usize,
    pub accepted_steps: usize,
    pub moved_nodes: usize,
    pub initial_distance_estimate: f64,
    pub final_distance_estimate: f64,
}

struct Face {
    nodes: [usize; 6],
    level: f64,
    area: f64,
}

fn root(parent: &mut [usize], mut n: usize) -> usize {
    while parent[n] != n {
        parent[n] = parent[parent[n]];
        n = parent[n];
    }
    n
}

fn basis() -> Vec<[f64; 6]> {
    let mut out = Vec::new();
    for i in 0..=6 {
        for j in 0..=6 - i {
            let [a, b, c] = [i as f64 / 6., j as f64 / 6., (6 - i - j) as f64 / 6.];
            out.push([
                a * (2. * a - 1.),
                b * (2. * b - 1.),
                c * (2. * c - 1.),
                4. * a * b,
                4. * b * c,
                4. * c * a,
            ]);
        }
    }
    out
}

fn observations(
    faces: &[&Face],
    points: &[Point],
    basis: &[[f64; 6]],
    g: &impl Geometry,
) -> Result<(Vec<f64>, Vec<Point>), String> {
    let bounds = g.bounds();
    let scale = norm(sub(bounds[1], bounds[0]));
    let mut locations = Vec::with_capacity(faces.len() * basis.len());
    let mut locks = Vec::with_capacity(locations.capacity());
    for face in faces {
        for weights in basis {
            let mut p: Point = std::array::from_fn(|a| {
                (0..6).map(|i| weights[i] * points[face.nodes[i]][a]).sum()
            });
            let mut locked = [false; 3];
            for a in 0..3 {
                for bound in bounds {
                    if (p[a] - bound[a]).abs() <= 1e-12 * scale {
                        p[a] = bound[a];
                        locked[a] = true;
                    }
                }
            }
            locations.push(p);
            locks.push(locked);
        }
    }
    g.check()?;
    let values = g.values(&locations);
    let mut gradients = g.gradients(&locations);
    let mut residuals = Vec::with_capacity(locations.len());
    for i in 0..locations.len() {
        let residual = values[i] - faces[i / basis.len()].level;
        let allowed: Point =
            std::array::from_fn(|a| if locks[i][a] { 0. } else { gradients[i][a] });
        let length = norm(allowed);
        if length > 1e-14 {
            residuals.push(residual / length);
            gradients[i] = gradients[i].map(|v| v / length);
        } else if residual.abs() <= 1e-12 {
            residuals.push(0.);
            gradients[i] = [0.; 3];
        } else {
            return Err("curved fitting encountered a stationary off-wall sample".into());
        }
    }
    Ok((residuals, gradients))
}

fn maximum(values: &[f64]) -> f64 {
    values.iter().map(|v| v.abs()).fold(0., f64::max)
}

fn validity_samples() -> Vec<[f64; 4]> {
    let mut out = samples();
    // Positive five-point Gauss-Legendre grid under the tetrahedral Duffy map.
    let gauss = [
        0.046910077030668,
        0.230765344947158,
        0.5,
        0.769234655052842,
        0.953089922969332,
    ];
    for u in gauss {
        for v in gauss {
            for w in gauss {
                let r = u;
                let s = (1. - u) * v;
                let t = (1. - u) * (1. - v) * w;
                out.push([1. - r - s - t, r, s, t]);
            }
        }
    }
    out
}

pub(crate) fn improve(
    mesh: &mut Quadratic,
    g: &impl Geometry,
    threshold: f64,
    periodic: [bool; 3],
) -> Result<FitStats, String> {
    let original = mesh.points.clone();
    let n = original.len();
    let mut parent: Vec<_> = (0..n).collect();
    for pairs in node_pairs(&original, &mesh.surface, g.bounds(), periodic)? {
        for [a, b] in pairs {
            let a = root(&mut parent, a);
            let b = root(&mut parent, b);
            parent[a.max(b)] = a.min(b);
        }
    }
    let classes: Vec<_> = (0..n).map(|i| root(&mut parent, i)).collect();
    let faces: Vec<_> = mesh
        .surface
        .iter()
        .filter(|f| !g.cap([original[f[0]], original[f[1]], original[f[2]]]))
        .map(|&nodes| {
            let center =
                std::array::from_fn(|a| (0..3).map(|i| original[nodes[i]][a]).sum::<f64>() / 3.);
            let a = sub(original[nodes[1]], original[nodes[0]]);
            let b = sub(original[nodes[2]], original[nodes[0]]);
            let cross = [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ];
            Face {
                nodes,
                level: g.level(g.value(center), threshold),
                area: norm(cross) / 2.,
            }
        })
        .collect();
    if faces.is_empty() {
        return Ok(FitStats::default());
    }
    let mut mids: Vec<_> = faces
        .iter()
        .flat_map(|f| f.nodes[3..].iter().copied())
        .collect();
    mids.sort_unstable();
    mids.dedup();
    let mut direction = vec![[0.; 3]; n];
    let locations: Vec<_> = mids.iter().map(|&m| original[m]).collect();
    for (&m, gradient) in mids.iter().zip(g.gradients(&locations)) {
        let locked = g.locked(original[m]);
        for a in 0..3 {
            if !locked[a] {
                direction[classes[m]][a] += gradient[a];
            }
        }
    }
    for d in &mut direction {
        let length = norm(*d);
        if length > 0. {
            *d = d.map(|v| v / length);
        }
    }
    let mut limits = vec![f64::INFINITY; n];
    for face in &faces {
        for (k, [a, b]) in [[0, 1], [1, 2], [2, 0]].into_iter().enumerate() {
            let group = classes[face.nodes[3 + k]];
            limits[group] = limits[group]
                .min(0.02 * norm(sub(original[face.nodes[a]], original[face.nodes[b]])));
        }
    }
    let basis = basis();
    let all: Vec<_> = faces.iter().collect();
    let (initial, _) = observations(&all, &original, &basis, g)?;
    let initial_max = maximum(&initial);
    let mut stats = FitStats {
        initial_distance_estimate: initial_max,
        final_distance_estimate: initial_max,
        ..FitStats::default()
    };
    if initial_max <= 1e-14 * norm(sub(g.bounds()[1], g.bounds()[0])) {
        return Ok(stats);
    }
    let mut active = vec![false; n];
    for (face, values) in faces.iter().zip(initial.chunks_exact(basis.len())) {
        if maximum(values) >= 0.5 * initial_max {
            for &m in &face.nodes[3..] {
                active[classes[m]] = true;
            }
        }
    }
    stats.selected_groups = active.iter().filter(|&&v| v).count();
    let support: Vec<_> = faces
        .iter()
        .filter(|f| f.nodes[3..].iter().any(|&m| active[classes[m]]))
        .collect();
    stats.affected_faces = support.len();
    let affected: Vec<_> = mesh
        .tets
        .iter()
        .filter(|q| q[4..].iter().any(|&m| active[classes[m]]))
        .collect();
    let probes = validity_samples();
    let objective = |values: &[f64]| -> f64 {
        values
            .iter()
            .enumerate()
            .map(|(i, &r)| {
                support[i / basis.len()].area * (r / initial_max).powi(8) / basis.len() as f64
            })
            .sum()
    };
    let mut offsets = vec![0.; n];
    for _ in 0..6 {
        g.check()?;
        let (residuals, normals) = observations(&support, &mesh.points, &basis, g)?;
        let before = objective(&residuals);
        let rmax = maximum(&residuals);
        let mut rhs = vec![0.; n];
        let mut diagonal = vec![0.; n];
        for (i, &r) in residuals.iter().enumerate() {
            let face = support[i / basis.len()];
            let weight =
                face.area * (r.abs() / rmax.max(1e-30)).max(1e-3).powi(6) / basis.len() as f64;
            for k in 0..3 {
                let group = classes[face.nodes[3 + k]];
                if active[group] {
                    let coefficient =
                        basis[i % basis.len()][3 + k] * dot(normals[i], direction[group]);
                    rhs[group] -= weight * coefficient * r;
                    diagonal[group] += weight * coefficient * coefficient;
                }
            }
        }
        let floor = diagonal.iter().copied().fold(0., f64::max) * 1e-10;
        let step: Vec<_> = (0..n)
            .map(|i| {
                if active[i] && diagonal[i] > 0. {
                    (offsets[i] + rhs[i] / diagonal[i].max(floor)).clamp(-limits[i], limits[i])
                        - offsets[i]
                } else {
                    0.
                }
            })
            .collect();
        let mut accepted = false;
        for fraction in [1., 0.5, 0.25, 0.125] {
            let mut trial = original.clone();
            for &m in &mids {
                let group = classes[m];
                if !active[group] {
                    continue;
                }
                let locked = g.locked(original[m]);
                for a in 0..3 {
                    if !locked[a] {
                        // Copy common tangential coordinates, avoiding drift between aliases.
                        trial[m][a] = original[group][a]
                            + direction[group][a] * (offsets[group] + fraction * step[group]);
                    }
                }
            }
            let (values, _) = observations(&support, &trial, &basis, g)?;
            if objective(&values) >= before * (1. - 1e-10)
                || maximum(&values) > rmax + 1e-12 * norm(sub(g.bounds()[1], g.bounds()[0]))
            {
                continue;
            }
            let valid = affected.iter().all(|q| {
                let coordinates = q.map(|m| trial[m]);
                let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
                for &p in &probes {
                    let d = jacobian(&coordinates, p);
                    if !d.is_finite() {
                        return false;
                    }
                    low = low.min(d);
                    high = high.max(d);
                }
                low.is_finite() && high.is_finite() && low > 1e-3 * high && high > 0.
            });
            if !valid {
                continue;
            }
            for i in 0..n {
                offsets[i] += fraction * step[i];
            }
            mesh.points = trial;
            stats.accepted_steps += 1;
            accepted = true;
            break;
        }
        if !accepted {
            break;
        }
    }
    stats.final_distance_estimate = maximum(&observations(&all, &mesh.points, &basis, g)?.0);
    stats.moved_nodes = mids
        .iter()
        .filter(|&&m| mesh.points[m] != original[m])
        .count();
    if stats.moved_nodes > 0 {
        let probes = samples();
        mesh.minimum_jacobian_ratio = mesh
            .tets
            .iter()
            .map(|q| {
                let x = q.map(|m| mesh.points[m]);
                let values: Vec<_> = probes.iter().map(|&l| jacobian(&x, l)).collect();
                values.iter().copied().fold(f64::INFINITY, f64::min)
                    / values.iter().copied().fold(f64::NEG_INFINITY, f64::max)
            })
            .fold(1., f64::min);
    }
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Target {
        plane: bool,
    }
    impl Geometry for Target {
        fn value(&self, p: Point) -> f64 {
            if self.plane {
                p[2] - 0.2
            } else {
                dot(p, p) - 0.09
            }
        }
        fn gradient(&self, p: Point) -> Point {
            if self.plane {
                [0., 0., 1.]
            } else {
                p.map(|x| 2. * x)
            }
        }
        fn rank(&self, _: Point) -> [i64; 3] {
            [0; 3]
        }
        fn bounds(&self) -> [Point; 2] {
            [[-0.5; 3], [0.5; 3]]
        }
        fn band(&self) -> bool {
            false
        }
    }

    #[test]
    fn curved_fit_controls_exact_sphere_plane_and_face_order() {
        let mut points = vec![[0.2, 0., 1.], [0., 0.2, 1.], [-0.2, -0.2, 1.]];
        for p in &mut points {
            let length = norm(*p);
            *p = p.map(|v| v * 0.3 / length);
        }
        for [a, b] in [[0, 1], [1, 2], [2, 0]] {
            let mid: Point = std::array::from_fn(|i| points[a][i] + points[b][i]);
            points.push(mid.map(|v| v * 0.3 / norm(mid)));
        }
        let make = |points| Quadratic {
            points,
            tets: vec![],
            surface: vec![[0, 1, 2, 3, 4, 5]],
            curved_edges: 3,
            straightened_edges: 0,
            partially_curved_edges: 0,
            minimum_jacobian_ratio: 1.,
            fit: None,
        };
        let mut sphere = make(points.clone());
        let mut permuted = make(points.clone());
        permuted.surface[0] = [1, 2, 0, 4, 5, 3];
        let before = points.clone();
        let a = improve(&mut sphere, &Target { plane: false }, 0., [false; 3]).unwrap();
        improve(&mut permuted, &Target { plane: false }, 0., [false; 3]).unwrap();
        assert!(a.moved_nodes > 0);
        assert_eq!(sphere.points[..3], before[..3]);
        for (a, b) in sphere.points.iter().zip(&permuted.points) {
            assert!(norm(sub(*a, *b)) < 1e-12);
        }
        let error = |p: &[Point]| {
            let mut maximum = 0_f64;
            for i in 0..=41 {
                for j in 0..=41 - i {
                    let [a, b, c] = [i as f64 / 41., j as f64 / 41., (41 - i - j) as f64 / 41.];
                    let weights = [
                        a * (2. * a - 1.),
                        b * (2. * b - 1.),
                        c * (2. * c - 1.),
                        4. * a * b,
                        4. * b * c,
                        4. * c * a,
                    ];
                    let sample =
                        std::array::from_fn(|axis| (0..6).map(|n| weights[n] * p[n][axis]).sum());
                    maximum = maximum.max((norm(sample) - 0.3).abs());
                }
            }
            maximum
        };
        assert!(error(&sphere.points) < error(&before));
        for p in &mut points {
            p[2] = 0.2;
        }
        let mut plane = make(points.clone());
        assert_eq!(
            improve(&mut plane, &Target { plane: true }, 0., [false; 3])
                .unwrap()
                .moved_nodes,
            0
        );
        assert_eq!(points, plane.points);
    }
}
