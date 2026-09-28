//! Opt-in research variants. No effect without optimizer-experiments feature.
use super::*;
use crate::cross;
use std::sync::OnceLock;

fn mode() -> &'static str {
    static MODE: OnceLock<String> = OnceLock::new();
    MODE.get_or_init(|| std::env::var("MESHER_OPT_EXPERIMENT").unwrap_or_default())
}

/// Derivative of sum(1/q²), where q is squared MMG quality. Periodic
/// copies in one tetrahedron contribute to the same translation derivative.
fn penalty_gradient(mesh: &Mesh, group: &Group) -> Option<Point> {
    let mut gradient = [0.; 3];
    for &ci in &group.1 {
        let ids = mesh.tets[ci];
        let p = ids.map(|i| mesh.points[i]);
        let d = determinant(p);
        if d <= 0. {
            return None;
        }
        let u = sub(p[1], p[0]);
        let v = sub(p[2], p[0]);
        let w = sub(p[3], p[0]);
        let d1 = cross(v, w);
        let d2 = cross(w, u);
        let d3 = cross(u, v);
        let ds = [std::array::from_fn(|a| -d1[a] - d2[a] - d3[a]), d1, d2, d3];
        let mut length = 0.;
        let mut dl = [0.; 3];
        let mut dd = [0.; 3];
        for i in 0..4 {
            if group.0.contains(&ids[i]) {
                for a in 0..3 {
                    dd[a] += ds[i][a];
                }
            }
            for j in i + 1..4 {
                let edge = sub(p[i], p[j]);
                length += dot(edge, edge);
                let moving =
                    i32::from(group.0.contains(&ids[i])) - i32::from(group.0.contains(&ids[j]));
                for a in 0..3 {
                    dl[a] += 2. * f64::from(moving) * edge[a];
                }
            }
        }
        let q = 432. * d * d / length.powi(3);
        let cost = 1. / (q * q);
        for a in 0..3 {
            gradient[a] += cost * (6. * dl[a] / length - 4. * dd[a] / d);
        }
    }
    gradient.iter().all(|v| v.is_finite()).then_some(gradient)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn direction(
    mesh: &Mesh,
    group: &Group,
    threshold: f64,
    h: f64,
    weight: f64,
    level: Option<f64>,
    geometry: &impl crate::geometry::Geometry,
) -> Option<Point> {
    if !mode().contains("analytic") {
        return None;
    }
    let origin = mesh.points[group.0[0]];
    let locked = geometry.locked(origin);
    let mut g = penalty_gradient(mesh, group)?;
    for a in 0..3 {
        if locked[a] {
            g[a] = 0.;
        }
    }
    if level.is_some() {
        let mut normal = geometry.gradient(origin);
        for a in 0..3 {
            if locked[a] {
                normal[a] = 0.;
            }
        }
        let n2 = dot(normal, normal);
        if n2 > 1e-16 {
            let component = dot(g, normal) / n2;
            for a in 0..3 {
                g[a] -= normal[a] * component;
            }
        }
    }
    let eps = h * 0.002;
    let mut direction = g.map(|v| -2. * eps * v);
    // Surface error depends on field curvature. Retain the existing projected
    // finite differences for that term, while avoiding repeated tetrahedron work.
    if !group.2.is_empty() && weight > 0. {
        let surface = (group.0.clone(), Vec::new(), group.2.clone());
        let cost = |p| evaluate(mesh, &surface, threshold, h, weight, p, geometry).1;
        let old = cost(origin);
        for a in 0..3 {
            let mut lo = origin;
            lo[a] -= eps;
            let mut hi = origin;
            hi[a] += eps;
            let c0 = project(lo, origin, level, geometry).map_or(old, cost);
            let c1 = project(hi, origin, level, geometry).map_or(old, cost);
            if !c0.is_finite() || !c1.is_finite() {
                return None;
            }
            direction[a] += c0 - c1;
        }
    }
    direction.iter().all(|v| v.is_finite()).then_some(direction)
}

pub(super) struct Active {
    enabled: bool,
    dirty: Vec<bool>,
    adjacent: Vec<Vec<usize>>,
    skipped: usize,
    visited: usize,
}
impl Active {
    pub(super) fn new(mesh: &Mesh, groups: &[Group]) -> Self {
        Self::with_enabled(mesh, groups, mode().contains("active"))
    }
    fn with_enabled(mesh: &Mesh, groups: &[Group], enabled: bool) -> Self {
        let mut adjacent = vec![Vec::new(); if enabled { groups.len() } else { 0 }];
        if enabled {
            let mut owner = vec![0; mesh.points.len()];
            for (gi, group) in groups.iter().enumerate() {
                for &i in &group.0 {
                    owner[i] = gi;
                }
            }
            for (gi, group) in groups.iter().enumerate() {
                for &ci in &group.1 {
                    for &i in &mesh.tets[ci] {
                        adjacent[gi].push(owner[i]);
                    }
                }
                adjacent[gi].sort_unstable();
                adjacent[gi].dedup();
            }
        }
        Self {
            enabled,
            dirty: vec![true; groups.len()],
            adjacent,
            skipped: 0,
            visited: 0,
        }
    }
    pub(super) fn begin_pass(&mut self, pass: usize) {
        // Proposal step length changes at pass 3, so previous rejections expire.
        if self.enabled && pass == 3 {
            self.dirty.fill(true);
        }
    }
    pub(super) fn visit(&mut self, i: usize) -> bool {
        if self.enabled && !self.dirty[i] {
            self.skipped += 1;
            false
        } else {
            self.visited += 1;
            true
        }
    }
    pub(super) fn update(&mut self, i: usize, delta: Point) {
        if !self.enabled {
            return;
        }
        if delta == [0.; 3] {
            self.dirty[i] = false;
        } else {
            for &j in &self.adjacent[i] {
                self.dirty[j] = true;
            }
        }
    }
}
impl Drop for Active {
    fn drop(&mut self) {
        if self.enabled && std::env::var_os("MESHER_CPU_PROFILE").is_some() {
            eprintln!("ACTIVE visited={} skipped={}", self.visited, self.skipped);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn active_set_invalidates_periodic_neighbours_and_step_changes() {
        let mesh = Mesh {
            points: vec![[0.; 3]; 6],
            tets: vec![[0, 1, 2, 3], [4, 1, 2, 3]],
            surface: vec![],
        };
        let groups = vec![
            (vec![0, 4], vec![0, 1], vec![]),
            (vec![1], vec![0, 1], vec![]),
            (vec![2], vec![0, 1], vec![]),
            (vec![3], vec![0, 1], vec![]),
            (vec![5], vec![], vec![]),
        ];
        let mut active = Active::with_enabled(&mesh, &groups, true);
        active.update(0, [0.; 3]);
        active.update(4, [0.; 3]);
        assert!(!active.visit(0));
        assert!(!active.visit(4));
        active.update(1, [0.01, 0., 0.]);
        assert!(active.visit(0));
        assert!(!active.visit(4));
        active.begin_pass(3);
        assert!(active.visit(4));
    }
    #[test]
    fn penalty_derivative_matches_numerical_orbit_translation() {
        let mesh = Mesh {
            points: vec![
                [0.1, 0.2, 0.1],
                [0.8, 0.1, 0.2],
                [0.2, 0.9, 0.1],
                [0.2, 0.1, 0.9],
            ],
            tets: vec![[0, 1, 2, 3]],
            surface: vec![],
        };
        for members in [vec![0], vec![1], vec![2, 3], vec![0, 1, 2, 3]] {
            let group = (members, vec![0], vec![]);
            let g = penalty_gradient(&mesh, &group).unwrap();
            let origin = mesh.points[group.0[0]];
            for a in 0..3 {
                let mut lo = origin;
                let mut hi = origin;
                lo[a] -= 1e-6;
                hi[a] += 1e-6;
                let cost = |p| evaluate(&mesh, &group, 0.5, 1., 0., p, &crate::geometry::Legacy).1;
                let fd = (cost(hi) - cost(lo)) / 2e-6;
                assert!(
                    (g[a] - fd).abs() < 1e-5 * (1. + fd.abs()),
                    "analytic={} fd={}",
                    g[a],
                    fd
                );
            }
        }
    }
}
