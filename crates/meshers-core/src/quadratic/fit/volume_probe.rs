//! Private saved-mesh payoff experiment. No production meshing path calls this.
use super::{
    Face, Geometry, Point, Quadratic, basis, jacobian, maximum, node_pairs, root, validity_samples,
};
use crate::{cross, dot, norm, sub};
use serde::{Deserialize, Serialize};
use std::{
    f64::consts::TAU,
    fs::File,
    io::{BufReader, BufWriter},
    process::Command,
    time::Instant,
};

struct Fixture {
    name: String,
    periodic: [bool; 3],
    coefficients: Vec<[f64; 5]>,
    threshold: f64,
    knots: Vec<f64>,
    widths: Vec<f64>,
}

const FLOAT_ENCODING: &str = "ieee754-binary64-u64";
fn finite_from_bits(bits: u64) -> f64 {
    let value = f64::from_bits(bits);
    assert!(value.is_finite(), "nonfinite scientific input");
    value
}
fn decode_points(bits: Vec<[u64; 3]>) -> Vec<Point> {
    bits.into_iter().map(|p| p.map(finite_from_bits)).collect()
}
fn encode_points(points: &[Point]) -> Vec<[u64; 3]> {
    points
        .iter()
        .map(|p| {
            assert!(p.iter().all(|v| v.is_finite()));
            p.map(f64::to_bits)
        })
        .collect()
}

fn sorted_pairs(mut pairs: [Vec<[usize; 2]>; 3]) -> [Vec<[usize; 2]>; 3] {
    for axis in &mut pairs {
        axis.sort_unstable();
    }
    pairs
}

#[derive(Deserialize)]
struct FixtureBits {
    name: String,
    periodic: [bool; 3],
    #[serde(default)]
    coefficients_bits: Vec<[u64; 5]>,
    #[serde(default)]
    threshold_bits: u64,
    #[serde(default)]
    knots_bits: Vec<u64>,
    #[serde(default)]
    widths_bits: Vec<u64>,
}
impl From<FixtureBits> for Fixture {
    fn from(bits: FixtureBits) -> Self {
        Self {
            name: bits.name,
            periodic: bits.periodic,
            coefficients: bits
                .coefficients_bits
                .into_iter()
                .map(|row| row.map(finite_from_bits))
                .collect(),
            threshold: finite_from_bits(bits.threshold_bits),
            knots: bits.knots_bits.into_iter().map(finite_from_bits).collect(),
            widths: bits.widths_bits.into_iter().map(finite_from_bits).collect(),
        }
    }
}

impl Fixture {
    fn evaluate(&self, mut p: Point, gradient: bool) -> (f64, Point) {
        for (axis, x) in p.iter_mut().enumerate() {
            // Match implicit::Context: periodic values use the canonical trace,
            // while gradients retain the physical closed-cell side. Equal cap
            // values do not imply equal normal derivatives for the lattice.
            if !gradient && self.periodic[axis] && *x == 0.5 {
                *x = -0.5;
            }
        }
        if self.name == "plane" {
            return (p[2] - 0.2, [0., 0., 1.]);
        }
        if self.name == "sphere" {
            return (dot(p, p) - 0.09, p.map(|v| 2. * v));
        }
        if self.name.starts_with("spinodoid") {
            let (mut value, mut normal) = (self.threshold, [0.; 3]);
            for &[a, b, c, amplitude, phase] in &self.coefficients {
                let angle = TAU * (a * (p[0] + 0.5) + b * (p[1] + 0.5) + c * (p[2] + 0.5)) + phase;
                value -= amplitude * angle.cos();
                if gradient {
                    let factor = amplitude * TAU * angle.sin();
                    for (g, k) in normal.iter_mut().zip([a, b, c]) {
                        *g += factor * k;
                    }
                }
            }
            return (value, normal);
        }
        if self.name == "smooth-lattice" {
            let (mut value, mut normal) = (0., [0.; 3]);
            let mut first = true;
            // Frozen ordered union: x outermost, then y, then z. Tie derivatives
            // match the existing symmetric compiled field; no union reordering.
            for x in [-0.5, 0.5] {
                for y in [-0.5, 0.5] {
                    for z in [-0.5, 0.5] {
                        let corner = [x, y, z];
                        let t = (dot(p, corner) / 0.75).clamp(0., 1.);
                        let delta = std::array::from_fn(|i| p[i] - t * corner[i]);
                        let radius = norm(delta);
                        let next = if radius > 0. {
                            delta.map(|v| v / radius)
                        } else {
                            [0.; 3]
                        };
                        if first {
                            value = radius;
                            normal = next;
                            first = false;
                            continue;
                        }
                        let difference = value - radius;
                        let blend = ((0.08 - difference.abs()) / 0.08).max(0.);
                        let sign = if difference > 0. {
                            1.
                        } else if difference < 0. {
                            -1.
                        } else {
                            0.
                        };
                        let fraction = 0.5 - 0.5 * sign * (1. - blend * blend);
                        normal = std::array::from_fn(|i| {
                            fraction * normal[i] + (1. - fraction) * next[i]
                        });
                        value =
                            (value + radius - difference.abs()) / 2. - blend.powi(3) * 0.08 / 6.;
                    }
                }
            }
            value -= 0.12;
            return if value > 0.02 {
                (0.02, [0.; 3])
            } else {
                (value, normal)
            };
        }
        let [x, y, z] = p.map(|v| v * TAU);
        let (sx, cx, sy, cy, sz, cz) = (x.sin(), x.cos(), y.sin(), y.cos(), z.sin(), z.cos());
        if self.name == "schwarz-solid" {
            return (cx + cy + cz, [-TAU * sx, -TAU * sy, -TAU * sz]);
        }
        if self.name == "schwarz-d-solid" {
            return (
                sx * sy * sz + sx * cy * cz + cx * sy * cz + cx * cy * sz,
                [
                    cx * sy * sz + cx * cy * cz - sx * sy * cz - sx * cy * sz,
                    sx * cy * sz - sx * sy * cz + cx * cy * cz - cx * sy * sz,
                    sx * sy * cz - sx * cy * sz - cx * sy * sz + cx * cy * cz,
                ]
                .map(|v| TAU * v),
            );
        }
        let raw = sx * cy + sy * cz + sz * cx;
        let (half, slope) = match self.name.as_str() {
            "gyroid" => (0.55, 0.),
            "graded-gyroid-linear" => (0.55 * (1. + 0.5 * p[2]), 0.275),
            "graded-gyroid-periodic" => (0.55 * (1. + 0.25 * cz), -0.1375 * TAU * sz),
            "graded-gyroid-density-linear" => {
                if p[2] <= self.knots[0] {
                    (self.widths[0], 0.)
                } else if p[2] >= *self.knots.last().unwrap() {
                    (*self.widths.last().unwrap(), 0.)
                } else {
                    let i = self.knots.partition_point(|&k| k <= p[2]) - 1;
                    let slope =
                        (self.widths[i + 1] - self.widths[i]) / (self.knots[i + 1] - self.knots[i]);
                    (self.widths[i] + slope * (p[2] - self.knots[i]), slope)
                }
            }
            _ => panic!("unknown frozen fixture {}", self.name),
        };
        let mut normal =
            [cx * cy - sz * sx, cy * cz - sx * sy, cz * cx - sy * sz].map(|v| TAU * v / half);
        normal[2] -= raw * slope / (half * half);
        (raw / half, normal)
    }
}

impl Geometry for Fixture {
    fn value(&self, p: Point) -> f64 {
        self.evaluate(p, false).0
    }
    fn gradient(&self, p: Point) -> Point {
        self.evaluate(p, true).1
    }
    fn rank(&self, _: Point) -> [i64; 3] {
        [0; 3]
    }
    fn bounds(&self) -> [Point; 2] {
        [[-0.5; 3], [0.5; 3]]
    }
    fn band(&self) -> bool {
        self.name.contains("gyroid")
    }
}

struct Rule {
    shape: [f64; 6],
    du: [f64; 6],
    dv: [f64; 6],
    weight: f64,
}
fn triangle(a: f64, b: f64, c: f64, weight: f64) -> Rule {
    Rule {
        shape: [
            a * (2. * a - 1.),
            b * (2. * b - 1.),
            c * (2. * c - 1.),
            4. * a * b,
            4. * b * c,
            4. * c * a,
        ],
        du: [1. - 4. * a, 4. * b - 1., 0., 4. * (a - b), 4. * c, -4. * c],
        dv: [1. - 4. * a, 0., 4. * c - 1., -4. * b, 4. * b, 4. * (a - c)],
        weight,
    }
}
fn gauss(degree: usize) -> Vec<Rule> {
    let (nodes, weights): (&[f64], &[f64]) = if degree == 6 {
        (
            &[
                0.033765242898424,
                0.169395306766868,
                0.380690406958402,
                0.619309593041598,
                0.830604693233132,
                0.966234757101576,
            ],
            &[
                0.085662246189585,
                0.180380786524069,
                0.233956967286346,
                0.233956967286346,
                0.180380786524069,
                0.085662246189585,
            ],
        )
    } else {
        assert_eq!(degree, 8);
        (
            &[
                0.019855071751232,
                0.101666761293187,
                0.237233795041836,
                0.408282678752175,
                0.591717321247825,
                0.762766204958164,
                0.898333238706813,
                0.980144928248768,
            ],
            &[
                0.050614268145188,
                0.111190517226687,
                0.156853322938944,
                0.181341891689181,
                0.181341891689181,
                0.156853322938944,
                0.111190517226687,
                0.050614268145188,
            ],
        )
    };
    let mut rules = Vec::new();
    for (&u, &wu) in nodes.iter().zip(weights) {
        for (&v, &wv) in nodes.iter().zip(weights) {
            rules.push(triangle(
                1. - u - (1. - u) * v,
                u,
                (1. - u) * v,
                wu * wv * (1. - u),
            ));
        }
    }
    rules
}

struct Observation {
    residual: f64,
    normal: Point,
    weight: f64,
    cap_residual: f64,
}
fn observe(
    faces: &[Face],
    points: &[Point],
    rule: &[Rule],
    g: &impl Geometry,
) -> Result<Vec<Observation>, String> {
    let mut out = Vec::with_capacity(faces.len() * rule.len());
    for face in faces {
        for sample in rule {
            let mut p: Point = std::array::from_fn(|a| {
                (0..6)
                    .map(|i| sample.shape[i] * points[face.nodes[i]][a])
                    .sum()
            });
            let mut locked = [false; 3];
            for a in 0..3 {
                for bound in g.bounds() {
                    if (p[a] - bound[a]).abs() <= 1e-12 * 3_f64.sqrt() {
                        p[a] = bound[a];
                        locked[a] = true;
                    }
                }
            }
            let value = g.value(p) - face.level;
            let gradient = g.gradient(p);
            let length = norm(gradient);
            if !length.is_finite() || length <= 1e-10 || !value.is_finite() {
                return Err("stationary or nonfinite wall sample".into());
            }
            let allowed = norm(std::array::from_fn(|a| {
                if locked[a] { 0. } else { gradient[a] }
            }));
            let cap_residual = if allowed > 1e-14 {
                value / allowed
            } else if value.abs() <= 1e-12 {
                0.
            } else {
                return Err("stationary off-wall cap sample".into());
            };
            let sign = if g.band() && face.level < 0. { -1. } else { 1. };
            let weight = if sample.weight == 0. {
                0.
            } else {
                let du = std::array::from_fn(|a| {
                    (0..6)
                        .map(|i| sample.du[i] * points[face.nodes[i]][a])
                        .sum()
                });
                let dv = std::array::from_fn(|a| {
                    (0..6)
                        .map(|i| sample.dv[i] * points[face.nodes[i]][a])
                        .sum()
                });
                sample.weight * norm(cross(du, dv))
            };
            out.push(Observation {
                residual: sign * value / length,
                normal: gradient.map(|v| sign * v / length),
                weight,
                cap_residual,
            });
        }
    }
    Ok(out)
}
fn bias(values: &[Observation]) -> f64 {
    values.iter().map(|s| s.residual * s.weight).sum()
}
fn envelope(values: &[Observation]) -> (f64, f64) {
    values.iter().fold((0_f64, 0_f64), |(lo, hi), s| {
        (lo.min(s.residual), hi.max(s.residual))
    })
}
fn within_envelope(values: &[Observation], bounds: (f64, f64)) -> bool {
    values
        .iter()
        .all(|s| s.residual >= bounds.0 - 1e-12 && s.residual <= bounds.1 + 1e-12)
}

#[derive(Default, Serialize)]
struct Timing {
    setup_s: f64,
    accumulation_s: f64,
    guards_s: f64,
    jacobians_s: f64,
    total_s: f64,
}
#[derive(Serialize)]
struct Attempt {
    fraction: f64,
    bias: f64,
    signed_guard: bool,
    cap_guard: bool,
    jacobian_guard: bool,
    minimum_ratio: f64,
    accepted: bool,
}
#[derive(Serialize)]
struct Stats {
    timing: Timing,
    before6: f64,
    before8: f64,
    after6: f64,
    quadrature_stable: bool,
    attempts: Vec<Attempt>,
    wall_faces: usize,
    moved_nodes: usize,
    affected_tets: usize,
    budget_ratio: f64,
    linear_unclipped_bias: f64,
    linear_clipped_bias: f64,
    cap_max_before: f64,
    cap_max_after: f64,
    periodic_error_before: f64,
    periodic_error_after: f64,
    periodic_offset_change: f64,
}

fn periodic_error(points: &[Point], pairs: &[Vec<[usize; 2]>; 3], bounds: [Point; 2]) -> f64 {
    let mut error = 0_f64;
    for (axis, rows) in pairs.iter().enumerate() {
        for &[a, b] in rows {
            for (coordinate, (&lo, &hi)) in points[a].iter().zip(&points[b]).enumerate() {
                let shift = if coordinate == axis {
                    bounds[1][axis] - bounds[0][axis]
                } else {
                    0.
                };
                error = error.max((hi - lo - shift).abs());
            }
        }
    }
    error
}

fn paired_offset_change(before: &[Point], after: &[Point], pairs: &[Vec<[usize; 2]>; 3]) -> f64 {
    let mut change = 0_f64;
    for rows in pairs {
        for &[a, b] in rows {
            for coordinate in 0..3 {
                change = change.max(
                    ((after[b][coordinate] - after[a][coordinate])
                        - (before[b][coordinate] - before[a][coordinate]))
                        .abs(),
                );
            }
        }
    }
    change
}

fn translated_trial(
    current: &[Point],
    mids: &[usize],
    classes: &[usize],
    direction: &[Point],
    step: &[f64],
    fraction: f64,
    g: &impl Geometry,
) -> Vec<Point> {
    let mut trial = current.to_vec();
    for &m in mids {
        let group = classes[m];
        let locked = g.locked(current[m]);
        for a in 0..3 {
            let delta = fraction * step[group] * direction[group][a];
            if !locked[a] && delta != 0. {
                // Preserve each fitted alias's existing coordinates. The common
                // requested displacement preserves its translation in real
                // arithmetic; floating-point addition may round by an ulp.
                // Skipping zero delta also preserves signed zero exactly.
                trial[m][a] = current[m][a] + delta;
            }
        }
    }
    trial
}

fn correct(
    mesh: &mut Quadratic,
    original: &[Point],
    g: &impl Geometry,
    periodic: [bool; 3],
) -> Result<Stats, String> {
    let total = Instant::now();
    let mut timing = Timing::default();
    let stage = Instant::now();
    let current = mesh.points.clone();
    let n = current.len();
    assert_eq!(original.len(), n);
    let pairs = node_pairs(&current, &mesh.surface, g.bounds(), periodic)?;
    let periodic_error_before = periodic_error(&current, &pairs, g.bounds());
    // The unchanged independent admission tolerance for these unit-cell cases.
    assert!(periodic_error_before <= 1e-10);
    let mut parent: Vec<_> = (0..n).collect();
    for axis in &pairs {
        for &[a, b] in axis {
            let a = root(&mut parent, a);
            let b = root(&mut parent, b);
            parent[a.max(b)] = a.min(b);
        }
    }
    let classes: Vec<_> = (0..n).map(|i| root(&mut parent, i)).collect();
    let faces: Vec<_> = mesh
        .surface
        .iter()
        .filter(|f| !g.cap([current[f[0]], current[f[1]], current[f[2]]]))
        .map(|&nodes| {
            let center =
                std::array::from_fn(|a| (0..3).map(|i| current[nodes[i]][a]).sum::<f64>() / 3.);
            Face {
                nodes,
                level: g.level(g.value(center), 1.),
                area: 0.,
            }
        })
        .collect();
    let mut mids: Vec<_> = faces
        .iter()
        .flat_map(|f| f.nodes[3..].iter().copied())
        .collect();
    mids.sort_unstable();
    mids.dedup();
    // Pre-fit and fitted aliases can differ by a few ulps while passing the
    // existing periodic gate. Each retains its own position and movement budget.
    // Their locked axes must nevertheless describe the same periodic class.
    for &m in &mids {
        assert_eq!(g.locked(current[m]), g.locked(current[classes[m]]));
    }
    let mut direction = vec![[0.; 3]; n];
    let mut limit = vec![f64::INFINITY; n];
    let mut used = vec![0_f64; n];
    for &m in &mids {
        let gradient = g.gradient(current[m]);
        let locked = g.locked(current[m]);
        for a in 0..3 {
            if !locked[a] {
                direction[classes[m]][a] += gradient[a];
            }
        }
        used[classes[m]] = used[classes[m]].max(norm(sub(current[m], original[m])));
    }
    for d in &mut direction {
        let length = norm(*d);
        if length > 0. {
            *d = d.map(|v| v / length);
        }
    }
    for face in &faces {
        for (k, [a, b]) in [[0, 1], [1, 2], [2, 0]].into_iter().enumerate() {
            let group = classes[face.nodes[3 + k]];
            limit[group] = limit[group]
                .min(0.02 * norm(sub(original[face.nodes[a]], original[face.nodes[b]])));
        }
    }
    let remaining: Vec<_> = (0..n).map(|i| (limit[i] - used[i]).max(0.)).collect();
    let rule6 = gauss(6);
    let rule8 = gauss(8);
    let grid: Vec<_> = basis()
        .into_iter()
        .map(|shape| Rule {
            shape,
            du: [0.; 6],
            dv: [0.; 6],
            weight: 0.,
        })
        .collect();
    timing.setup_s = stage.elapsed().as_secs_f64();
    let stage = Instant::now();
    let data = observe(&faces, &current, &rule6, g)?;
    let before6 = bias(&data);
    let before8 = bias(&observe(&faces, &current, &rule8, g)?);
    let stable = (before6 - before8).abs() <= 0.1 * before8.abs();
    let mut gradient = vec![0.; n];
    let mut diagonal = vec![0_f64; n];
    for (face, values) in faces.iter().zip(data.chunks_exact(rule6.len())) {
        for (s, q) in values.iter().zip(&rule6) {
            for k in 0..3 {
                let group = classes[face.nodes[3 + k]];
                let coefficient = q.shape[3 + k] * dot(s.normal, direction[group]);
                gradient[group] += s.weight * coefficient;
                diagonal[group] += s.weight * coefficient * coefficient;
            }
        }
    }
    let floor = diagonal.iter().copied().fold(0., f64::max) * 1e-12;
    let denominator: f64 = (0..n)
        .filter(|&i| diagonal[i] > 0.)
        .map(|i| gradient[i] * gradient[i] / diagonal[i].max(floor))
        .sum();
    let mut step: Vec<_> = (0..n)
        .map(|i| {
            if diagonal[i] > 0. {
                -before6 * gradient[i] / diagonal[i].max(floor) / denominator.max(1e-30)
            } else {
                0.
            }
        })
        .collect();
    let linear_unclipped_bias =
        before6 + gradient.iter().zip(&step).map(|(a, b)| a * b).sum::<f64>();
    timing.accumulation_s += stage.elapsed().as_secs_f64();
    let stage = Instant::now();
    let initial = observe(&faces, &current, &grid, g)?;
    let bounds: Vec<_> = initial.chunks_exact(grid.len()).map(envelope).collect();
    let cap_before = initial
        .iter()
        .map(|s| s.cap_residual.abs())
        .fold(0., f64::max);
    let mut lower: Vec<_> = remaining.iter().map(|v| -v).collect();
    let mut upper = remaining;
    for (f, face) in faces.iter().enumerate() {
        for (i, q) in grid.iter().enumerate() {
            let s = &initial[f * grid.len() + i];
            for k in 0..3 {
                let group = classes[face.nodes[3 + k]];
                let coefficient = q.shape[3 + k] * dot(s.normal, direction[group]);
                if coefficient.abs() > 1e-14 {
                    let a = (bounds[f].0 - s.residual) / (3. * coefficient);
                    let b = (bounds[f].1 - s.residual) / (3. * coefficient);
                    lower[group] = lower[group].max(a.min(b));
                    upper[group] = upper[group].min(a.max(b));
                }
            }
        }
    }
    for i in 0..n {
        assert!(lower[i] <= 1e-12 && upper[i] >= -1e-12);
        step[i] = step[i].clamp(lower[i].min(0.), upper[i].max(0.));
    }
    let linear_clipped_bias = before6 + gradient.iter().zip(&step).map(|(a, b)| a * b).sum::<f64>();
    let affected: Vec<_> = mesh
        .tets
        .iter()
        .filter(|tet| tet[4..].iter().any(|&m| step[classes[m]] != 0.))
        .collect();
    let probes = validity_samples();
    let mut attempts = Vec::new();
    let (mut after6, mut cap_after) = (before6, cap_before);
    timing.guards_s += stage.elapsed().as_secs_f64();
    if stable && step.iter().any(|&s| s != 0.) {
        for fraction in [1., 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125] {
            let stage = Instant::now();
            let trial = translated_trial(&current, &mids, &classes, &direction, &step, fraction, g);
            timing.guards_s += stage.elapsed().as_secs_f64();
            let stage = Instant::now();
            let candidate = bias(&observe(&faces, &trial, &rule6, g)?);
            timing.accumulation_s += stage.elapsed().as_secs_f64();
            let stage = Instant::now();
            let values = observe(&faces, &trial, &grid, g)?;
            let signed_guard = values
                .chunks_exact(grid.len())
                .zip(&bounds)
                .all(|(v, &b)| within_envelope(v, b));
            let cap_max = values
                .iter()
                .map(|s| s.cap_residual.abs())
                .fold(0., f64::max);
            let cap_guard = cap_max <= cap_before + 1e-12 * 3_f64.sqrt();
            timing.guards_s += stage.elapsed().as_secs_f64();
            let stage = Instant::now();
            let mut ratio = 1_f64;
            let jacobian_guard = if signed_guard && cap_guard && candidate.abs() < before6.abs() {
                affected.iter().all(|tet| {
                    let coordinates = tet.map(|m| trial[m]);
                    let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
                    for &p in &probes {
                        let d = jacobian(&coordinates, p);
                        if !d.is_finite() {
                            return false;
                        }
                        low = low.min(d);
                        high = high.max(d);
                    }
                    if high > 0. {
                        ratio = ratio.min(low / high);
                    } else {
                        ratio = 0.;
                    }
                    low > 1e-3 * high && high > 0.
                })
            } else {
                false
            };
            timing.jacobians_s += stage.elapsed().as_secs_f64();
            let accepted =
                signed_guard && cap_guard && jacobian_guard && candidate.abs() < before6.abs();
            attempts.push(Attempt {
                fraction,
                bias: candidate,
                signed_guard,
                cap_guard,
                jacobian_guard,
                minimum_ratio: ratio,
                accepted,
            });
            if accepted {
                mesh.points = trial;
                after6 = candidate;
                cap_after = cap_max;
                break;
            }
        }
    }
    let stage = Instant::now();
    let mut budget_ratio = 0_f64;
    for &m in &mids {
        budget_ratio = budget_ratio.max(norm(sub(mesh.points[m], original[m])) / limit[classes[m]]);
    }
    assert!(budget_ratio <= 1. + 1e-10);
    for f in &mesh.surface {
        for &m in &f[..3] {
            assert_eq!(
                mesh.points[m].map(f64::to_bits),
                current[m].map(f64::to_bits)
            );
        }
    }
    for tet in &mesh.tets {
        for &m in &tet[..4] {
            assert_eq!(
                mesh.points[m].map(f64::to_bits),
                current[m].map(f64::to_bits)
            );
        }
    }
    let mut wall_midpoint = vec![false; n];
    for &m in &mids {
        wall_midpoint[m] = true;
    }
    for m in 0..n {
        if !wall_midpoint[m] {
            assert_eq!(
                mesh.points[m].map(f64::to_bits),
                current[m].map(f64::to_bits)
            );
        }
        for a in 0..3 {
            if g.locked(current[m])[a] {
                assert_eq!(mesh.points[m][a].to_bits(), current[m][a].to_bits());
            }
        }
    }
    assert_eq!(
        sorted_pairs(node_pairs(
            &mesh.points,
            &mesh.surface,
            g.bounds(),
            periodic
        )?),
        sorted_pairs(pairs.clone()),
        "periodic pairing changed, independent of coordinate-dependent list order"
    );
    let periodic_error_after = periodic_error(&mesh.points, &pairs, g.bounds());
    assert!(periodic_error_after <= 1e-10);
    let periodic_offset_change = paired_offset_change(&current, &mesh.points, &pairs);
    let moved_nodes = mids
        .iter()
        .filter(|&&m| mesh.points[m] != current[m])
        .count();
    timing.guards_s += stage.elapsed().as_secs_f64();
    timing.total_s = total.elapsed().as_secs_f64();
    Ok(Stats {
        timing,
        before6,
        before8,
        after6,
        quadrature_stable: stable,
        attempts,
        wall_faces: faces.len(),
        moved_nodes,
        affected_tets: affected.len(),
        budget_ratio,
        linear_unclipped_bias,
        linear_clipped_bias,
        cap_max_before: cap_before,
        cap_max_after: cap_after,
        periodic_error_before,
        periodic_error_after,
        periodic_offset_change,
    })
}

struct Input {
    fixture: Fixture,
    points: Vec<Point>,
    original: Vec<Point>,
    tetrahedra: Vec<[usize; 10]>,
    surface: Vec<[usize; 6]>,
    probes: Vec<Point>,
    values: Vec<f64>,
    gradients: Vec<Point>,
}

#[derive(Deserialize)]
struct InputBits {
    float_encoding: String,
    fixture: FixtureBits,
    points_bits: Vec<[u64; 3]>,
    original_bits: Vec<[u64; 3]>,
    tetrahedra: Vec<[usize; 10]>,
    surface: Vec<[usize; 6]>,
    probes_bits: Vec<[u64; 3]>,
    values_bits: Vec<u64>,
    gradients_bits: Vec<[u64; 3]>,
}
impl From<InputBits> for Input {
    fn from(bits: InputBits) -> Self {
        assert_eq!(bits.float_encoding, FLOAT_ENCODING);
        Self {
            fixture: bits.fixture.into(),
            points: decode_points(bits.points_bits),
            original: decode_points(bits.original_bits),
            tetrahedra: bits.tetrahedra,
            surface: bits.surface,
            probes: decode_points(bits.probes_bits),
            values: bits.values_bits.into_iter().map(finite_from_bits).collect(),
            gradients: decode_points(bits.gradients_bits),
        }
    }
}
fn mesh(input: &Input) -> Quadratic {
    Quadratic {
        points: input.points.clone(),
        tets: input.tetrahedra.clone(),
        surface: input.surface.clone(),
        curved_edges: 0,
        straightened_edges: 0,
        partially_curved_edges: 0,
        minimum_jacobian_ratio: 1.,
        fit: None,
    }
}
fn progress(fraction: f64, message: &str) {
    if let Ok(part) = std::env::var("MESHERS_VOLUME_PROGRESS_PART") {
        let code = Command::new("pueue-progress")
            .args([
                "--part",
                &part,
                &fraction.to_string(),
                message,
                "--value",
                &format!("repeats_done={}", 4. * fraction),
                "--value",
                "threads=1",
            ])
            .status()
            .unwrap();
        assert!(code.success(), "standard progress helper failed");
    }
}

#[test]
#[ignore = "bounded saved-mesh experiment; parent owns the shared queue"]
fn volume_bias_payoff() {
    let input_bits: InputBits = serde_json::from_reader(BufReader::new(
        File::open(std::env::var("MESHERS_VOLUME_INPUT").unwrap()).unwrap(),
    ))
    .unwrap();
    let input: Input = input_bits.into();
    assert!(!input.probes.is_empty());
    assert_eq!(input.probes.len(), input.values.len());
    assert_eq!(input.probes.len(), input.gradients.len());
    assert_eq!(input.points.len(), input.original.len());
    for tet in &input.tetrahedra {
        for &m in &tet[..4] {
            assert_eq!(input.points[m], input.original[m]);
        }
    }
    let mut value_error = 0_f64;
    let mut gradient_error = 0_f64;
    for ((&p, &expected), &normal) in input.probes.iter().zip(&input.values).zip(&input.gradients) {
        let (value, gradient) = input.fixture.evaluate(p, true);
        value_error = value_error.max((value - expected).abs());
        gradient_error = gradient_error.max(norm(sub(gradient, normal)));
        assert!((value - expected).abs() <= 1e-11 * (1. + expected.abs()));
        assert!(
            norm(sub(gradient, normal)) <= 2e-9 * (1. + norm(normal)),
            "gradient mismatch at {p:?}: {gradient:?} vs {normal:?}"
        );
    }
    progress(0., "field identities passed; excluded warmup");
    let mut repeats = Vec::new();
    let mut exported: Option<Vec<Point>> = None;
    for repeat in 0..4 {
        let mut candidate = mesh(&input);
        let stats = correct(
            &mut candidate,
            &input.original,
            &input.fixture,
            input.fixture.periodic,
        )
        .unwrap();
        assert!(
            [
                stats.timing.setup_s,
                stats.timing.accumulation_s,
                stats.timing.guards_s,
                stats.timing.jacobians_s,
                stats.timing.total_s,
            ]
            .iter()
            .all(|value| value.is_finite() && *value >= 0.)
        );
        if let Some(points) = &exported {
            assert_eq!(encode_points(&candidate.points), encode_points(points));
        }
        exported = Some(candidate.points);
        if repeat > 0 {
            repeats.push(stats);
        }
        progress(
            (repeat + 1) as f64 / 4.,
            if repeat == 0 {
                "warmup complete"
            } else {
                "timed native repeat complete"
            },
        );
    }
    let output = serde_json::json!({"case":input.fixture.name,"points_bits":encode_points(&exported.unwrap()),
        "float_encoding":FLOAT_ENCODING,"repeats":repeats,
        "identity":{"probes":input.probes.len(),"maximum_value_error":value_error,"maximum_gradient_error":gradient_error},
        "warmups":1,"timed_repetitions":3,"all_four_outputs_identical":true,
        "additional_guards_vs_python":"Existing native cap-constrained maximum and sampled min/max Jacobian ratio >0.001 are added to the Python signed-envelope policy; identical candidates are not assumed.",
        "periodic_motion":"Common class displacement added to each fitted alias's own position; zero components and locks unchanged bitwise. Existing offsets are preserved in real arithmetic, with recorded addition-rounding changes. No coordinate reconciliation; unchanged1e-10 unit-cell coordinate gate plus native pair/facet/class checks.",
        "evaluator":"direct Rust analytic fixture adapter; not production JIT workflow timing"});
    serde_json::to_writer(
        BufWriter::new(File::create(std::env::var("MESHERS_VOLUME_OUTPUT").unwrap()).unwrap()),
        &output,
    )
    .unwrap();
}

fn transport_controls() {
    use std::io::Cursor;
    // Job780: node9009/x, decimal -0.23958333333333334, was read one ulp
    // below Python's frozen input. The global serde parser remains untouched.
    let authoritative = 13821171956431104683_u64;
    let historical_decimal_parse = 13821171956431104684_u64;
    let legacy: f64 = serde_json::from_reader(Cursor::new(b"-0.23958333333333334")).unwrap();
    assert!([authoritative, historical_decimal_parse].contains(&legacy.to_bits()));
    let bits = [
        0,
        0x8000_0000_0000_0000,
        1,
        0x000f_ffff_ffff_ffff,
        0x0010_0000_0000_0000,
        0x7fef_ffff_ffff_ffff,
        0xffef_ffff_ffff_ffff,
        authoritative,
    ];
    let bytes = serde_json::to_vec(&bits).unwrap();
    let parsed: Vec<u64> = serde_json::from_reader(Cursor::new(bytes)).unwrap();
    let decoded: Vec<f64> = parsed.into_iter().map(finite_from_bits).collect();
    assert_eq!(
        decoded.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        bits
    );
    assert_ne!(decoded.last().unwrap().to_bits(), historical_decimal_parse);
    let points = vec![
        [decoded[0], decoded[1], decoded[7]],
        [decoded[2], decoded[5], decoded[6]],
    ];
    let bytes = serde_json::to_vec(&encode_points(&points)).unwrap();
    let parsed: Vec<[u64; 3]> = serde_json::from_reader(Cursor::new(bytes)).unwrap();
    assert_eq!(
        encode_points(&decode_points(parsed)),
        encode_points(&points)
    );
}

fn periodic_motion_controls(g: &impl Geometry) {
    // Frozen Gyroid current aliases: root9859 and node10907 differ by one ulp.
    // Their pre-fit arrays also differ. Moving a root coordinate onto its alias
    // would erase this real input difference, including for a zero step.
    let current = vec![
        [-0.5, -0.0, f64::from_bits(4599090790396799705)],
        [0.5, -0.0, f64::from_bits(4599090790396799704)],
        [-0.5, -0.0, 0.125],
        [0.1, -0.0, 0.2],
    ];
    let mids = [0, 1, 3];
    let classes = [0, 0, 2, 3];
    let directions = [[1., 0., 1.], [0.; 3], [0.; 3], [1., 1., 1.]];
    let zero = translated_trial(&current, &mids, &classes, &directions, &[0.; 4], 1., g);
    assert_eq!(encode_points(&zero), encode_points(&current));
    let pairs = [vec![[0, 1]], vec![], vec![]];
    let step = [2_f64.powi(-20), 0., 0., 0.];
    let moved = translated_trial(&current, &mids, &classes, &directions, &step, 1., g);
    assert_eq!(
        (moved[1][2] - moved[0][2]).to_bits(),
        (current[1][2] - current[0][2]).to_bits()
    );
    assert_ne!(moved[1][2].to_bits(), moved[0][2].to_bits());
    for m in 0..4 {
        for a in 0..3 {
            if a != 2 || m >= 2 {
                assert_eq!(moved[m][a].to_bits(), current[m][a].to_bits());
            }
        }
    }
    assert_eq!(paired_offset_change(&current, &moved, &pairs), 0.);
    let rounded = translated_trial(
        &current,
        &mids,
        &classes,
        &directions,
        &[1.23456789e-6, 0., 0., 0.],
        1.,
        g,
    );
    assert!(paired_offset_change(&current, &rounded, &pairs) <= 8. * f64::EPSILON);
    assert!(periodic_error(&rounded, &pairs, g.bounds()) <= 1e-10);
}

fn lattice_cap_controls(mut g: Fixture) {
    g.periodic = [true; 3];
    // Six original compiled-field probes from the frozen failed study, including
    // the positive-x cap that exposed accidental gradient canonicalization.
    // native-input.json SHA256: 65fa512f29236be4cf95db1e64dfc271051c31219cf0dea0d15ab6b5f29d8caa
    let controls: [([u64; 3], u64, [u64; 3]); 6] = [
        (
            [
                13826050856027422720,
                13824149630989893632,
                4600016382552178688,
            ],
            13801643307645574208,
            [
                13828612655006171804,
                4596227833288637540,
                13826847701756646278,
            ],
        ),
        (
            [
                4602678819172646912,
                4600020601048924160,
                13824967167730253824,
            ],
            13803211501767418048,
            [
                4604118576313143262,
                13828270228154903452,
                13814273377413815784,
            ],
        ),
        (
            [
                13823750935937024000,
                13826050856027422720,
                13823195263990235136,
            ],
            13786297933658096640,
            [
                4598563453706143594,
                13828776706095552898,
                4602958752346658253,
            ],
        ),
        (
            [
                4599934358827040768,
                4602678819172646912,
                4601090902201466880,
            ],
            13801345722656417216,
            [
                13827557914138032000,
                4604850536748709784,
                13813353733241569482,
            ],
        ),
        (
            [
                13823272564660633600,
                4602004259192438784,
                13826050856027422720,
            ],
            13797213648467566208,
            [
                4605327593636451832,
                4597543211521861556,
                13826605832683396106,
            ],
        ),
        (
            [
                13824149630989893632,
                4600016382552178688,
                4602678819172646912,
            ],
            13801643307645574208,
            [
                4596227833288637542,
                13826847701756646278,
                4605240618151395998,
            ],
        ),
    ];
    for (point, value, gradient) in controls {
        let p = point.map(finite_from_bits);
        let expected = gradient.map(finite_from_bits);
        let value = finite_from_bits(value);
        assert!((g.value(p) - value).abs() <= 1e-11 * (1. + value.abs()));
        assert!(norm(sub(g.gradient(p), expected)) <= 2e-9 * (1. + norm(expected)));
    }
}

#[test]
fn volume_probe_controls() {
    transport_controls();
    // Pair lists follow spatial order; tangential motion may reorder them.
    // Compare every directed pair including multiplicity, not just pair counts.
    let pairs = [vec![[0, 2], [1, 3]], vec![], vec![]];
    assert_eq!(
        sorted_pairs(pairs.clone()),
        sorted_pairs([vec![[1, 3], [0, 2]], vec![], vec![]])
    );
    assert_ne!(
        sorted_pairs(pairs.clone()),
        sorted_pairs([vec![[0, 3], [1, 2]], vec![], vec![]])
    );
    assert_ne!(
        sorted_pairs(pairs),
        sorted_pairs([vec![[0, 2], [0, 2]], vec![], vec![]])
    );
    let make = |name: &str| Fixture {
        name: name.into(),
        periodic: [false; 3],
        coefficients: vec![],
        threshold: 0.,
        knots: vec![],
        widths: vec![],
    };
    periodic_motion_controls(&make("plane"));
    lattice_cap_controls(make("smooth-lattice"));
    for degree in [6, 8] {
        assert!((gauss(degree).iter().map(|r| r.weight).sum::<f64>() - 0.5).abs() < 1e-14);
    }
    let observed = |r| Observation {
        residual: r,
        normal: [0.; 3],
        weight: 0.,
        cap_residual: 0.,
    };
    let bounds = envelope(&[observed(-0.7), observed(1.)]);
    assert!(maximum(&[-0.8, 0.9]) < maximum(&[-0.7, 1.]));
    assert!(!within_envelope(&[observed(-0.8), observed(0.9)], bounds));
    assert!(within_envelope(&[observed(-0.6), observed(0.8)], bounds));
    let mut points = vec![[0.2, 0., 1.], [0., 0.2, 1.], [-0.2, -0.2, 1.]];
    for p in &mut points {
        let length = norm(*p);
        *p = p.map(|v| v * 0.3 / length);
    }
    for [a, b] in [[0, 1], [1, 2], [2, 0]] {
        let p = std::array::from_fn(|i| points[a][i] + points[b][i]);
        points.push(p.map(|v| v * 0.3 / norm(p)));
    }
    let input = Input {
        fixture: make("sphere"),
        points: points.clone(),
        original: points.clone(),
        tetrahedra: vec![],
        surface: vec![[0, 1, 2, 3, 4, 5]],
        probes: vec![],
        values: vec![],
        gradients: vec![],
    };
    let mut curved = mesh(&input);
    // A small outward midpoint bias gives the envelope room to move toward the
    // exact sphere, as after the fitted saved meshes. Corners remain exact.
    for p in &mut curved.points[3..] {
        *p = p.map(|v| v * 1.0001);
    }
    let stats = correct(&mut curved, &points, &input.fixture, [false; 3]).unwrap();
    assert!(stats.after6.abs() < stats.before6.abs());
    assert_eq!(curved.points[..3], points[..3]);
    for p in &mut points {
        p[2] = 0.2;
    }
    let mut plane = mesh(&input);
    plane.points = points.clone();
    assert_eq!(
        correct(&mut plane, &points, &make("plane"), [false; 3])
            .unwrap()
            .moved_nodes,
        0
    );
    assert_eq!(plane.points, points);
}
