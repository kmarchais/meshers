//! Experimental simultaneous intersection of independently sampled sublevel sets.
//!
//! Each field is interpolated linearly inside a background tetrahedron. Cutting
//! retains every constraint's identity until all cuts are complete; only then
//! is the convex cell tetrahedralized. This preserves shared feature curves.
//! A supplied coordinate map is applied to the background before cutting, not
//! to the final thin elements. Optional joint background snapping enforces all
//! nearby constraints together; it is not a guarantee of element quality.
mod repair;
use super::{Options, ScalarField};
use crate::{Mesh, Point, Tet, determinant, dot, norm, quality, sub};
use std::collections::{HashMap, HashSet};

/// Experimental mesh plus active surface constraints at its vertices.
pub struct Output {
    pub mesh: Mesh,
    /// Bit j denotes membership in the piecewise-linear boundary of field j.
    pub constraints: Vec<u8>,
    /// Parameter coordinates, useful for identifying periodic end surfaces.
    pub parameters: Vec<Point>,
    /// Quality before and after the integrated improvement stage.
    pub quality: QualityDiagnostics,
}

/// Measured quality and accepted local operations, not solver qualification.
#[derive(Default, serde::Serialize)]
pub struct QualityDiagnostics {
    /// Background vertices whose joint snaps were backed out for topology.
    pub reverted_snap_vertices: usize,
    pub initial_minimum_quality: f64,
    pub initial_elements_below_01: usize,
    pub minimum_quality: f64,
    pub elements_below_01: usize,
    pub collapsed_vertices: usize,
    pub accepted_vertex_moves: usize,
    pub accepted_reconnections: usize,
}

struct Cutter<'a> {
    options: &'a Options,
    points: Vec<Point>,
    parameters: Vec<Point>,
    values: Vec<Vec<f64>>,
    masks: Vec<u8>,
    edges: HashMap<(usize, usize, usize), usize>,
}
impl Cutter<'_> {
    fn rank(&self, i: usize) -> [i64; 3] {
        std::array::from_fn(|a| {
            let [lo, hi] = self.options.bounds;
            let x = self.parameters[i][a];
            let x = if self.options.periodic[a] && x == hi[a] {
                lo[a]
            } else {
                x
            };
            ((x - lo[a]) / (hi[a] - lo[a]) * 1e11).round() as i64
        })
    }
    fn cut(&mut self, a: usize, b: usize, field: usize) -> usize {
        if self.values[a][field] == 0. {
            return a;
        }
        if self.values[b][field] == 0. {
            return b;
        }
        let (a, b) = if a < b { (a, b) } else { (b, a) };
        if let Some(&i) = self.edges.get(&(a, b, field)) {
            return i;
        }
        let va = self.values[a][field].abs();
        let vb = self.values[b][field].abs();
        let t = if va > vb {
            1. / (1. + vb / va)
        } else {
            (va / vb) / (1. + va / vb)
        };
        let interpolate = |x: Point, y: Point| std::array::from_fn(|k| x[k] + t * (y[k] - x[k]));
        let p = interpolate(self.points[a], self.points[b]);
        let u = interpolate(self.parameters[a], self.parameters[b]);
        let mut values: Vec<_> = self.values[a]
            .iter()
            .zip(&self.values[b])
            .map(|(x, y)| (1. - t) * x + t * y)
            .collect();
        let mut mask = (self.masks[a] & self.masks[b]) | (1 << field);
        for (j, value) in values.iter_mut().enumerate() {
            // Coincident cuts can differ by a few interpolation ulps. Keeping
            // those as distinct vertices creates zero-volume wedges and holes.
            let roundoff =
                32. * f64::EPSILON * self.values[a][j].abs().max(self.values[b][j].abs());
            if mask & (1 << j) != 0 || value.abs() <= roundoff {
                *value = 0.;
                mask |= 1 << j;
            }
        }
        let i = self.points.len();
        self.points.push(p);
        self.parameters.push(u);
        self.values.push(values);
        self.masks.push(mask);
        self.edges.insert((a, b, field), i);
        i
    }
    fn clip(&mut self, faces: Vec<Vec<usize>>, field: usize) -> Result<Vec<Vec<usize>>, String> {
        let mut out = Vec::new();
        let mut cap_edges = Vec::new();
        for face in faces {
            let mut poly = Vec::new();
            let mut cuts = Vec::new();
            for k in 0..face.len() {
                let a = face[k];
                let b = face[(k + 1) % face.len()];
                let ia = self.values[a][field] <= 0.;
                let ib = self.values[b][field] <= 0.;
                if ia {
                    poly.push(a);
                }
                if ia != ib {
                    let p = self.cut(a, b, field);
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
                cap_edges.push([cuts[0], cuts[1]]);
            }
        }
        if cap_edges.len() >= 3 {
            let first = cap_edges[0][0];
            let mut cap = vec![first];
            let mut current = first;
            let mut previous = usize::MAX;
            for _ in 0..cap_edges.len() {
                let next = cap_edges
                    .iter()
                    .find_map(|&[a, b]| {
                        if a == current && b != previous {
                            Some(b)
                        } else if b == current && a != previous {
                            Some(a)
                        } else {
                            None
                        }
                    })
                    .ok_or("intersection cap is not a cycle")?;
                if next == first {
                    break;
                }
                cap.push(next);
                previous = current;
                current = next;
            }
            if cap.len() != cap_edges.len() {
                return Err("intersection cap has disconnected pieces".into());
            }
            out.push(cap);
        }
        Ok(out)
    }
    fn triangles(&self, face: &[usize]) -> Vec<[usize; 3]> {
        let pivot = (0..face.len())
            .min_by(|&a, &b| self.compare_vertices(face[a], face[b]))
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
    // Rounded periodic ranks can coincide for distinct, very close cuts. Use
    // one total ordering on every face, including the cell's pulling anchor.
    fn compare_vertices(&self, a: usize, b: usize) -> std::cmp::Ordering {
        self.rank(a).cmp(&self.rank(b)).then_with(|| {
            for axis in 0..3 {
                let canonical = |i: usize| {
                    let p = self.parameters[i][axis];
                    if self.options.periodic[axis] && p == self.options.bounds[1][axis] {
                        self.options.bounds[0][axis]
                    } else {
                        p
                    }
                };
                let order = canonical(a).total_cmp(&canonical(b));
                if order != std::cmp::Ordering::Equal {
                    return order;
                }
            }
            a.cmp(&b)
        })
    }
    // Choose the best local filling without changing any shared-face diagonal.
    // A pulling anchor is admissible only when it agrees with the already
    // canonical face triangulations. An interior point is another conforming
    // candidate; it is committed only if it improves the worst element.
    fn tetrahedralize(&mut self, faces: &[Vec<usize>]) -> Result<Vec<Tet>, String> {
        let mut vertices: Vec<_> = faces.iter().flatten().copied().collect();
        vertices.sort_unstable();
        vertices.dedup();
        if vertices.is_empty() {
            return Ok(Vec::new());
        }
        vertices.sort_by(|&a, &b| self.compare_vertices(a, b));
        let triangles: Vec<_> = faces.iter().map(|f| self.triangles(f)).collect();
        let mut best = Vec::new();
        let mut score = -1.;
        for &anchor in &vertices {
            if faces
                .iter()
                .zip(&triangles)
                .any(|(f, tris)| f.contains(&anchor) && tris.iter().any(|t| !t.contains(&anchor)))
            {
                continue;
            }
            let candidate: Vec<_> = faces
                .iter()
                .zip(&triangles)
                .filter(|(f, _)| !f.contains(&anchor))
                .flat_map(|(_, tris)| tris.iter().map(|t| [anchor, t[0], t[1], t[2]]))
                .collect();
            let (candidate, quality) = orient_and_score(candidate, &self.points)?;
            if quality > score + 1e-12 {
                best = candidate;
                score = quality;
            }
        }
        if score < 0. {
            return Err("no conforming cell tetrahedralization".into());
        }
        if score < 0.1 && vertices.len() > 4 {
            let center: Point = std::array::from_fn(|a| {
                vertices.iter().map(|&i| self.points[i][a]).sum::<f64>() / vertices.len() as f64
            });
            let index = self.points.len();
            self.points.push(center);
            let candidate = triangles
                .iter()
                .flatten()
                .map(|t| [index, t[0], t[1], t[2]])
                .collect();
            let (candidate, quality) = orient_and_score(candidate, &self.points)?;
            if quality >= 0.1 && quality > score + 1e-12 {
                self.parameters.push(std::array::from_fn(|a| {
                    vertices.iter().map(|&i| self.parameters[i][a]).sum::<f64>()
                        / vertices.len() as f64
                }));
                self.values.push(
                    (0..self.values[vertices[0]].len())
                        .map(|a| {
                            vertices.iter().map(|&i| self.values[i][a]).sum::<f64>()
                                / vertices.len() as f64
                        })
                        .collect(),
                );
                self.masks.push(0);
                best = candidate;
            } else {
                self.points.pop();
            }
        }
        Ok(best)
    }
}

/// Intersect two to eight negative sublevel sets in one cell-construction pass.
///
/// `options` supplies the parameter box, grid, periodic ranking and element
/// budget, and optional joint snapping (0..=0.2). Region must be Negative.
/// Up to 20 optimization passes remove short edges, move constrained vertices
/// and reconnect interior tetrahedra. `minimum_quality` is an acceptance gate,
/// not a guarantee that the geometry admits that quality. Geometry tolerance
/// must be zero: adaptive physical error control is not implemented.
/// Callers must validate mapped periodicity and physical approximation error.
/// This API does not promise solver-quality tetrahedra.
pub fn generate(
    fields: &[&dyn ScalarField],
    map: impl Fn(Point) -> Point,
    options: Options,
) -> Result<Output, String> {
    generate_with_map(fields, &map, options)
}

fn generate_with_map(
    fields: &[&dyn ScalarField],
    map: &dyn Fn(Point) -> Point,
    options: Options,
) -> Result<Output, String> {
    if !(2..=8).contains(&fields.len())
        || (!options.snap.is_finite() || !(0.0..=0.2).contains(&options.snap))
        || options.optimize_passes > 20
        || options.geometry_tolerance != 0.
        || !options.minimum_quality.is_finite()
        || !(0.0..=1.).contains(&options.minimum_quality)
        || !matches!(options.region, super::Region::Negative)
    {
        return Err("intersection requires 2..8 fields, Negative region and snap 0..0.2 and passes 0..20, quality 0..1 and zero tolerance; geometric acceptance is caller-owned".into());
    }
    let [lo, hi] = options.bounds;
    if (0..3).any(|a| {
        !lo[a].is_finite()
            || !hi[a].is_finite()
            || hi[a] <= lo[a]
            || !(hi[a] - lo[a]).is_finite()
            || !(4..=128).contains(&options.cells[a])
    }) || options.max_tetrahedra == 0
    {
        return Err("invalid intersection bounds, grid or budget".into());
    }
    // Check each separate constraint, rather than only their maximum: an
    // inactive branch must not conceal a nonperiodic surface.
    for axis in 0..3 {
        if !options.periodic[axis] {
            continue;
        }
        for i in 0..11 {
            for j in 0..11 {
                let mut a = lo;
                a[(axis + 1) % 3] += (hi[(axis + 1) % 3] - lo[(axis + 1) % 3]) * i as f64 / 10.;
                a[(axis + 2) % 3] += (hi[(axis + 2) % 3] - lo[(axis + 2) % 3]) * j as f64 / 10.;
                let mut b = a;
                b[axis] = hi[axis];
                for field in fields {
                    let x = field.value(a);
                    let y = field.value(b);
                    if !x.is_finite()
                        || !y.is_finite()
                        || (x - y).abs() > 1e-10 * (1. + x.abs().max(y.abs()))
                    {
                        return Err(
                            "constraint is nonfinite or incompatible with periodic axis".into()
                        );
                    }
                }
            }
        }
    }
    // Reuse the validated BCC box topology. Fields are not meshed individually.
    let background = super::generate(
        &|_: Point| -1.,
        Options {
            geometry_tolerance: 1.,
            snap: 0.,
            optimize_passes: 0,
            minimum_quality: 0.,
            periodic: [false; 3],
            max_tetrahedra: usize::MAX,
            ..options
        },
    )
    .map_err(|e| e.to_string())?
    .mesh;
    let mut c = Cutter {
        options: &options,
        points: Vec::new(),
        parameters: background.points,
        values: Vec::new(),
        masks: Vec::new(),
        edges: HashMap::new(),
    };
    for &p in &c.parameters {
        let q = map(p);
        if q.iter().any(|x| !x.is_finite()) {
            return Err("nonfinite coordinate map".into());
        }
        c.points.push(q);
        let canonical = std::array::from_fn(|a| {
            if options.periodic[a] && p[a] == hi[a] {
                lo[a]
            } else {
                p[a]
            }
        });
        let values: Vec<_> = fields.iter().map(|f| f.value(canonical)).collect();
        if values.iter().any(|v| !v.is_finite()) {
            return Err("nonfinite constraint".into());
        }
        c.masks.push(
            values
                .iter()
                .enumerate()
                .fold(0, |mask, (j, &v)| mask | if v == 0. { 1 << j } else { 0 }),
        );
        c.values.push(values);
    }
    // Values at an analytically coincident background vertex can straddle zero
    // by roundoff (e.g. a Pythagorean point on a cylinder). Classify them once,
    // using each incident cell's field scale, before sharing any edge cuts.
    let mut coincident = vec![0_u8; c.values.len()];
    for t in &background.tets {
        for j in 0..fields.len() {
            let scale = t.iter().map(|&i| c.values[i][j].abs()).fold(0., f64::max);
            for &i in t {
                if c.values[i][j].abs() <= 32. * f64::EPSILON * scale {
                    coincident[i] |= 1 << j;
                }
            }
        }
    }
    for (i, mask) in coincident.into_iter().enumerate() {
        c.masks[i] |= mask;
        for j in 0..fields.len() {
            if mask & (1 << j) != 0 {
                c.values[i][j] = 0.;
            }
        }
    }
    let original_background = c.points.clone();
    let original_masks = c.masks.clone();
    if options.snap > 0. {
        snap_background(&mut c, fields, &map, &background.tets)?;
    }
    let snapped = c
        .points
        .iter()
        .zip(&original_background)
        .enumerate()
        .filter(|(i, (a, b))| a != b || c.masks[*i] != original_masks[*i])
        .count();
    drop(original_background);
    drop(original_masks);
    let mut tets: Vec<Tet> = Vec::new();
    for (index, t) in background.tets.into_iter().enumerate() {
        if index % 4096 == 0 {
            for field in fields {
                field.check()?;
            }
        }
        let background_det = determinant(t.map(|i| c.points[i]));
        if !background_det.is_finite() || background_det <= 0. {
            return Err("coordinate map inverted a background tetrahedron".into());
        }
        if (0..fields.len()).any(|j| t.iter().all(|&i| c.values[i][j] > 0.)) {
            continue;
        }
        let mut faces = vec![
            vec![t[0], t[1], t[2]],
            vec![t[0], t[3], t[1]],
            vec![t[0], t[2], t[3]],
            vec![t[1], t[3], t[2]],
        ];
        for j in 0..fields.len() {
            faces = c.clip(faces, j)?;
        }
        for tet in c.tetrahedralize(&faces)? {
            tets.push(tet);
            if tets.len() > options.max_tetrahedra {
                return Err("intersection tetrahedron budget exhausted".into());
            }
        }
    }
    let mut remap = vec![usize::MAX; c.points.len()];
    let mut points = Vec::new();
    let mut parameters = Vec::new();
    let mut constraints = Vec::new();
    for t in &mut tets {
        for i in t {
            if remap[*i] == usize::MAX {
                remap[*i] = points.len();
                points.push(c.points[*i]);
                parameters.push(c.parameters[*i]);
                constraints.push(c.masks[*i]);
            }
            *i = remap[*i];
        }
    }
    let boundary = crate::boundary(&tets).and_then(|surface| {
        if super::nonmanifold_vertices(&surface, points.len()).is_empty() {
            Ok(surface)
        } else {
            Err("intersection boundary is not manifold".into())
        }
    });
    let surface = match boundary {
        Ok(surface) => surface,
        Err(_) if snapped > 0 && options.snap > 0. => {
            // Snapping is an optional optimization, never grounds for accepting
            // a pinched solid. Retry the same fields, map, box and resolution.
            let mut result = generate_with_map(
                fields,
                map,
                Options {
                    snap: 0.,
                    ..options
                },
            )?;
            result.quality.reverted_snap_vertices = snapped;
            return Ok(result);
        }
        Err(error) => return Err(error),
    };
    let mut output = Output {
        mesh: Mesh {
            points,
            tets,
            surface,
        },
        parameters,
        constraints,
        quality: QualityDiagnostics::default(),
    };
    output.quality.initial_minimum_quality = 1.;
    for t in &output.mesh.tets {
        let q = quality(t.map(|i| output.mesh.points[i])).powf(1.5);
        output.quality.initial_minimum_quality = output.quality.initial_minimum_quality.min(q);
        output.quality.initial_elements_below_01 += usize::from(q < 0.1);
    }
    for _ in 0..options.optimize_passes {
        for field in fields {
            field.check()?;
        }
        let collapsed = repair::improve(
            &mut output,
            &Options {
                optimize_passes: 1,
                ..options
            },
        )?;
        let moved = repair::smooth(&mut output, &options, fields, &map)?;
        let flipped = repair::flip(&mut output)?;
        output.quality.collapsed_vertices += collapsed;
        output.quality.accepted_vertex_moves += moved;
        output.quality.accepted_reconnections += flipped;
        if collapsed == 0 && moved == 0 && flipped == 0 {
            break;
        }
    }
    if output.mesh.tets.len() > options.max_tetrahedra {
        return Err("quality improvement exhausted tetrahedron budget".into());
    }
    let mut minimum: f64 = 1.;
    for t in &output.mesh.tets {
        let p = t.map(|i| output.mesh.points[i]);
        let q = quality(p).powf(1.5);
        if !q.is_finite() || q <= 0. || determinant(p) <= 0. {
            return Err("invalid element after quality improvement".into());
        }
        minimum = minimum.min(q);
        output.quality.elements_below_01 += usize::from(q < 0.1);
    }
    output.quality.minimum_quality = minimum;
    if !output.mesh.tets.is_empty() && minimum < options.minimum_quality {
        return Err(format!(
            "minimum MMG quality {minimum:e} is below required {:e}; constraints were preserved",
            options.minimum_quality
        ));
    }
    Ok(output)
}

fn orient_and_score(tets: Vec<Tet>, points: &[Point]) -> Result<(Vec<Tet>, f64), String> {
    let mut result = Vec::new();
    let mut minimum: f64 = 1.;
    for mut t in tets {
        let p = t.map(|i| points[i]);
        let det = determinant(p);
        if !det.is_finite() {
            return Err("nonfinite intersection element".into());
        }
        if det == 0. {
            continue;
        }
        if det < 0. {
            t.swap(2, 3);
        }
        minimum = minimum.min(quality(p).powf(1.5));
        result.push(t);
    }
    Ok((result, minimum))
}

fn field_gradient(field: &dyn ScalarField, p: Point, step: f64) -> Point {
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

// Solve for the least-norm displacement satisfying the selected surface
// equations, with box-face coordinates locked. A singular junction is left
// unsnapped; no constraint is silently dropped to force a solution.
fn project_junction(
    fields: &[&dyn ScalarField],
    active: &[usize],
    original: Point,
    options: &Options,
    step: f64,
) -> Option<Point> {
    project_candidate(fields, active, original, original, options, step)
}

fn project_candidate(
    fields: &[&dyn ScalarField],
    active: &[usize],
    original: Point,
    candidate: Point,
    options: &Options,
    step: f64,
) -> Option<Point> {
    let mut p = candidate;
    let locked: [bool; 3] = std::array::from_fn(|a| {
        original[a] == options.bounds[0][a] || original[a] == options.bounds[1][a]
    });
    for a in 0..3 {
        if locked[a] {
            p[a] = original[a];
        }
    }
    for _ in 0..12 {
        let mut values: Vec<_> = active.iter().map(|&j| fields[j].value(p)).collect();
        if values.iter().all(|&v| v == 0.) {
            return Some(p);
        }
        let mut gradients: Vec<_> = active
            .iter()
            .map(|&j| field_gradient(fields[j], p, step))
            .collect();
        // A positive field multiplier must not change the surface or the solve.
        // Normalize before locking coordinates so residuals measure distance in
        // parameter space and Gram pivots measure angles, not field magnitudes.
        for (value, g) in values.iter_mut().zip(&mut gradients) {
            let magnitude = g[0].hypot(g[1]).hypot(g[2]);
            if !magnitude.is_finite() || magnitude == 0. || !value.is_finite() {
                return None;
            }
            *value /= magnitude;
            *g = g.map(|v| v / magnitude);
        }
        let extent = (0..3)
            .map(|a| options.bounds[1][a] - options.bounds[0][a])
            .fold(0., f64::max);
        if values.iter().all(|v| v.abs() < 1e-12 * extent) {
            return Some(p);
        }
        for g in &mut gradients {
            for a in 0..3 {
                if locked[a] {
                    g[a] = 0.;
                }
            }
        }
        let n = active.len();
        let mut matrix = vec![vec![0.; n + 1]; n];
        for i in 0..n {
            for j in 0..n {
                matrix[i][j] = dot(gradients[i], gradients[j]);
            }
            matrix[i][n] = values[i];
        }
        for k in 0..n {
            let pivot =
                (k..n).max_by(|&a, &b| matrix[a][k].abs().total_cmp(&matrix[b][k].abs()))?;
            matrix.swap(k, pivot);
            let d = matrix[k][k];
            if !d.is_finite() || d.abs() < 64. * f64::EPSILON {
                return None;
            }
            for value in &mut matrix[k][k..=n] {
                *value /= d;
            }
            let pivot_row = matrix[k].clone();
            for (i, row) in matrix.iter_mut().enumerate() {
                if i != k {
                    let factor = row[k];
                    for (value, pivot) in row[k..=n].iter_mut().zip(&pivot_row[k..=n]) {
                        *value -= factor * pivot;
                    }
                }
            }
        }
        for a in 0..3 {
            if !locked[a] {
                p[a] -= (0..n).map(|i| gradients[i][a] * matrix[i][n]).sum::<f64>();
            }
        }
        if p.iter().any(|x| !x.is_finite()) {
            return None;
        }
    }
    None
}

fn snap_background(
    c: &mut Cutter<'_>,
    fields: &[&dyn ScalarField],
    map: &impl Fn(Point) -> Point,
    tets: &[Tet],
) -> Result<(), String> {
    let o = c.options;
    let [lo, hi] = o.bounds;
    let h = (0..3)
        .map(|a| (hi[a] - lo[a]) / o.cells[a] as f64)
        .fold(f64::INFINITY, f64::min);
    let original = c.parameters.clone();
    let original_points = c.points.clone();
    let original_values = c.values.clone();
    let original_masks = c.masks.clone();
    let canonical = |p: Point| {
        std::array::from_fn(|a| {
            if o.periodic[a] && p[a] == hi[a] {
                lo[a]
            } else {
                p[a]
            }
        })
    };
    let rank = |p: Point| -> [i64; 3] {
        let p = canonical(p);
        std::array::from_fn(|a| ((p[a] - lo[a]) / (hi[a] - lo[a]) * 1e11).round() as i64)
    };
    let mut changed = vec![false; original.len()];
    for (i, &p) in original.iter().enumerate() {
        if i % 4096 == 0 {
            for f in fields {
                f.check()?;
            }
        }
        let p0 = canonical(p);
        let active: Vec<_> = fields
            .iter()
            .enumerate()
            .filter_map(|(j, f)| {
                let g = field_gradient(*f, p0, h * 1e-5);
                (c.values[i][j].abs() / norm(g) <= o.snap * h).then_some(j)
            })
            .collect();
        if active.is_empty() || active.len() > 3 {
            continue;
        }
        let Some(mut q) = project_junction(fields, &active, p0, o, h * 1e-5) else {
            continue;
        };
        if norm(sub(q, p0)) > o.snap * h || (0..3).any(|a| q[a] < lo[a] || q[a] > hi[a]) {
            continue;
        }
        let mut values: Vec<_> = fields.iter().map(|f| f.value(q)).collect();
        if values.iter().any(|v| !v.is_finite()) {
            continue;
        }
        let mut mask = 0;
        for j in active {
            values[j] = 0.;
            mask |= 1 << j;
        }
        for a in 0..3 {
            if o.periodic[a] && p[a] == hi[a] {
                q[a] = hi[a];
            }
        }
        let physical = map(q);
        if physical.iter().any(|x| !x.is_finite()) {
            continue;
        }
        changed[i] = q != p || values != c.values[i];
        c.parameters[i] = q;
        c.points[i] = physical;
        c.values[i] = values;
        c.masks[i] |= mask;
    }
    for _ in 0..16 {
        let mut undo = HashSet::new();
        for &t in tets {
            let p = t.map(|i| c.points[i]);
            if determinant(p) <= 0. || quality(p).powf(1.5) < 0.02 {
                for i in t {
                    if changed[i] {
                        undo.insert(rank(original[i]));
                    }
                }
            }
        }
        if undo.is_empty() {
            return Ok(());
        }
        for (i, &p) in original.iter().enumerate() {
            if changed[i] && undo.contains(&rank(p)) {
                c.parameters[i] = p;
                c.points[i] = original_points[i];
                c.values[i] = original_values[i].clone();
                c.masks[i] = original_masks[i];
                changed[i] = false;
            }
        }
    }
    Err("joint snapping did not stabilize background elements".into())
}
