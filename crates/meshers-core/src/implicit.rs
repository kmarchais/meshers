//! CPU meshing of bounded scalar fields. Geometry is evaluated in physical coordinates.
//! Sampling is not a proof that arbitrarily small components have been found.
pub mod intersection;

use crate::{Mesh, MeshingError, Point, determinant, dot, geometry::Geometry, norm, quality, sub};
use std::collections::{BTreeMap, HashMap};
use std::sync::atomic::{AtomicBool, Ordering};

/// A deterministic spatial function. Values and optional gradients must be finite.
/// The function must be defined in a small neighborhood of the box for numerical derivatives.
pub trait ScalarField: Sync {
    fn value(&self, point: Point) -> f64;
    /// Prefer bounded vector evaluation in optimization. Scalar implementations keep their path.
    fn prefers_batches(&self) -> bool {
        false
    }
    fn gradients(&self, points: &[Point]) -> Vec<Option<Point>> {
        points.iter().map(|&p| self.gradient(p)).collect()
    }
    fn values(&self, points: &[Point]) -> Vec<f64> {
        points.iter().map(|&p| self.value(p)).collect()
    }
    /// Cooperative interruption hook for language bindings or application cancellation.
    fn check(&self) -> Result<(), String> {
        Ok(())
    }
    fn gradient(&self, _point: Point) -> Option<Point> {
        None
    }
}
impl<F: Fn(Point) -> f64 + Sync> ScalarField for F {
    fn value(&self, point: Point) -> f64 {
        self(point)
    }
}

#[derive(Clone, Copy, Debug)]
pub enum Region {
    /// The negative sublevel set of the supplied function.
    Negative,
    /// A constant interval in raw field units. Width is upper minus lower.
    Band { lower: f64, upper: f64 },
}

#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Minimum and maximum corners in physical coordinates.
    pub bounds: [Point; 2],
    pub cells: [usize; 3],
    pub region: Region,
    /// Explicit translational constraints, one per box axis. Default: none.
    pub periodic: [bool; 3],
    /// Maximum sampled implicit-boundary distance estimate in physical units.
    pub geometry_tolerance: f64,
    /// Required minimum MMG element quality in `[0,1]`. Zero disables the quality gate.
    pub minimum_quality: f64,
    pub max_tetrahedra: usize,
    pub optimize_passes: usize,
    /// Maximum constrained lattice snap displacement as a fraction of the smallest cell spacing.
    pub snap: f64,
    /// Worker count. Defaults to one; 0 retains legacy serial ordering.
    pub threads: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            bounds: [[0.; 3], [1.; 3]],
            cells: [24; 3],
            region: Region::Negative,
            periodic: [false; 3],
            geometry_tolerance: 0.01,
            minimum_quality: 0.,
            max_tetrahedra: 2_000_000,
            optimize_passes: 4,
            snap: 0.2,
            threads: 1,
        }
    }
}

#[derive(Debug)]
pub struct Diagnostics {
    pub empty: bool,
    pub volume: f64,
    pub minimum_mmg_quality: f64,
    pub elements_below_quality_01: usize,
    pub reverted_snap_vertices: usize,
    pub maximum_sampled_surface_error: f64,
}
pub struct Output {
    pub mesh: Mesh,
    /// One tag per boundary triangle: 0 implicit, 1/2 x-/x+, 3/4 y-/y+, 5/6 z-/z+.
    pub boundary_tags: Vec<u8>,
    /// Opposite node indices per requested axis. Translations are box side lengths.
    pub periodic_pairs: [Vec<[usize; 2]>; 3],
    pub diagnostics: Diagnostics,
}

struct Context<'a, F: ScalarField> {
    field: &'a F,
    options: Options,
    center: f64,
    scale: f64,
    invalid: AtomicBool,
}
impl<F: ScalarField> Context<'_, F> {
    fn canonical(&self, mut p: Point) -> Point {
        for (a, x) in p.iter_mut().enumerate() {
            if self.options.periodic[a] && *x == self.options.bounds[1][a] {
                *x = self.options.bounds[0][a];
            }
        }
        p
    }
    fn checked(&self, value: f64) -> f64 {
        if value.is_finite() {
            value
        } else {
            self.invalid.store(true, Ordering::Relaxed);
            0.
        }
    }
    fn threshold(&self) -> f64 {
        match self.options.region {
            Region::Negative => 0.,
            Region::Band { lower, upper } => (upper - lower) / 2. / self.scale,
        }
    }
    fn distance(&self, p: Point) -> f64 {
        let residual = self.residual(p, self.threshold()).abs();
        if residual == 0. {
            return 0.;
        }
        let gradient = norm(self.gradient(p));
        if gradient > 0. {
            residual / gradient
        } else {
            f64::INFINITY
        }
    }
}
impl<F: ScalarField> Geometry for Context<'_, F> {
    fn batched(&self) -> bool {
        self.field.prefers_batches()
    }
    fn values(&self, points: &[Point]) -> Vec<f64> {
        let canonical: Vec<_> = points.iter().map(|&p| self.canonical(p)).collect();
        let values = self.field.values(&canonical);
        if values.len() != points.len() {
            self.invalid.store(true, Ordering::Relaxed);
            return vec![0.; points.len()];
        }
        values
            .into_iter()
            .map(|v| self.checked((v - self.center) / self.scale))
            .collect()
    }
    fn gradients(&self, points: &[Point]) -> Vec<Point> {
        let canonical: Vec<_> = points.iter().map(|&p| self.canonical(p)).collect();
        let gradients = self.field.gradients(&canonical);
        if gradients.len() != points.len() {
            self.invalid.store(true, Ordering::Relaxed);
            return vec![[0.; 3]; points.len()];
        }
        gradients
            .into_iter()
            .zip(points)
            .map(|(g, &p)| {
                g.map_or_else(
                    || self.gradient(p),
                    |g| g.map(|v| self.checked(v / self.scale)),
                )
            })
            .collect()
    }
    fn value(&self, p: Point) -> f64 {
        self.checked((self.field.value(self.canonical(p)) - self.center) / self.scale)
    }
    fn gradient(&self, p: Point) -> Point {
        let p = self.canonical(p);
        if let Some(g) = self.field.gradient(p) {
            return g.map(|v| self.checked(v / self.scale));
        }
        std::array::from_fn(|a| {
            let h = (self.options.bounds[1][a] - self.options.bounds[0][a]) * 1e-5;
            let mut lo = p;
            let mut hi = p;
            lo[a] -= h;
            hi[a] += h;
            self.checked((self.field.value(hi) - self.field.value(lo)) / (2. * h * self.scale))
        })
    }
    fn rank(&self, p: Point) -> [i64; 3] {
        let p = self.canonical(p);
        std::array::from_fn(|a| {
            ((p[a] - self.options.bounds[0][a])
                / (self.options.bounds[1][a] - self.options.bounds[0][a])
                * 1e11)
                .round() as i64
        })
    }
    fn bounds(&self) -> [Point; 2] {
        self.options.bounds
    }
    fn check(&self) -> Result<(), String> {
        self.field.check()
    }
    fn band(&self) -> bool {
        matches!(self.options.region, Region::Band { .. })
    }
}

/// Generate a closed solid boundary and conforming tetrahedra directly from a field.
/// Custom functions and grading use the same evaluator in roots, optimization and refinement.
/// Unmet sampled tolerance or invalid topology returns an error, never a partial success.
pub fn generate<F: ScalarField>(field: &F, options: Options) -> Result<Output, MeshingError> {
    let [lo, hi] = options.bounds;
    if !(0..3).all(|a| {
        lo[a].is_finite()
            && hi[a].is_finite()
            && hi[a] > lo[a]
            && (hi[a] - lo[a]).is_finite()
            && (hi[a] - lo[a]) * 1e-5 > f64::EPSILON * lo[a].abs().max(hi[a].abs())
            && (4..=128).contains(&options.cells[a])
    }) || !options.geometry_tolerance.is_finite()
        || options.geometry_tolerance <= 0.
        || !options.snap.is_finite()
        || !(0.0..=0.2).contains(&options.snap)
        || options.threads > 256
        || options.optimize_passes > 100
        || !options.minimum_quality.is_finite()
        || !(0.0..=1.).contains(&options.minimum_quality)
        || options.max_tetrahedra == 0
    {
        return Err(MeshingError::InvalidOptions("finite, representable positive bounds; cells 4..128; positive tolerance/budget; snap 0..0.2; threads <=256; passes <=100; minimum_quality in 0..1 required".into()));
    }
    let center = match options.region {
        Region::Negative => 0.,
        Region::Band { lower, upper } => {
            if !lower.is_finite()
                || !upper.is_finite()
                || lower >= upper
                || !(upper - lower).is_finite()
            {
                return Err(MeshingError::InvalidOptions(
                    "band endpoints must be finite and strictly increasing".into(),
                ));
            }
            lower / 2. + upper / 2.
        }
    };
    let mut scale: f64 = match options.region {
        Region::Negative => 0.,
        Region::Band { lower, upper } => (upper - lower) / 2.,
    };
    for x in [0., 0.23, 0.61, 1.] {
        for y in [0., 0.23, 0.61, 1.] {
            for z in [0., 0.23, 0.61, 1.] {
                let p = std::array::from_fn(|a| lo[a] + [x, y, z][a] * (hi[a] - lo[a]));
                let v = field.value(p) - center;
                if !v.is_finite() {
                    return Err(MeshingError::GenerationFailed(
                        "nonfinite field value".into(),
                    ));
                }
                scale = scale.max(v.abs());
            }
        }
    }
    if scale == 0. {
        return Err(MeshingError::GenerationFailed(
            "field is zero on all preflight samples; no resolved solid boundary".into(),
        ));
    }
    let context = Context {
        field,
        options,
        center,
        scale,
        invalid: AtomicBool::new(false),
    };
    // Inspect raw field traces before canonicalizing opposite faces.
    for axis in 0..3 {
        if !options.periodic[axis] {
            continue;
        }
        for i in 0..11 {
            for j in 0..11 {
                let mut p = lo;
                p[(axis + 1) % 3] += (hi[(axis + 1) % 3] - lo[(axis + 1) % 3]) * i as f64 / 10.;
                p[(axis + 2) % 3] += (hi[(axis + 2) % 3] - lo[(axis + 2) % 3]) * j as f64 / 10.;
                let mut q = p;
                q[axis] = hi[axis];
                let difference = (field.value(p) - field.value(q)) / scale;
                if !difference.is_finite() || difference.abs() > 1e-8 {
                    return Err(MeshingError::InvalidOptions(format!(
                        "field/grading is incompatible with periodic axis {axis}"
                    )));
                }
            }
        }
    }
    context.check().map_err(MeshingError::GenerationFailed)?;
    let result = build(&context);
    if context.invalid.load(Ordering::Relaxed) {
        return Err(MeshingError::GenerationFailed(
            "field or gradient returned nonfinite values".into(),
        ));
    }
    result.map_err(MeshingError::GenerationFailed)
}

// Every closed surface edge has two incident triangles, and every vertex link
// must be one cycle. Edge counts alone miss two shells touching at one vertex.
fn nonmanifold_vertices(surface: &[[usize; 3]], point_count: usize) -> Vec<usize> {
    let mut links = vec![Vec::<[usize; 2]>::new(); point_count];
    let mut edges = HashMap::<[usize; 2], usize>::new();
    for &[a, b, c] in surface {
        links[a].push([b, c]);
        links[b].push([c, a]);
        links[c].push([a, b]);
        for [i, j] in [[a, b], [b, c], [c, a]] {
            *edges.entry([i.min(j), i.max(j)]).or_default() += 1;
        }
    }
    let mut bad = vec![false; point_count];
    for (edge, count) in edges {
        if count != 2 {
            for v in edge {
                bad[v] = true;
            }
        }
    }
    for (v, link) in links.iter().enumerate() {
        if bad[v] || link.is_empty() {
            continue;
        }
        let mut reached = vec![link[0][0]];
        let mut cursor = 0;
        while cursor < reached.len() {
            let current = reached[cursor];
            cursor += 1;
            for &[a, b] in link {
                let next = if a == current {
                    b
                } else if b == current {
                    a
                } else {
                    continue;
                };
                if !reached.contains(&next) {
                    reached.push(next);
                }
            }
        }
        // Degree two is already established by the edge counts above.
        bad[v] = reached.len() != link.len();
    }
    bad.into_iter()
        .enumerate()
        .filter_map(|(i, bad)| bad.then_some(i))
        .collect()
}

fn finish<F: ScalarField>(mesh: Mesh, g: &Context<F>) -> Result<Output, String> {
    let mut volume = 0.;
    let mut minimum: f64 = 1.;
    let mut elements_below_quality_01 = 0;
    for t in &mesh.tets {
        let p = t.map(|i| mesh.points[i]);
        let det = determinant(p);
        if !det.is_finite() || det <= 0. {
            return Err("invalid or inverted tetrahedron".into());
        }
        volume += det / 6.;
        let q = quality(p).powf(1.5);
        if !q.is_finite() || q <= 0. {
            return Err("nonfinite or nonpositive tetrahedron quality".into());
        }
        minimum = minimum.min(q);
        elements_below_quality_01 += usize::from(q < 0.1);
    }
    let mut tags = Vec::new();
    let mut error: f64 = 0.;
    let mut samples = Vec::new();
    let [lo, hi] = g.bounds();
    for f in &mesh.surface {
        let p = f.map(|i| mesh.points[i]);
        let mut tag = 0;
        for a in 0..3 {
            if p.iter().all(|v| v[a] == lo[a]) {
                tag = 1 + 2 * a as u8;
            }
            if p.iter().all(|v| v[a] == hi[a]) {
                tag = 2 + 2 * a as u8;
            }
        }
        tags.push(tag);
        if tag == 0 {
            for sample in [
                p[0],
                p[1],
                p[2],
                std::array::from_fn(|a| (p[0][a] + p[1][a] + p[2][a]) / 3.),
                std::array::from_fn(|a| (p[0][a] + p[1][a]) / 2.),
                std::array::from_fn(|a| (p[1][a] + p[2][a]) / 2.),
                std::array::from_fn(|a| (p[2][a] + p[0][a]) / 2.),
            ] {
                if g.batched() {
                    samples.push(sample);
                    if samples.len() == 4096 {
                        error = error.max(distance_batch(g, &samples)?);
                        samples.clear();
                    }
                } else {
                    error = error.max(g.distance(sample));
                }
            }
        }
    }
    if !samples.is_empty() {
        error = error.max(distance_batch(g, &samples)?);
    }
    if !nonmanifold_vertices(&mesh.surface, mesh.points.len()).is_empty() {
        return Err("solid boundary is not a closed two-manifold at its edges or vertices".into());
    }
    if !error.is_finite() || error > g.options.geometry_tolerance {
        return Err(format!(
            "sampled geometry error {error:e} exceeds tolerance {:e}; increase resolution or budget",
            g.options.geometry_tolerance
        ));
    }
    let mut pairs: [Vec<[usize; 2]>; 3] = std::array::from_fn(|_| Vec::new());
    for (a, pair) in pairs.iter_mut().enumerate() {
        if !g.options.periodic[a] {
            continue;
        }
        let mut maps = [BTreeMap::new(), BTreeMap::new()];
        // Preserve other coordinates here: pairing on x must not collapse y/z corners.
        let key = |p: Point| {
            std::array::from_fn::<_, 3, _>(|axis| {
                if axis == a {
                    0
                } else {
                    ((p[axis] - lo[axis]) / (hi[axis] - lo[axis]) * 1e11).round() as i64
                }
            })
        };
        for (i, &p) in mesh.points.iter().enumerate() {
            for (side, value) in [lo[a], hi[a]].iter().enumerate() {
                if p[a] == *value && maps[side].insert(key(p), i).is_some() {
                    return Err("ambiguous periodic node match".into());
                }
            }
        }
        if maps[0].keys().ne(maps[1].keys()) {
            return Err(format!("periodic node mismatch on axis {a}"));
        }
        for (k, &i) in &maps[0] {
            pair.push([i, maps[1][k]]);
        }
        let mut faces = [Vec::new(), Vec::new()];
        for (f, &tag) in mesh.surface.iter().zip(&tags) {
            for (side, list) in faces.iter_mut().enumerate() {
                if tag == 1 + 2 * a as u8 + side as u8 {
                    let mut keys = f.map(|i| key(mesh.points[i]));
                    keys.sort_unstable();
                    list.push(keys);
                }
            }
        }
        faces.iter_mut().for_each(|f| f.sort_unstable());
        if faces[0] != faces[1] {
            return Err(format!("periodic triangle mismatch on axis {a}"));
        }
    }
    let empty = mesh.tets.is_empty();
    if !empty && minimum < g.options.minimum_quality {
        return Err(format!(
            "minimum MMG quality {minimum:e} is below required {:e}; mesh rejected ({} elements below 0.1). Improve resolution/optimization or use a quality remesher",
            g.options.minimum_quality, elements_below_quality_01
        ));
    }
    Ok(Output {
        mesh,
        boundary_tags: tags,
        periodic_pairs: pairs,
        diagnostics: Diagnostics {
            empty,
            volume,
            minimum_mmg_quality: if empty { 0. } else { minimum },
            maximum_sampled_surface_error: error,
            elements_below_quality_01,
            reverted_snap_vertices: 0,
        },
    })
}

fn build<F: ScalarField>(g: &Context<F>) -> Result<Output, String> {
    let o = g.options;
    let [nx, ny, nz] = o.cells;
    let [lo, hi] = o.bounds;
    let step: Point = std::array::from_fn(|a| (hi[a] - lo[a]) / o.cells[a] as f64);
    let h = step.into_iter().fold(f64::INFINITY, f64::min);
    let threshold = g.threshold();
    let mut b = crate::Builder {
        geometry: g,
        points: Vec::new(),
        values: Vec::new(),
        edges: HashMap::new(),
        threshold,
    };
    for z in 0..=nz {
        for y in 0..=ny {
            for x in 0..=nx {
                let ijk = [x, y, z];
                let p = std::array::from_fn(|a| {
                    if ijk[a] == o.cells[a] {
                        hi[a]
                    } else {
                        lo[a] + ijk[a] as f64 * step[a]
                    }
                });
                b.add(p, 0.);
            }
        }
    }
    for z in 0..nz {
        g.check()?;
        for y in 0..ny {
            for x in 0..nx {
                let p = std::array::from_fn(|a| lo[a] + ([x, y, z][a] as f64 + 0.5) * step[a]);
                b.add(p, 0.);
            }
        }
    }
    for (points, values) in b.points.chunks(4096).zip(b.values.chunks_mut(4096)) {
        g.check()?;
        let canonical: Vec<_> = points.iter().map(|&p| g.canonical(p)).collect();
        let sampled = g.field.values(&canonical);
        if sampled.len() != values.len() {
            return Err("field batch has wrong length".into());
        }
        for (target, value) in values.iter_mut().zip(sampled) {
            *target = g.checked((value - g.center) / g.scale);
        }
        if g.invalid.load(Ordering::Relaxed) {
            return Err("nonfinite field batch".into());
        }
    }
    let original_points = b.points.clone();
    let original_values = b.values.clone();
    if o.snap > 0. && g.batched() {
        snap_batch(&mut b.points, &mut b.values, g, threshold, o.snap * h)?;
    } else if o.snap > 0. {
        for i in 0..b.points.len() {
            if i % 256 == 0 {
                g.check()?;
            }
            let original = b.points[i];
            let value = b.values[i];
            let level = g.level(value, threshold);
            let mut grad = g.gradient(original);
            for (a, component) in grad.iter_mut().enumerate() {
                if g.locked(original)[a] {
                    *component = 0.;
                }
            }
            if norm(grad) < 1e-12 || (value - level).abs() / norm(grad) > o.snap * h {
                continue;
            }
            let mut p = original;
            for _ in 0..8 {
                let f = g.value(p) - level;
                let mut grad = g.gradient(p);
                for (a, component) in grad.iter_mut().enumerate() {
                    if g.locked(original)[a] {
                        *component = 0.;
                    }
                }
                let g2 = dot(grad, grad);
                if g2 < 1e-20 {
                    break;
                }
                p = std::array::from_fn(|a| p[a] - f * grad[a] / g2);
            }
            if g.contains(p)
                && norm(sub(p, original)) <= o.snap * h
                && (g.value(p) - level).abs() < 1e-11
            {
                b.points[i] = p;
                b.values[i] = level;
            }
        }
    }
    let mut mesh = None;
    let mut reverted_snap_vertices = 0;
    for _ in 0..8 {
        b.points.truncate(original_points.len());
        b.values.truncate(original_values.len());
        b.edges.clear();
        let (candidate, sources) = clip_lattice(&mut b, g)?;
        let bad = nonmanifold_vertices(&candidate.surface, candidate.points.len());
        if bad.is_empty() {
            mesh = Some(candidate);
            break;
        }
        let mut undo = std::collections::BTreeSet::new();
        for v in bad {
            let i = sources[v];
            if i < original_points.len()
                && (b.points[i] != original_points[i] || b.values[i] != original_values[i])
            {
                undo.insert(g.rank(original_points[i]));
            }
        }
        if undo.is_empty() {
            return Err(
                "solid boundary is not a closed two-manifold at its edges or vertices".into(),
            );
        }
        // Roll back all periodic copies together, before recomputing shared cuts.
        for i in 0..original_points.len() {
            if undo.contains(&g.rank(original_points[i]))
                && (b.points[i] != original_points[i] || b.values[i] != original_values[i])
            {
                reverted_snap_vertices += 1;
                b.points[i] = original_points[i];
                b.values[i] = original_values[i];
            }
        }
    }
    let mut mesh =
        mesh.ok_or("safe snapping could not resolve boundary topology in eight rounds")?;
    drop(b);
    drop(original_points);
    drop(original_values);
    if mesh.tets.is_empty() {
        return finish(mesh, g);
    }
    crate::optimize::optimize_for(
        &mut mesh,
        threshold,
        h,
        o.optimize_passes,
        3.,
        o.threads,
        None,
        g,
    )?;
    crate::adaptive::refine_for(
        &mut mesh,
        threshold,
        o.geometry_tolerance,
        o.max_tetrahedra,
        g,
    )?;
    crate::optimize::optimize_for(
        &mut mesh,
        threshold,
        h,
        o.optimize_passes,
        3.,
        o.threads,
        None,
        g,
    )?;
    let mut output = finish(mesh, g)?;
    output.diagnostics.reverted_snap_vertices = reverted_snap_vertices;
    Ok(output)
}

fn clip_lattice<F: ScalarField>(
    b: &mut crate::Builder<'_, Context<'_, F>>,
    g: &Context<F>,
) -> Result<(Mesh, Vec<usize>), String> {
    let o = g.options;
    let [nx, ny, nz] = o.cells;
    let [lo, hi] = o.bounds;
    let h = (0..3)
        .map(|a| (hi[a] - lo[a]) / o.cells[a] as f64)
        .fold(f64::INFINITY, f64::min);
    let threshold = g.threshold();
    let center_start = (nx + 1) * (ny + 1) * (nz + 1);
    let id = |x: usize, y: usize, z: usize| x + (nx + 1) * (y + (ny + 1) * z);
    let cube_faces = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [3, 7, 6, 2],
        [0, 4, 7, 3],
        [1, 2, 6, 5],
    ];
    let mut tets = Vec::new();
    for z in 0..nz {
        g.check()?;
        for y in 0..ny {
            for x in 0..nx {
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
                let c = center_start + x + nx * (y + ny * z);
                // True BCC interior: adjacent cell centers plus one shared
                // face edge form a disphenoid. Boundary pyramids end at the cube.
                for (fi, face) in cube_faces.iter().enumerate() {
                    let face = face.map(|i| corners[i]);
                    let coord = [z, z, y, y, x, x][fi];
                    let positive = fi % 2 == 1;
                    let boundary_face = if positive {
                        coord == [nz, nz, ny, ny, nx, nx][fi] - 1
                    } else {
                        coord == 0
                    };
                    let mut background = Vec::with_capacity(4);
                    if boundary_face {
                        for tri in b.triangulate(&face) {
                            background.push([c, tri[0], tri[1], tri[2]]);
                        }
                    } else if positive {
                        let neighbor = c + [nx * ny, nx * ny, nx, nx, 1, 1][fi];
                        for j in 0..4 {
                            background.push([c, neighbor, face[j], face[(j + 1) % 4]]);
                        }
                    } else {
                        continue;
                    }
                    for tet in background {
                        if tet.iter().all(|&i| b.values[i] > threshold + 1e-12)
                            || (g.band() && tet.iter().all(|&i| b.values[i] < -threshold - 1e-12))
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
                        let faces = if g.band() { b.clip(faces, -1) } else { faces };
                        let Some(anchor) = faces
                            .iter()
                            .flatten()
                            .copied()
                            .min_by_key(|&i| g.rank(b.points[i]))
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
                                if det.abs() < h * h * h * 1e-14 {
                                    continue;
                                }
                                if det < 0. {
                                    t.swap(2, 3);
                                }
                                tets.push(t);
                                if tets.len() > o.max_tetrahedra {
                                    return Err("tetrahedron budget exhausted".into());
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    let mut remap = vec![usize::MAX; b.points.len()];
    let mut points = Vec::new();
    let mut sources = Vec::new();
    for t in &mut tets {
        for i in t {
            if remap[*i] == usize::MAX {
                remap[*i] = points.len();
                points.push(b.points[*i]);
                sources.push(*i);
            }
            *i = remap[*i];
        }
    }
    let surface = crate::boundary(&tets)?;
    let mesh = Mesh {
        points,
        tets,
        surface,
    };
    Ok((mesh, sources))
}

fn snap_batch<F: ScalarField>(
    points: &mut [Point],
    values: &mut [f64],
    g: &Context<F>,
    threshold: f64,
    limit: f64,
) -> Result<(), String> {
    for (points, values) in points.chunks_mut(4096).zip(values.chunks_mut(4096)) {
        g.check()?;
        let originals = points.to_vec();
        let levels: Vec<_> = values.iter().map(|&v| g.level(v, threshold)).collect();
        let gradients = g.gradients(points);
        let mut selected = Vec::new();
        for (i, mut grad) in gradients.into_iter().enumerate() {
            for (a, v) in grad.iter_mut().enumerate() {
                if g.locked(originals[i])[a] {
                    *v = 0.;
                }
            }
            if norm(grad) >= 1e-12 && (values[i] - levels[i]).abs() / norm(grad) <= limit {
                selected.push(i);
            }
        }
        let mut active = selected.clone();
        for _ in 0..8 {
            if active.is_empty() {
                break;
            }
            g.check()?;
            let samples: Vec<_> = active.iter().map(|&i| points[i]).collect();
            let field = g.values(&samples);
            let gradients = g.gradients(&samples);
            let mut next = Vec::new();
            for ((i, value), mut grad) in active.into_iter().zip(field).zip(gradients) {
                for (a, v) in grad.iter_mut().enumerate() {
                    if g.locked(originals[i])[a] {
                        *v = 0.;
                    }
                }
                let g2 = dot(grad, grad);
                if g2 < 1e-20 {
                    continue;
                }
                points[i] =
                    std::array::from_fn(|a| points[i][a] - (value - levels[i]) * grad[a] / g2);
                next.push(i);
            }
            active = next;
        }
        let samples: Vec<_> = selected.iter().map(|&i| points[i]).collect();
        for (i, v) in selected.into_iter().zip(g.values(&samples)) {
            if g.contains(points[i])
                && norm(sub(points[i], originals[i])) <= limit
                && (v - levels[i]).abs() < 1e-11
            {
                values[i] = levels[i];
            } else {
                points[i] = originals[i];
            }
        }
    }
    g.check()
}

fn distance_batch<F: ScalarField>(g: &Context<F>, points: &[Point]) -> Result<f64, String> {
    g.check()?;
    let values = g.values(points);
    let gradients = g.gradients(points);
    let mut error = 0_f64;
    for (v, grad) in values.into_iter().zip(gradients) {
        let residual = if g.band() { v.abs() - g.threshold() } else { v };
        let n = norm(grad);
        let distance = if residual == 0. {
            0.
        } else if n > 0. {
            residual.abs() / n
        } else {
            f64::INFINITY
        };
        error = error.max(distance);
    }
    g.check()?;
    Ok(error)
}

#[cfg(test)]
mod topology_tests {
    use super::nonmanifold_vertices;

    #[test]
    fn catches_vertex_pinches_even_when_all_edges_have_two_faces() {
        let surface = crate::boundary(&[[0, 1, 2, 3], [0, 4, 5, 6]]).unwrap();
        assert_eq!(nonmanifold_vertices(&surface, 7), vec![0]);
        let separated = crate::boundary(&[[0, 1, 2, 3], [4, 5, 6, 7]]).unwrap();
        assert!(nonmanifold_vertices(&separated, 8).is_empty());
    }
}
