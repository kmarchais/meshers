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
    /// Optional length-ranked interior bisection after optimization. Zero disables it.
    pub interior_refine_nodes: usize,
    /// Maximum constrained lattice snap displacement as a fraction of the smallest cell spacing.
    pub snap: f64,
    /// Worker count. Defaults to one; 0 retains legacy serial ordering.
    pub threads: usize,
    /// Size-driven adaptation with exact boundary projection instead of fixed
    /// surface refinement. `None` keeps the lattice pipeline.
    pub adapt: Option<crate::tet_adapt::Options>,
    /// Element order: 1 (linear) or 2 (10-node, wall edges curved onto the field).
    pub order: u8,
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
            interior_refine_nodes: 0,
            snap: 0.2,
            threads: 1,
            adapt: None,
            order: 1,
        }
    }
}

#[derive(Debug)]
pub struct Diagnostics {
    pub retained_constraints: usize,
    pub polyhedral_constraints: usize,
    pub polyhedral_regions: usize,
    pub empty: bool,
    pub volume: f64,
    pub minimum_mmg_quality: f64,
    pub elements_below_quality_01: usize,
    pub reverted_snap_vertices: usize,
    pub inserted_interior_nodes: usize,
    pub maximum_sampled_surface_error: f64,
    pub adapt: Option<crate::tet_adapt::Stats>,
}
pub struct QuadraticOutput {
    pub elements: crate::quadratic::Quadratic,
    pub periodic_pairs: [Vec<[usize; 2]>; 3],
}
pub struct Output {
    /// Linear mesh; with `order: 2` its points are the vertices of `quadratic`.
    pub quadratic: Option<QuadraticOutput>,
    pub mesh: Mesh,
    /// One tag per boundary triangle: 0 implicit, 1/2 x-/x+, 3/4 y-/y+, 5/6 z-/z+.
    pub boundary_tags: Vec<u8>,
    /// Opposite node indices per requested axis. Translations are box side lengths.
    pub periodic_pairs: [Vec<[usize; 2]>; 3],
    pub diagnostics: Diagnostics,
}

struct Context<'a, F: ScalarField> {
    field: &'a F,
    constraints: Option<&'a [&'a dyn ScalarField]>,
    halfspaces: Option<&'a [[f64; 4]]>,
    halfspace_regions: Option<&'a [intersection::ConstraintMask]>,
    size_field: Option<&'a dyn ScalarField>,
    options: Options,
    center: f64,
    scale: f64,
    optimizer_length: f64,
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
    /// First-order distance to the wall, moving only along unlocked axes.
    fn distance_in(&self, p: Point, locked: [bool; 3]) -> f64 {
        let residual = self.residual(p, self.threshold()).abs();
        if residual == 0. {
            return 0.;
        }
        let mut gradient = self.gradient(p);
        for a in 0..3 {
            if locked[a] {
                gradient[a] = 0.;
            }
        }
        let gradient = norm(gradient);
        if gradient > 0. {
            residual / gradient
        } else {
            f64::INFINITY
        }
    }
}
impl<F: ScalarField> Geometry for Context<'_, F> {
    fn constraints(&self) -> Option<&[&dyn ScalarField]> {
        self.constraints
    }
    fn halfspaces(&self) -> Option<&[[f64; 4]]> {
        self.halfspaces
    }
    fn halfspace_regions(&self) -> Option<&[intersection::ConstraintMask]> {
        self.halfspace_regions
    }
    fn optimizer_length(&self) -> f64 {
        self.optimizer_length
    }
    fn minimum_quality(&self) -> f64 {
        self.options.minimum_quality
    }
    fn refinement_quality_floor(&self) -> f64 {
        // Experimental: allow curvature refinement down to the requested
        // FEA shape gate, retaining the historical ceiling for higher targets.
        if self.options.minimum_quality > 0. {
            self.options.minimum_quality.clamp(0.1, 0.25)
        } else {
            0.25
        }
    }
    fn geometry_tolerance(&self) -> Option<f64> {
        Some(self.options.geometry_tolerance)
    }
    fn batched(&self) -> bool {
        self.field.prefers_batches()
    }
    fn target_lengths(&self, points: &[Point]) -> Option<Vec<f64>> {
        self.size_field.map(|field| {
            let canonical: Vec<_> = points.iter().map(|&point| self.canonical(point)).collect();
            field.values(&canonical)
        })
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
        // Equal periodic traces can meet with different physical normal derivatives.
        // Keep the closed-cell side when evaluating normals, including exact caps.
        let gradients = self.field.gradients(points);
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
        self.field.check()?;
        if let Some(fields) = self.constraints {
            for field in fields {
                field.check()?;
            }
        }
        if let Some(field) = self.size_field {
            field.check()?;
        }
        Ok(())
    }
    fn band(&self) -> bool {
        matches!(self.options.region, Region::Band { .. })
    }
}

/// Generate a closed solid boundary and conforming tetrahedra directly from a field.
/// Custom functions and grading use the same evaluator in roots, optimization and refinement.
/// Unmet sampled tolerance or invalid topology returns an error, never a partial success.
pub fn generate<F: ScalarField>(field: &F, options: Options) -> Result<Output, MeshingError> {
    generate_with_size(field, None, options)
}

/// Generate with an optional positive spatial target length for tet adaptation.
/// Geometry sizing and gradation remain active. Periodic size traces must match.
pub fn generate_with_size<F: ScalarField>(
    field: &F,
    size_field: Option<&dyn ScalarField>,
    options: Options,
) -> Result<Output, MeshingError> {
    generate_with_features(field, size_field, options, None)
}

/// Whether finite affine hints certify equal complete max-field traces at the ends.
/// This is conservative: each restricted affine equation must have an exact match.
/// A matching trace does not certify that a generated seed has matching cap facets.
pub fn periodic_affine_traces(
    planes: &[[f64; 4]],
    bounds: [Point; 2],
    periodic: [bool; 3],
) -> bool {
    if !(2..=intersection::MAX_CONSTRAINTS).contains(&planes.len()) {
        return false;
    }
    let full = intersection::ConstraintMask::MAX >> (intersection::MAX_CONSTRAINTS - planes.len());
    periodic_affine_region_traces(planes, bounds, periodic, &[full])
}

/// Conservative complete-region trace certificate for a union of intersections.
/// Each conjunction of normalized restricted planes must have an exact end match.
/// Periodic seed nodes and cap facets still require independent validation.
pub fn periodic_affine_region_traces(
    planes: &[[f64; 4]],
    bounds: [Point; 2],
    periodic: [bool; 3],
    regions: &[intersection::ConstraintMask],
) -> bool {
    if !(2..=intersection::MAX_CONSTRAINTS).contains(&planes.len())
        || regions.is_empty()
        || regions.len() > 32
        || (0..3).any(|a| {
            !bounds[0][a].is_finite() || !bounds[1][a].is_finite() || bounds[0][a] >= bounds[1][a]
        })
    {
        return false;
    }
    let allowed =
        intersection::ConstraintMask::MAX >> (intersection::MAX_CONSTRAINTS - planes.len());
    if regions.iter().any(|&m| m == 0 || m & !allowed != 0) {
        return false;
    }
    crate::polyhedral::normalized(planes).is_some_and(|planes| {
        crate::polyhedral::periodic_traces(&planes, bounds, periodic, Some(regions))
    })
}

/// Generate with optional equivalent max-of-affine hints for adaptive Tet4.
/// The complete plane set must have matching traces on requested periodic ends.
/// Hints must describe the same negative sublevel set. Hints
/// whose normalized coefficients are unrepresentable retain scalar behavior.
pub fn generate_with_features<F: ScalarField>(
    field: &F,
    size_field: Option<&dyn ScalarField>,
    options: Options,
    halfspaces: Option<&[[f64; 4]]>,
) -> Result<Output, MeshingError> {
    generate_with_regions(field, size_field, options, halfspaces, None)
}

/// Generate a Tet4 union of intersections using equivalent affine hints.
/// Periodic regions require adaptation and certified complete end traces.
/// Region masks select nonpositive planes; callers own hint equivalence.
pub fn generate_with_regions<F: ScalarField>(
    field: &F,
    size_field: Option<&dyn ScalarField>,
    options: Options,
    halfspaces: Option<&[[f64; 4]]>,
    regions: Option<&[intersection::ConstraintMask]>,
) -> Result<Output, MeshingError> {
    generate_internal(field, size_field, options, halfspaces, regions, None)
}

/// Adapt an equivalent retained intersection. Periodic and quadratic feature
/// contexts are not admitted by this interface.
pub fn generate_with_constraints<F: ScalarField>(
    field: &F,
    size_field: Option<&dyn ScalarField>,
    options: Options,
    constraints: &[&dyn ScalarField],
) -> Result<Output, MeshingError> {
    generate_internal(field, size_field, options, None, None, Some(constraints))
}

fn generate_internal<F: ScalarField>(
    field: &F,
    size_field: Option<&dyn ScalarField>,
    options: Options,
    halfspaces: Option<&[[f64; 4]]>,
    regions: Option<&[intersection::ConstraintMask]>,
    constraints: Option<&[&dyn ScalarField]>,
) -> Result<Output, MeshingError> {
    if let Some(fields) = constraints
        && (!(2..=intersection::MAX_CONSTRAINTS).contains(&fields.len())
            || options.adapt.is_none()
            || options.periodic.iter().any(|&v| v)
            || options.order != 1
            || options.threads != 1
            || options.interior_refine_nodes != 0
            || options.optimize_passes > 20
            || !matches!(options.region, Region::Negative))
    {
        return Err(MeshingError::InvalidOptions(
            "retained adaptive constraints require bounded nonperiodic Tet4 and one worker".into(),
        ));
    }
    if let Some(groups) = regions {
        let count = halfspaces.map_or(0, <[_]>::len);
        if !(2..=intersection::MAX_CONSTRAINTS).contains(&count)
            || groups.is_empty()
            || groups.len() > 32
            || (options.periodic.iter().any(|&axis| axis) && options.adapt.is_none())
        {
            return Err(MeshingError::InvalidOptions(
                "bounded affine planes and 1..32 regions required; periodic regions need adaptation".into(),
            ));
        }
        let allowed = intersection::ConstraintMask::MAX >> (intersection::MAX_CONSTRAINTS - count);
        if groups.iter().any(|&m| m == 0 || m & !allowed != 0) {
            return Err(MeshingError::InvalidOptions(
                "nonempty valid region masks required".into(),
            ));
        }
    }
    if let Some(planes) = halfspaces
        && (!(2..=intersection::MAX_CONSTRAINTS).contains(&planes.len())
            || planes
                .iter()
                .any(|p| p.iter().any(|x| !x.is_finite()) || p[..3] == [0.; 3])
            || (options.adapt.is_none() && regions.is_none())
            || !matches!(options.region, Region::Negative)
            || options.order != 1
            || options.threads != 1
            || options.interior_refine_nodes != 0
            || options.optimize_passes > 20)
    {
        return Err(MeshingError::InvalidOptions("affine features require bounded finite planes, adaptive Tet4, one worker and no interior-refinement request".into()));
    }
    let halfspaces = halfspaces.and_then(crate::polyhedral::normalized);
    if let Some(planes) = &halfspaces
        && !crate::polyhedral::periodic_traces(planes, options.bounds, options.periodic, regions)
    {
        return Err(MeshingError::InvalidOptions(
            "periodic features require matching affine end traces".into(),
        ));
    }
    if size_field.is_some() && options.adapt.is_none() {
        return Err(MeshingError::InvalidOptions(
            "size_field requires adaptation".into(),
        ));
    }
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
        constraints,
        halfspaces: halfspaces.as_deref(),
        halfspace_regions: halfspaces.as_ref().and(regions),
        size_field,
        options,
        center,
        scale,
        optimizer_length: (0..3).map(|a| hi[a] - lo[a]).product::<f64>().cbrt(),
        invalid: AtomicBool::new(false),
    };
    // Inspect raw field traces before canonicalizing opposite faces.
    for axis in 0..3 {
        if !options.periodic[axis] {
            continue;
        }
        let mut size_samples = Vec::new();
        for i in 0..11 {
            for j in 0..11 {
                let mut p = lo;
                p[(axis + 1) % 3] += (hi[(axis + 1) % 3] - lo[(axis + 1) % 3]) * i as f64 / 10.;
                p[(axis + 2) % 3] += (hi[(axis + 2) % 3] - lo[(axis + 2) % 3]) * j as f64 / 10.;
                let mut q = p;
                q[axis] = hi[axis];
                if size_field.is_some() {
                    size_samples.extend([p, q]);
                }
                let difference = (field.value(p) - field.value(q)) / scale;
                if !difference.is_finite() || difference.abs() > 1e-8 {
                    return Err(MeshingError::InvalidOptions(format!(
                        "field/grading is incompatible with periodic axis {axis}"
                    )));
                }
            }
        }
        if let Some(size) = size_field {
            let values = size.values(&size_samples);
            if values.len() != size_samples.len()
                || values
                    .iter()
                    .any(|value| !value.is_finite() || *value <= 0.)
            {
                return Err(MeshingError::InvalidOptions(
                    "size_field must return finite positive lengths".into(),
                ));
            }
            if values
                .chunks_exact(2)
                .any(|pair| (pair[0] - pair[1]).abs() > 1e-8 * pair[0].max(pair[1]))
            {
                return Err(MeshingError::InvalidOptions(format!(
                    "size_field is incompatible with periodic axis {axis}"
                )));
            }
        }
    }
    context.check().map_err(MeshingError::GenerationFailed)?;
    if !(1..=2).contains(&options.order) {
        return Err(MeshingError::InvalidOptions("order must be 1 or 2".into()));
    }
    let result = build(&context).and_then(|mut output| {
        if options.order == 2 && !output.mesh.tets.is_empty() {
            let elements = crate::quadratic::elevate(
                &output.mesh,
                &output.boundary_tags,
                &context,
                context.threshold(),
                options.periodic,
            )?;
            let periodic_pairs = crate::quadratic::node_pairs(
                &elements.points,
                &elements.surface,
                options.bounds,
                options.periodic,
            )?;
            output.quadratic = Some(QuadraticOutput {
                elements,
                periodic_pairs,
            });
        }
        Ok(output)
    });
    if context.invalid.load(Ordering::Relaxed) {
        return Err(MeshingError::GenerationFailed(
            "field or gradient returned nonfinite values".into(),
        ));
    }
    result.map_err(MeshingError::GenerationFailed)
}

// Every closed surface edge has two incident triangles, and every vertex link
// must be one cycle. Edge counts alone miss two shells touching at one vertex.
pub(crate) fn nonmanifold_vertices(surface: &[[usize; 3]], point_count: usize) -> Vec<usize> {
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

/// Match every material cap node and facet, then verify complete node classes.
/// A wall may touch a box edge/corner without a positive-area cap on that
/// axis. Its empty opposite image needs no node; all existing images must
/// still be connected by the returned translational pairs.
pub(crate) fn periodic_pairs<const N: usize>(
    points: &[Point],
    surface: &[[usize; N]],
    bounds: [Point; 2],
    periodic: [bool; 3],
) -> Result<[Vec<[usize; 2]>; 3], String> {
    fn root(parent: &mut [usize], mut node: usize) -> usize {
        while parent[node] != node {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        node
    }
    let [lo, hi] = bounds;
    let mut pairs: [Vec<[usize; 2]>; 3] = std::array::from_fn(|_| Vec::new());
    let mut parent: Vec<_> = (0..points.len()).collect();
    for (a, pair) in pairs.iter_mut().enumerate() {
        if !periodic[a] {
            continue;
        }
        // Transverse box corners remain distinct within each face.
        let key = |p: Point| {
            std::array::from_fn::<_, 3, _>(|axis| {
                if axis == a {
                    0
                } else {
                    ((p[axis] - lo[axis]) / (hi[axis] - lo[axis]) * 1e11).round() as i64
                }
            })
        };
        let mut maps = [BTreeMap::new(), BTreeMap::new()];
        let mut cap_nodes = [
            std::collections::BTreeSet::new(),
            std::collections::BTreeSet::new(),
        ];
        let mut faces = [Vec::new(), Vec::new()];
        for f in surface {
            for (side, &value) in [lo[a], hi[a]].iter().enumerate() {
                if f[..3].iter().all(|&i| points[i][a] == value) {
                    if !f.iter().all(|&i| points[i][a] == value) {
                        return Err(format!("periodic cap node leaves axis {a} plane"));
                    }
                    cap_nodes[side].extend(f.iter().copied());
                    let mut keys = f.map(|i| key(points[i]));
                    keys.sort_unstable();
                    faces[side].push(keys);
                }
            }
        }
        for (i, &p) in points.iter().enumerate() {
            for (side, &value) in [lo[a], hi[a]].iter().enumerate() {
                if p[a] == value && maps[side].insert(key(p), i).is_some() {
                    return Err("ambiguous periodic node match".into());
                }
            }
        }
        for side in 0..2 {
            if maps[side]
                .iter()
                .any(|(k, i)| !maps[1 - side].contains_key(k) && cap_nodes[side].contains(i))
            {
                return Err(format!("periodic cap node mismatch on axis {a}"));
            }
        }
        for (k, &i) in &maps[0] {
            if let Some(&j) = maps[1].get(k) {
                pair.push([i, j]);
                let ri = root(&mut parent, i);
                let rj = root(&mut parent, j);
                parent[rj] = ri;
            }
        }
        faces.iter_mut().for_each(|f| f.sort_unstable());
        if faces[0] != faces[1] {
            return Err(format!("periodic triangle mismatch on axis {a}"));
        }
    }
    let mut classes = BTreeMap::new();
    for (i, p) in points.iter().enumerate() {
        if !(0..3).any(|a| periodic[a] && (p[a] == lo[a] || p[a] == hi[a])) {
            continue;
        }
        let key = std::array::from_fn::<_, 3, _>(|a| {
            let x = if periodic[a] && p[a] == hi[a] {
                lo[a]
            } else {
                p[a]
            };
            ((x - lo[a]) / (hi[a] - lo[a]) * 1e11).round() as i64
        });
        let representative = root(&mut parent, i);
        if classes
            .insert(key, representative)
            .is_some_and(|old| old != representative)
        {
            return Err("disconnected periodic edge/corner node class".into());
        }
    }
    Ok(pairs)
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
            // A chord lying in a box plane must follow the wall/plane curve inside
            // that plane. The full gradient measures the distance to the unclipped
            // wall, which can pass much closer outside the box, so it under-reports
            // curve error. Only the adaptive path, which refines curve chords by this
            // measure, applies it; the lattice pipeline keeps its released check.
            let in_plane = g.options.adapt.is_some();
            let lock = |q: &[Point]| -> [bool; 3] {
                std::array::from_fn(|a| {
                    in_plane && (q.iter().all(|v| v[a] == lo[a]) || q.iter().all(|v| v[a] == hi[a]))
                })
            };
            for sample in [
                (p[0], lock(&p[..1])),
                (p[1], lock(&p[1..2])),
                (p[2], lock(&p[2..])),
                (
                    std::array::from_fn(|a| (p[0][a] + p[1][a] + p[2][a]) / 3.),
                    [false; 3],
                ),
                (
                    std::array::from_fn(|a| (p[0][a] + p[1][a]) / 2.),
                    lock(&[p[0], p[1]]),
                ),
                (
                    std::array::from_fn(|a| (p[1][a] + p[2][a]) / 2.),
                    lock(&[p[1], p[2]]),
                ),
                (
                    std::array::from_fn(|a| (p[2][a] + p[0][a]) / 2.),
                    lock(&[p[2], p[0]]),
                ),
            ] {
                if g.batched() {
                    samples.push(sample);
                    if samples.len() == 4096 {
                        error = error.max(distance_batch(g, &samples)?);
                        samples.clear();
                    }
                } else {
                    error = error.max(g.distance_in(sample.0, sample.1));
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
    let pairs = periodic_pairs(&mesh.points, &mesh.surface, g.bounds(), g.options.periodic)?;
    let empty = mesh.tets.is_empty();
    if !empty && minimum < g.options.minimum_quality {
        return Err(format!(
            "minimum MMG quality {minimum:e} is below required {:e}; mesh rejected ({} elements below 0.1). Improve resolution/optimization or use a quality remesher",
            g.options.minimum_quality, elements_below_quality_01
        ));
    }
    Ok(Output {
        quadratic: None,
        mesh,
        boundary_tags: tags,
        periodic_pairs: pairs,
        diagnostics: Diagnostics {
            retained_constraints: g.constraints.map_or(0, <[_]>::len),
            polyhedral_constraints: g.halfspaces.map_or(0, <[_]>::len),
            polyhedral_regions: g.halfspace_regions.map_or(0, <[_]>::len),
            empty,
            volume,
            minimum_mmg_quality: if empty { 0. } else { minimum },
            maximum_sampled_surface_error: error,
            elements_below_quality_01,
            reverted_snap_vertices: 0,
            inserted_interior_nodes: 0,
            adapt: None,
        },
    })
}

fn build<F: ScalarField>(g: &Context<F>) -> Result<Output, String> {
    let o = g.options;
    if let Some(fields) = g.constraints {
        let seed = intersection::generate(
            fields,
            |p| p,
            Options {
                adapt: None,
                geometry_tolerance: 0.,
                minimum_quality: 0.,
                ..o
            },
        )?;
        let mut mesh = seed.mesh;
        if mesh.tets.is_empty() {
            return finish(mesh, g);
        }
        let stats = crate::tet_adapt::adapt_with_polish(
            &mut mesh,
            g,
            0.,
            o.periodic,
            &o.adapt.expect("validated retained adaptation"),
            o.max_tetrahedra,
            o.optimize_passes,
            Some(&seed.constraints),
        )?;
        let mut output = finish(mesh, g)?;
        output.diagnostics.adapt = Some(stats);
        output.diagnostics.reverted_snap_vertices = seed.quality.reverted_snap_vertices;
        return Ok(output);
    }
    if let Some(planes) = g.halfspaces {
        struct Affine<'a, F: ScalarField> {
            plane: [f64; 4],
            context: &'a Context<'a, F>,
        }
        impl<F: ScalarField> ScalarField for Affine<'_, F> {
            fn value(&self, p: Point) -> f64 {
                dot([self.plane[0], self.plane[1], self.plane[2]], p) + self.plane[3]
            }
            fn gradient(&self, _: Point) -> Option<Point> {
                Some([self.plane[0], self.plane[1], self.plane[2]])
            }
            fn check(&self) -> Result<(), String> {
                self.context.check()
            }
        }
        let fields: Vec<_> = planes
            .iter()
            .map(|&plane| Affine { plane, context: g })
            .collect();
        let references: Vec<&dyn ScalarField> =
            fields.iter().map(|f| f as &dyn ScalarField).collect();
        let permuting = planes
            .iter()
            .any(|p| (0..3).any(|a| o.periodic[a] && p[a] != 0.));
        let seed_options = Options {
            periodic: if permuting { [false; 3] } else { o.periodic },
            optimize_passes: if permuting { 0 } else { o.optimize_passes },
            adapt: None,
            geometry_tolerance: 0.,
            minimum_quality: 0.,
            ..o
        };
        let seed = if let Some(regions) = g.halfspace_regions {
            intersection::generate_union(&references, regions, seed_options)?
        } else {
            intersection::generate(&references, |p| p, seed_options)?
        };
        if permuting {
            periodic_pairs(&seed.mesh.points, &seed.mesh.surface, o.bounds, o.periodic)?;
        }
        let reverted_snap_vertices = seed.quality.reverted_snap_vertices;
        let mut mesh = seed.mesh;
        if mesh.tets.is_empty() || o.adapt.is_none() {
            let mut output = finish(mesh, g)?;
            output.diagnostics.reverted_snap_vertices = reverted_snap_vertices;
            return Ok(output);
        }
        let stats = crate::tet_adapt::adapt_with_polish(
            &mut mesh,
            g,
            0.,
            o.periodic,
            &o.adapt.expect("validated feature adaptation"),
            o.max_tetrahedra,
            o.optimize_passes,
            None,
        )?;
        let mut output = finish(mesh, g)?;
        output.diagnostics.reverted_snap_vertices = reverted_snap_vertices;
        output.diagnostics.adapt = Some(stats);
        return Ok(output);
    }
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
        // A tiny sampled component can disappear entirely if every interior
        // sample snaps onto its wall. Surface manifold tests cannot see an
        // absent component. Restore discarded interior samples before accepting.
        let mut represented = vec![false; original_points.len()];
        for &source in &sources {
            if source < represented.len() {
                represented[source] = true;
            }
        }
        let lost: Vec<_> = (0..original_points.len())
            .filter(|&i| {
                !represented[i]
                    && (if g.band() {
                        original_values[i].abs() - threshold
                    } else {
                        original_values[i]
                    }) < -1e-12
                    && (b.points[i] != original_points[i] || b.values[i] != original_values[i])
            })
            .collect();
        if bad.is_empty() && lost.is_empty() {
            mesh = Some(candidate);
            break;
        }
        let mut undo = std::collections::BTreeSet::new();
        for i in lost {
            undo.insert(g.rank(original_points[i]));
        }
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
    if let Some(adapt) = o.adapt {
        let mut stats = crate::tet_adapt::adapt(
            &mut mesh,
            g,
            threshold,
            o.periodic,
            &adapt,
            o.max_tetrahedra,
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
        if o.optimize_passes > 0 {
            // Quality improvement can make a blocked wall insertion feasible.
            // Close geometry after that pass without repeating size adaptation.
            let closing = crate::tet_adapt::adapt(
                &mut mesh,
                g,
                threshold,
                o.periodic,
                &crate::tet_adapt::Options {
                    iterations: 0,
                    ..adapt
                },
                o.max_tetrahedra,
            )?;
            stats.splits += closing.splits;
            stats.collapses += closing.collapses;
            stats.swaps += closing.swaps;
            stats.moves += closing.moves;
            stats.error_splits += closing.error_splits;
            stats.maximum_error = closing.maximum_error;
            for (total, elapsed) in stats.seconds.iter_mut().zip(closing.seconds) {
                *total += elapsed;
            }
        }
        let inserted_interior_nodes = crate::adaptive::refine_interior(
            &mut mesh,
            o.interior_refine_nodes,
            o.minimum_quality.max(0.1),
            o.max_tetrahedra,
            g,
        )?;
        let mut output = finish(mesh, g)?;
        output.diagnostics.reverted_snap_vertices = reverted_snap_vertices;
        output.diagnostics.inserted_interior_nodes = inserted_interior_nodes;
        output.diagnostics.adapt = Some(stats);
        return Ok(output);
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
    let inserted_interior_nodes = crate::adaptive::refine_interior(
        &mut mesh,
        o.interior_refine_nodes,
        o.minimum_quality.max(0.1),
        o.max_tetrahedra,
        g,
    )?;
    let mut output = finish(mesh, g)?;
    output.diagnostics.reverted_snap_vertices = reverted_snap_vertices;
    output.diagnostics.inserted_interior_nodes = inserted_interior_nodes;
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

fn distance_batch<F: ScalarField>(
    g: &Context<F>,
    samples: &[(Point, [bool; 3])],
) -> Result<f64, String> {
    g.check()?;
    let points: Vec<Point> = samples.iter().map(|s| s.0).collect();
    let values = g.values(&points);
    let gradients = g.gradients(&points);
    let mut error = 0_f64;
    for ((v, mut grad), (_, locked)) in values.into_iter().zip(gradients).zip(samples) {
        let residual = if g.band() { v.abs() - g.threshold() } else { v };
        for a in 0..3 {
            if locked[a] {
                grad[a] = 0.;
            }
        }
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
    fn schwarz_d_zero_level_periodic_regression() {
        let field = |p: crate::Point| {
            let [x, y, z] = p.map(|a| std::f64::consts::TAU * a);
            x.sin() * y.sin() * z.sin()
                + x.sin() * y.cos() * z.cos()
                + x.cos() * y.sin() * z.cos()
                + x.cos() * y.cos() * z.sin()
        };
        let output = super::generate(
            &field,
            super::Options {
                bounds: [[-0.5; 3], [0.5; 3]],
                cells: [4; 3],
                periodic: [true; 3],
                snap: 0.,
                optimize_passes: 0,
                geometry_tolerance: 1.,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(output.periodic_pairs.iter().all(|pairs| !pairs.is_empty()));
        assert!(output.mesh.points.contains(&[-0.5, 0.5, 0.5]));
        assert!(!output.mesh.points.contains(&[0.5, 0.5, 0.5]));
    }

    #[test]
    fn periodic_contact_nodes_still_require_complete_existing_image_classes() {
        let surface: Vec<[usize; 3]> = Vec::new();
        let mut points = vec![[-0.5, -0.5, 0.1], [0.5, 0.5, 0.1]];
        let check = |points: &[crate::Point]| {
            super::periodic_pairs(points, &surface, [[-0.5; 3], [0.5; 3]], [true, true, false])
        };
        assert!(
            check(&points)
                .unwrap_err()
                .contains("disconnected periodic")
        );
        points.push([0.5, -0.5, 0.1]);
        let pairs = check(&points).unwrap();
        assert_eq!(pairs[0], vec![[0, 2]]);
        assert_eq!(pairs[1], vec![[2, 1]]);
    }

    #[test]
    fn periodic_quadratic_caps_require_every_midpoint_on_the_matching_trace() {
        let lower = [
            [-0.5, -0.2, -0.2],
            [-0.5, 0.2, -0.2],
            [-0.5, 0., 0.2],
            [-0.5, 0., -0.2],
            [-0.5, 0.1, 0.],
            [-0.5, -0.1, 0.],
        ];
        let mut points = lower.to_vec();
        points.extend(lower.map(|p| [0.5, p[1], p[2]]));
        let surface = [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9, 10, 11]];
        let check = |points: &[crate::Point]| {
            super::periodic_pairs(
                points,
                &surface,
                [[-0.5; 3], [0.5; 3]],
                [true, false, false],
            )
        };
        assert_eq!(check(&points).unwrap()[0].len(), 6);
        points[9][1] += 0.001;
        assert!(check(&points).unwrap_err().contains("cap node mismatch"));
        points[9][1] -= 0.001;
        points[9][0] -= 0.001;
        assert!(check(&points).unwrap_err().contains("cap node leaves"));
    }

    #[test]
    fn periodic_cap_gradients_keep_the_physical_side() {
        use super::{Context, Options, ScalarField};
        use crate::geometry::Geometry;
        use std::sync::atomic::AtomicBool;

        struct Taper(bool);
        impl ScalarField for Taper {
            fn value(&self, p: crate::Point) -> f64 {
                p[0] * p[0] + p[1] * p[1] - 0.04 - 0.1 * p[2] * p[2]
            }
            fn gradient(&self, p: crate::Point) -> Option<crate::Point> {
                self.0.then_some([2. * p[0], 2. * p[1], -0.2 * p[2]])
            }
        }
        for analytical in [true, false] {
            let field = Taper(analytical);
            let context = Context {
                field: &field,
                constraints: None,
                halfspaces: None,
                halfspace_regions: None,
                size_field: None,
                options: Options {
                    bounds: [[-0.5; 3], [0.5; 3]],
                    periodic: [false, false, true],
                    ..Default::default()
                },
                center: 0.,
                scale: 1.,
                optimizer_length: 1.,
                invalid: AtomicBool::new(false),
            };
            let points = [[0.2, 0., -0.5], [0.2, 0., 0.5]];
            assert_eq!(context.value(points[0]), context.value(points[1]));
            for (p, batched) in points.into_iter().zip(context.gradients(&points)) {
                let expected = [0.4, 0., -0.2 * p[2]];
                for measured in [context.gradient(p), batched] {
                    for (value, target) in measured.into_iter().zip(expected) {
                        assert!(
                            (value - target).abs() < 1e-9,
                            "{measured:?} != {expected:?}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn periodic_seed_requires_both_bijective_nodes_and_matching_caps() {
        let mut mesh = crate::Mesh {
            points: vec![
                [-0.2, -0.2, -0.5],
                [0.2, -0.2, -0.5],
                [0.2, 0.2, -0.5],
                [-0.2, 0.2, -0.5],
                [-0.2, -0.2, 0.5],
                [0.2, -0.2, 0.5],
                [0.2, 0.2, 0.5],
                [-0.2, 0.2, 0.5],
            ],
            tets: Vec::new(),
            surface: vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
        };
        let check = |mesh: &crate::Mesh| {
            super::periodic_pairs(
                &mesh.points,
                &mesh.surface,
                [[-0.5; 3], [0.5; 3]],
                [false, false, true],
            )
        };
        assert_eq!(check(&mesh).unwrap()[2].len(), 4);
        mesh.surface[2] = [4, 5, 7];
        mesh.surface[3] = [5, 6, 7];
        assert!(check(&mesh).unwrap_err().contains("triangle mismatch"));
        mesh.points[4][0] += 0.01;
        assert!(check(&mesh).unwrap_err().contains("node mismatch"));
    }

    #[test]
    fn complete_affine_trace_certificate_keeps_plane_multiplicities() {
        let planes = [
            [1., 0., 0.3, -0.35],
            [1., 0., -0.3, -0.35],
            [-1., 0., 0.3, -0.35],
            [-1., 0., -0.3, -0.35],
            [0., 1., 0., -0.2],
            [0., -1., 0., -0.2],
        ];
        let check = |p: &[[f64; 4]], b| super::periodic_affine_traces(p, b, [false, false, true]);
        assert!(check(&planes, [[-0.5; 3], [0.5; 3]]));
        assert!(!check(&planes, [[-0.5; 3], [0.6; 3]]));
        let mut wrong = planes;
        wrong[1] = wrong[0];
        assert!(!check(&wrong, [[-0.5; 3], [0.5; 3]]));
        wrong = planes;
        wrong[1][3] -= 0.01;
        assert!(!check(&wrong, [[-0.5; 3], [0.5; 3]]));
        wrong = planes;
        wrong[0][3] = f64::NAN;
        assert!(!check(&wrong, [[-0.5; 3], [0.5; 3]]));
    }

    #[test]
    fn complete_region_traces_preserve_conjunctions_and_branch_multiplicities() {
        let planes = [[1., 0., 0.3, -0.35], [1., -0., -0.3, -0.35]];
        let check = |groups: &[super::intersection::ConstraintMask]| {
            super::periodic_affine_region_traces(
                &planes,
                [[-0.5; 3], [0.5; 3]],
                [false, false, true],
                groups,
            )
        };
        assert!(super::periodic_affine_traces(
            &planes,
            [[-0.5; 3], [0.5; 3]],
            [false, false, true]
        ));
        assert!(check(&[3]));
        assert!(check(&[1, 2]));
        assert!(!check(&[1]));
        assert!(!check(&[1, 1, 2]));
        assert!(!check(&[]));
        assert!(!check(&[0]));
        assert!(!check(&[4]));
        let huge = [[0., 0., 1., 1e308], [0., 0., -1., 1e308]];
        assert!(!super::periodic_affine_region_traces(
            &huge,
            [[-0.5, -0.5, -1e308], [0.5, 0.5, 1e308]],
            [false, false, true],
            &[1, 2]
        ));
    }

    #[test]
    fn periodic_regions_cannot_skip_coupled_adaptation() {
        let planes = [[1., 0., 0.3, -0.35], [1., 0., -0.3, -0.35]];
        let result = super::generate_with_regions(
            &|p: crate::Point| p[0] - 0.3 * p[2].abs() - 0.35,
            None,
            super::Options {
                periodic: [false, false, true],
                ..Default::default()
            },
            Some(&planes),
            Some(&[1, 2]),
        );
        assert!(
            matches!(result, Err(crate::MeshingError::InvalidOptions(message))
            if message.contains("periodic regions need adaptation"))
        );
    }

    #[test]
    fn catches_vertex_pinches_even_when_all_edges_have_two_faces() {
        let surface = crate::boundary(&[[0, 1, 2, 3], [0, 4, 5, 6]]).unwrap();
        assert_eq!(nonmanifold_vertices(&surface, 7), vec![0]);
        let separated = crate::boundary(&[[0, 1, 2, 3], [4, 5, 6, 7]]).unwrap();
        assert!(nonmanifold_vertices(&separated, 8).is_empty());
    }
}

#[cfg(test)]
mod periodic_feature_research;
