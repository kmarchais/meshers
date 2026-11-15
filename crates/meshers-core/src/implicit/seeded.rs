//! Adapt a supplied solid mesh while retaining its named physical constraints.
use super::*;

struct Maximum<'a>(&'a [&'a dyn ScalarField]);
impl ScalarField for Maximum<'_> {
    fn value(&self, p: Point) -> f64 {
        let mut maximum = f64::NEG_INFINITY;
        for field in self.0 {
            let value = field.value(p);
            if !value.is_finite() {
                return f64::NAN;
            }
            maximum = maximum.max(value);
        }
        maximum
    }
    fn gradient(&self, p: Point) -> Option<Point> {
        self.0
            .iter()
            .max_by(|a, b| a.value(p).total_cmp(&b.value(p)))?
            .gradient(p)
    }
    fn check(&self) -> Result<(), String> {
        for field in self.0 {
            field.check()?;
        }
        Ok(())
    }
}

/// Adapt a conforming nonperiodic Tet4 seed using the same physical fields.
/// Mask bits identify the supplied fields. The seed must pass a separate sampled
/// `seed_tolerance`; final quality and geometry use `options` without relaxation.
/// The caller owns the seed's topology and resolution of disconnected features.
pub fn adapt_intersection(
    fields: &[&dyn ScalarField],
    seed: &Mesh,
    masks: &[intersection::ConstraintMask],
    options: Options,
    seed_tolerance: f64,
) -> Result<(Output, Vec<intersection::ConstraintMask>), MeshingError> {
    let invalid = |message: &str| MeshingError::InvalidOptions(message.to_owned());
    let failed = MeshingError::GenerationFailed;
    let [lo, hi] = options.bounds;
    if !(2..=intersection::MAX_CONSTRAINTS).contains(&fields.len())
        || options.periodic.iter().any(|&v| v)
        || options.order != 1
        || options.threads != 1
        || options.optimize_passes > 20
        || options.interior_refine_nodes != 0
        || options.curved_fit
        || !matches!(options.region, Region::Negative)
        || options.adapt.is_none()
        || !seed_tolerance.is_finite()
        || seed_tolerance <= 0.
        || !options.geometry_tolerance.is_finite()
        || options.geometry_tolerance <= 0.
        || !options.minimum_quality.is_finite()
        || !(0.0..=1.0).contains(&options.minimum_quality)
        || options.max_tetrahedra == 0
        || (0..3).any(|a| {
            !lo[a].is_finite()
                || !hi[a].is_finite()
                || hi[a] <= lo[a]
                || !(hi[a] - lo[a]).is_finite()
        })
    {
        return Err(invalid(
            "seed adaptation requires finite nonperiodic bounds, named constraints, Tet4, one worker and positive tolerances/budget",
        ));
    }
    if seed.points.len() > options.max_tetrahedra
        || seed.tets.len() > options.max_tetrahedra
        || masks.len() != seed.points.len()
    {
        return Err(invalid(
            "seed counts or mask length exceed the declared bounds",
        ));
    }
    let allowed =
        intersection::ConstraintMask::MAX >> (intersection::MAX_CONSTRAINTS - fields.len());
    if masks.iter().any(|m| m & !allowed != 0)
        || seed
            .points
            .iter()
            .any(|p| (0..3).any(|a| !p[a].is_finite() || p[a] < lo[a] || p[a] > hi[a]))
    {
        return Err(invalid(
            "seed has invalid constraint bits or coordinates outside its bounds",
        ));
    }
    for t in &seed.tets {
        if t.iter().any(|&v| v >= seed.points.len()) {
            return Err(invalid("seed has an invalid tetrahedron index"));
        }
        let det = determinant(t.map(|v| seed.points[v]));
        if !det.is_finite() || det <= 0. {
            return Err(invalid("seed has an invalid or inverted tetrahedron"));
        }
    }
    let surface = crate::boundary(&seed.tets).map_err(failed)?;
    if !nonmanifold_vertices(&surface, seed.points.len()).is_empty() {
        return Err(invalid("seed boundary is not a closed vertex manifold"));
    }
    let field = Maximum(fields);
    let mut scale: f64 = 0.;
    for x in [0., 0.23, 0.61, 1.] {
        for y in [0., 0.23, 0.61, 1.] {
            for z in [0., 0.23, 0.61, 1.] {
                let p = std::array::from_fn(|a| lo[a] + [x, y, z][a] * (hi[a] - lo[a]));
                let value = field.value(p);
                if !value.is_finite() {
                    return Err(invalid("nonfinite field in seed preflight"));
                }
                scale = scale.max(value.abs());
            }
        }
    }
    if scale == 0. {
        return Err(invalid("seed field is unresolved on preflight samples"));
    }
    let g = Context {
        field: &field,
        constraints: Some(fields),
        halfspaces: None,
        halfspace_regions: None,
        size_field: None,
        background: None,
        #[cfg(feature = "experimental-quads")]
        linear_observer: None,
        options,
        center: 0.,
        scale,
        optimizer_length: (0..3).map(|a| hi[a] - lo[a]).product::<f64>().cbrt(),
        invalid: AtomicBool::new(false),
    };
    let weights = [
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
        [0.5, 0.5, 0.],
        [0.5, 0., 0.5],
        [0., 0.5, 0.5],
        [1. / 3.; 3],
    ];
    for (i, face) in surface.iter().enumerate() {
        if i % 512 == 0 {
            g.check().map_err(failed)?;
        }
        let active = face
            .iter()
            .fold(u64::MAX, |mask, &v| mask & u64::from(masks[v]));
        let box_cap = (0..3).any(|a| {
            face.iter().all(|&v| seed.points[v][a] == lo[a])
                || face.iter().all(|&v| seed.points[v][a] == hi[a])
        });
        if active == 0 && !box_cap {
            return Err(invalid(
                "seed boundary facet has no named supporting constraint",
            ));
        }
        for weight in weights {
            let p =
                std::array::from_fn(|a| (0..3).map(|k| weight[k] * seed.points[face[k]][a]).sum());
            let distance = if active == 0 {
                crate::projection::positive_constraint_distance(fields, p, options.bounds)
            } else {
                crate::projection::constraint_distance(
                    fields,
                    active,
                    p,
                    [false; 3],
                    options.bounds,
                )
            };
            if !distance.is_finite() || distance > seed_tolerance {
                return Err(invalid(
                    "seed fails its declared sampled physical geometry tolerance",
                ));
            }
        }
    }
    let mut mesh = Mesh {
        points: seed.points.clone(),
        tets: seed.tets.clone(),
        surface,
    };
    let (stats, masks) = crate::tet_adapt::adapt_with_constraints(
        &mut mesh,
        &g,
        0.,
        [false; 3],
        &options.adapt.unwrap(),
        options.max_tetrahedra,
        options.optimize_passes,
        Some(masks),
    )
    .map_err(failed)?;
    let masks =
        masks.ok_or_else(|| failed("adaptation lost retained constraint identities".into()))?;
    let mut out = finish_with_constraints(mesh, &g, Some(&masks)).map_err(failed)?;
    g.check().map_err(failed)?;
    if g.invalid.load(Ordering::Relaxed) {
        return Err(failed("nonfinite field during seed adaptation".into()));
    }
    out.diagnostics.adapt = Some(stats);
    Ok((out, masks))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn named_seed_adaptation_preserves_labels_and_rejects_an_unclipped_box() {
        let x = |p: Point| p[0] - 0.5;
        let y = |p: Point| p[1] - 0.75;
        let seed_options = Options {
            cells: [4; 3],
            snap: 0.,
            optimize_passes: 0,
            geometry_tolerance: 0.,
            ..Default::default()
        };
        let seed = intersection::generate(&[&x, &y], |p| p, seed_options).unwrap();
        let points = seed.mesh.points.clone();
        let tets = seed.mesh.tets.clone();
        let options = Options {
            geometry_tolerance: 1e-8,
            minimum_quality: 0.1,
            optimize_passes: 1,
            max_tetrahedra: 10000,
            adapt: Some(crate::tet_adapt::Options {
                tolerance: 1e-8,
                min_length: 0.1,
                max_length: 0.4,
                gradation: 0.6,
                iterations: 2,
                quality_guard: 0.1,
                ..Default::default()
            }),
            ..Default::default()
        };
        let (out, masks) =
            adapt_intersection(&[&x, &y], &seed.mesh, &seed.constraints, options, 1e-8).unwrap();
        assert_eq!(masks.len(), out.mesh.points.len());
        assert!((out.diagnostics.volume - 0.375).abs() < 1e-9);
        assert!(out.diagnostics.minimum_mmg_quality >= 0.1);
        for (&p, &m) in out.mesh.points.iter().zip(&masks) {
            if m & 1 != 0 {
                assert!(x(p).abs() < 1e-10);
            }
            if m & 2 != 0 {
                assert!(y(p).abs() < 1e-10);
            }
        }
        assert_eq!(points, seed.mesh.points);
        assert_eq!(tets, seed.mesh.tets);
        let negative = |_: Point| -1.;
        let cube = intersection::generate(&[&negative, &negative], |p| p, seed_options).unwrap();
        let error = adapt_intersection(&[&x, &y], &cube.mesh, &cube.constraints, options, 1e-8)
            .err()
            .unwrap();
        assert!(error.to_string().contains("seed fails"));
    }
    #[test]
    fn seed_preflight_rejects_overflowing_determinants() {
        let field = |_: Point| -1.;
        let seed = Mesh {
            points: vec![[0.; 3], [1e200, 0., 0.], [0., 1e200, 0.], [0., 0., 1e200]],
            tets: vec![[0, 1, 2, 3]],
            surface: vec![],
        };
        let options = Options {
            bounds: [[0.; 3], [1e200; 3]],
            adapt: Some(Default::default()),
            ..Default::default()
        };
        let error = adapt_intersection(&[&field, &field], &seed, &[0; 4], options, 0.01)
            .err()
            .unwrap();
        assert!(
            error
                .to_string()
                .contains("invalid or inverted tetrahedron")
        );
    }
}
