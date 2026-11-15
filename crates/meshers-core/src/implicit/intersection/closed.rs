//! Identify coincident mapped periodic caps before final quality repair.
use super::*;
use std::collections::BTreeMap;

pub fn generate_closed(
    fields: &[&dyn ScalarField],
    map: impl Fn(Point) -> Point,
    options: Options,
    closed: [bool; 3],
) -> Result<Output, String> {
    if !closed.iter().any(|&v| v) {
        return generate(fields, map, options);
    }
    if (0..3).any(|a| closed[a] && !options.periodic[a]) {
        return Err("closed axes must be periodic".into());
    }
    let mut out = generate(
        fields,
        &map,
        Options {
            minimum_quality: 0.,
            ..options
        },
    )?;
    for axis in 0..3 {
        if !closed[axis] {
            continue;
        }
        let [lo, hi] = options.bounds;
        let key = |p: Point| -> [i64; 3] {
            std::array::from_fn(|a| {
                if a == axis {
                    0
                } else {
                    ((p[a] - lo[a]) / (hi[a] - lo[a]) * 1e11).round() as i64
                }
            })
        };
        let mut low = BTreeMap::new();
        let mut high = BTreeMap::new();
        for (v, &p) in out.parameters.iter().enumerate() {
            let table = if p[axis] == lo[axis] {
                &mut low
            } else if p[axis] == hi[axis] {
                &mut high
            } else {
                continue;
            };
            if table.insert(key(p), v).is_some() {
                return Err("ambiguous closed cap node".into());
            }
        }
        if low.keys().ne(high.keys()) {
            return Err("closed cap parameter nodes do not match".into());
        }
        let extent = (0..3)
            .map(|a| {
                let low = out
                    .mesh
                    .points
                    .iter()
                    .map(|p| p[a])
                    .fold(f64::INFINITY, f64::min);
                let high = out
                    .mesh
                    .points
                    .iter()
                    .map(|p| p[a])
                    .fold(f64::NEG_INFINITY, f64::max);
                high - low
            })
            .fold(0., f64::max);
        let magnitude = out
            .mesh
            .points
            .iter()
            .flatten()
            .map(|v| v.abs())
            .fold(0., f64::max);
        let coincidence = extent * 1e-10 + 8. * f64::EPSILON * magnitude;
        let mut index: Vec<_> = (0..out.mesh.points.len()).collect();
        for (key, &a) in &low {
            let b = high[key];
            if out.constraints[a] != out.constraints[b]
                || norm(sub(out.mesh.points[a], out.mesh.points[b])) > coincidence
                || (0..3).any(|k| {
                    k != axis
                        && (out.parameters[a][k] - out.parameters[b][k]).abs()
                            > (hi[k] - lo[k]) * 1e-11
                })
            {
                return Err("closed cap nodes or constraints do not coincide".into());
            }
            index[b] = a;
        }
        let cap = |bound: f64| -> HashSet<[usize; 3]> {
            out.mesh
                .surface
                .iter()
                .filter(|f| f.iter().all(|&v| out.parameters[v][axis] == bound))
                .map(|f| {
                    let mut f = f.map(|v| index[v]);
                    f.sort_unstable();
                    f
                })
                .collect()
        };
        if cap(lo[axis]) != cap(hi[axis]) {
            return Err("closed cap triangles do not match".into());
        }
        let mut points = Vec::new();
        let mut parameters = Vec::new();
        let mut constraints = Vec::new();
        let mut compact = vec![0; index.len()];
        for v in 0..index.len() {
            if index[v] == v {
                compact[v] = points.len();
                points.push(out.mesh.points[v]);
                parameters.push(out.parameters[v]);
                constraints.push(out.constraints[v]);
            }
        }
        for t in &mut out.mesh.tets {
            *t = t.map(|v| compact[index[v]]);
            if determinant(t.map(|v| points[v])) <= 0. {
                return Err("closed seam creates an invalid cell".into());
            }
        }
        let surface = crate::boundary(&out.mesh.tets)?;
        if !super::super::nonmanifold_vertices(&surface, points.len()).is_empty() {
            return Err("closed seam boundary is not manifold".into());
        }
        out.mesh.points = points;
        out.mesh.surface = surface;
        out.parameters = parameters;
        out.constraints = constraints;
    }
    let wrap = |mut p: Point| {
        for a in 0..3 {
            if closed[a] {
                let lo = options.bounds[0][a];
                let period = options.bounds[1][a] - lo;
                p[a] = lo + (p[a] - lo).rem_euclid(period);
            }
        }
        p
    };
    struct Wrapped<'a> {
        field: &'a dyn ScalarField,
        wrap: &'a (dyn Fn(Point) -> Point + Sync),
    }
    impl ScalarField for Wrapped<'_> {
        fn value(&self, p: Point) -> f64 {
            self.field.value((self.wrap)(p))
        }
        fn gradient(&self, p: Point) -> Option<Point> {
            self.field.gradient((self.wrap)(p))
        }
        fn check(&self) -> Result<(), String> {
            self.field.check()
        }
    }
    let wrapped: Vec<_> = fields
        .iter()
        .map(|&field| Wrapped { field, wrap: &wrap })
        .collect();
    let refs: Vec<&dyn ScalarField> = wrapped.iter().map(|f| f as &dyn ScalarField).collect();
    let mut repair = options;
    for (a, &is_closed) in closed.iter().enumerate() {
        if is_closed {
            let half = (options.bounds[1][a] - options.bounds[0][a]) * 0.5;
            repair.bounds[0][a] -= half;
            repair.bounds[1][a] += half;
            repair.cells[a] *= 2;
            repair.periodic[a] = false;
        }
    }
    // The widened chart only removes artificial cap locks; all fields and the
    // physical map are wrapped onto the unchanged original periodic domain.
    let mut out = super::polish_output(out, &refs, &|p| map(wrap(p)), repair, None)?;
    for p in &mut out.parameters {
        *p = wrap(*p);
    }
    Ok(out)
}
