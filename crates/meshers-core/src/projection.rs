//! Projection of a point onto an implicit level set with box coordinates held fixed.
use crate::{Point, dot, geometry::Geometry};

/// Joint Newton projection with normalized equations and rank-aware affine
/// solves. Constraint identities come from the mesh, not proximity guesses.
pub(crate) fn onto_constraints(
    fields: &[&dyn crate::implicit::ScalarField],
    active: u64,
    mut p: Point,
    locked: [bool; 3],
    bounds: [Point; 2],
) -> Option<Point> {
    let scale = (0..3)
        .map(|a| bounds[1][a] - bounds[0][a])
        .fold(0., f64::max);
    for _ in 0..24 {
        let mut planes = Vec::new();
        let mut error = 0f64;
        for (j, field) in fields.iter().enumerate() {
            if active & (1u64 << j) == 0 {
                continue;
            }
            let value = field.value(p);
            let gradient = field.gradient(p).unwrap_or_else(|| {
                std::array::from_fn(|a| {
                    let mut lo = p;
                    let mut hi = p;
                    lo[a] -= scale * 1e-5;
                    hi[a] += scale * 1e-5;
                    (field.value(hi) - field.value(lo)) / (2e-5 * scale)
                })
            });
            let length = crate::norm(gradient);
            if !value.is_finite() || !length.is_finite() {
                return None;
            }
            if length == 0. {
                if value != 0. {
                    return None;
                }
                continue;
            }
            let n = gradient.map(|x| x / length);
            let residual = value / length;
            error = error.max(residual.abs());
            planes.push([n[0], n[1], n[2], residual - dot(n, p)]);
        }
        if error <= scale * 1e-12 {
            return Some(p);
        }
        let mask = (1u64 << planes.len()) - 1;
        let q = crate::polyhedral::project(&planes, p, mask, locked, scale * 1e-12, None)?;
        let step = crate::norm(crate::sub(q, p));
        let fraction = (0.05 * scale / step).min(1.);
        p = std::array::from_fn(|a| {
            if locked[a] {
                p[a]
            } else {
                p[a] + fraction * (q[a] - p[a])
            }
        });
        if p.iter().any(|x| !x.is_finite()) {
            return None;
        }
    }
    None
}

fn free_gradient(g: &impl Geometry, p: Point, locked: [bool; 3]) -> Point {
    let mut d = g.gradient(p);
    for a in 0..3 {
        if locked[a] {
            d[a] = 0.;
        }
    }
    d
}

/// Newton steps along the free gradient. When they do not converge (gradient kinks at
/// creases of min/max fields, walls nearly tangent to a held box plane), bracket the
/// level along the starting free gradient within `reach` and bisect, which converges
/// for any continuous field. The result is not clipped to the box.
pub(crate) fn onto_level(
    g: &impl Geometry,
    start: Point,
    level: f64,
    locked: [bool; 3],
    reach: f64,
) -> Option<Point> {
    let [lo, hi] = g.bounds();
    let scale = (0..3).map(|a| hi[a] - lo[a]).fold(0., f64::max);
    let mut p = start;
    for _ in 0..30 {
        let r = g.value(p) - level;
        let grad = free_gradient(g, p, locked);
        let g2 = dot(grad, grad);
        if !r.is_finite() || g2 <= 1e-300 {
            break;
        }
        if r.abs() / g2.sqrt() <= scale * 1e-13 {
            return Some(p);
        }
        let step = (r / g2).clamp(-0.05 * scale / g2.sqrt(), 0.05 * scale / g2.sqrt());
        p = std::array::from_fn(|a| p[a] - step * grad[a]);
    }
    let r0 = g.value(start) - level;
    if r0 == 0. {
        return Some(start);
    }
    let grad = free_gradient(g, start, locked);
    let norm = dot(grad, grad).sqrt();
    if !r0.is_finite() || norm <= 1e-300 {
        return None;
    }
    let direction: Point = grad.map(|v| -r0.signum() * v / norm);
    let at = |t: f64| -> Point { std::array::from_fn(|a| start[a] + t * direction[a]) };
    let (mut inside, mut outside) = (0., (r0.abs() / norm).max(scale * 1e-9));
    loop {
        let r = g.value(at(outside)) - level;
        if !r.is_finite() {
            return None;
        }
        if r.signum() != r0.signum() {
            break;
        }
        inside = outside;
        outside *= 2.;
        if outside > reach {
            return None;
        }
    }
    for _ in 0..80 {
        let mid = 0.5 * (inside + outside);
        let r = g.value(at(mid)) - level;
        if r == 0. || (outside - inside) <= scale * 1e-14 {
            return Some(at(mid));
        }
        if r.signum() == r0.signum() {
            inside = mid;
        } else {
            outside = mid;
        }
    }
    Some(at(0.5 * (inside + outside)))
}

#[cfg(test)]
mod tests {
    use super::onto_constraints;
    use crate::Point;

    #[test]
    fn joint_nonlinear_projection_handles_scaling_and_redundant_locks() {
        let cylinder = |p: Point| 1e6 * (p[0].hypot(p[1]) - 0.2);
        let cap = |p: Point| 1e-4 * (p[2] - 0.3);
        let bounds = [[-0.5; 3], [0.5; 3]];
        let q = onto_constraints(
            &[&cylinder, &cap],
            3,
            [0.18, 0.03, 0.3],
            [false, false, true],
            bounds,
        )
        .unwrap();
        assert!((q[0].hypot(q[1]) - 0.2).abs() < 1e-12);
        assert_eq!(q[2], 0.3);
        assert!(
            onto_constraints(
                &[&cylinder, &cap],
                3,
                [0.18, 0.03, 0.4],
                [false, false, true],
                bounds
            )
            .is_none()
        );
    }
}
