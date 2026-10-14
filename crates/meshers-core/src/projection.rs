//! Projection of a point onto an implicit level set with box coordinates held fixed.
use crate::{Point, dot, geometry::Geometry};

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
