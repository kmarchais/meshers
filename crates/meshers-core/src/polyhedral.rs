//! Exact affine strata used by the adaptive polyhedral path.
use crate::{Point, dot, norm};

pub(crate) fn normalized(raw: &[[f64; 4]]) -> Option<Vec<[f64; 4]>> {
    raw.iter()
        .map(|p| {
            let scale = p[..3].iter().map(|x| x.abs()).fold(0., f64::max);
            let q = p.map(|x| x / scale);
            let length = norm([q[0], q[1], q[2]]);
            let result = q.map(|x| x / length);
            result.iter().all(|x| x.is_finite()).then_some(result)
        })
        .collect()
}

pub(crate) fn active(planes: &[[f64; 4]], p: Point, tolerance: f64) -> u64 {
    planes.iter().enumerate().fold(0, |mask, (i, plane)| {
        let residual = dot([plane[0], plane[1], plane[2]], p) + plane[3];
        mask | if residual.abs() <= tolerance {
            1u64 << i
        } else {
            0
        }
    })
}

/// Number of independent affine normals, including held box coordinates.
pub(crate) fn constraint_rank(planes: &[[f64; 4]], mask: u64, locked: [bool; 3]) -> usize {
    let mut basis = [[0.; 3]; 3];
    let mut rank = 0;
    let axes = (0..3).filter(|&a| locked[a]).map(|a| {
        let mut n = [0.; 3];
        n[a] = 1.;
        n
    });
    let normals = planes.iter().enumerate().filter_map(|(i, plane)| {
        (mask & (1u64 << i) != 0).then_some([plane[0], plane[1], plane[2]])
    });
    for mut n in axes.chain(normals) {
        for q in &basis[..rank] {
            let factor = dot(n, *q);
            n = std::array::from_fn(|a| n[a] - factor * q[a]);
        }
        let length = norm(n);
        if length > 1e-10 {
            basis[rank] = n.map(|x| x / length);
            rank += 1;
            if rank == 3 {
                break;
            }
        }
    }
    rank
}

/// Orthogonal projection onto active planes and held coordinates. Redundant
/// equations are checked for consistency rather than counted as extra rank.
pub(crate) fn project(
    planes: &[[f64; 4]],
    mut p: Point,
    mask: u64,
    locked: [bool; 3],
    tolerance: f64,
    regions: Option<&[crate::implicit::intersection::ConstraintMask]>,
) -> Option<Point> {
    let mut equations = Vec::new();
    for a in 0..3 {
        if locked[a] {
            let mut n = [0.; 3];
            n[a] = 1.;
            equations.push((n, -p[a]));
        }
    }
    for (i, plane) in planes.iter().enumerate() {
        if mask & (1u64 << i) != 0 {
            equations.push(([plane[0], plane[1], plane[2]], plane[3]));
        }
    }
    let mut basis: Vec<(Point, f64)> = Vec::new();
    for (mut n, mut d) in equations {
        for &(q, e) in &basis {
            let factor = dot(n, q);
            n = std::array::from_fn(|a| n[a] - factor * q[a]);
            d -= factor * e;
        }
        let length = norm(n);
        if length <= 1e-10 {
            if d.abs() > tolerance {
                return None;
            }
        } else {
            basis.push((n.map(|x| x / length), d / length));
        }
    }
    let original = p;
    for (n, d) in basis {
        let residual = dot(n, p) + d;
        p = std::array::from_fn(|a| p[a] - residual * n[a]);
    }
    for a in 0..3 {
        if locked[a] {
            p[a] = original[a];
        }
    }
    let residuals: Vec<_> = planes
        .iter()
        .map(|plane| dot([plane[0], plane[1], plane[2]], p) + plane[3])
        .collect();
    let inside = regions.map_or_else(
        || residuals.iter().all(|&r| r <= tolerance),
        |groups| {
            groups.iter().any(|group| {
                residuals
                    .iter()
                    .enumerate()
                    .all(|(i, &r)| group & (1 << i) == 0 || r <= tolerance)
            })
        },
    );
    let valid = p.iter().all(|x| x.is_finite())
        && inside
        && planes.iter().enumerate().all(|(i, plane)| {
            let residual = dot([plane[0], plane[1], plane[2]], p) + plane[3];
            mask & (1u64 << i) == 0 || residual.abs() <= tolerance
        });
    valid.then_some(p)
}

#[cfg(test)]
mod tests {
    use super::{active, normalized, project};

    #[test]
    fn affine_rank_ignores_redundant_planes_and_box_locks() {
        let planes = normalized(&[
            [1., 0., 0., 0.],
            [0., 1., 0., 0.],
            [0., 0., 1., 0.],
            [1., 1., 0., 0.],
        ])
        .unwrap();
        assert_eq!(super::constraint_rank(&planes, 3, [false; 3]), 2);
        assert_eq!(super::constraint_rank(&planes, 3, [false, false, true]), 3);
        assert_eq!(super::constraint_rank(&planes, 5, [false, false, true]), 2);
        assert_eq!(super::constraint_rank(&planes, 11, [false; 3]), 2);
    }

    #[test]
    fn periodic_hints_reject_noninvariant_planes() {
        use crate::implicit::{Options, generate_with_features};
        let result = generate_with_features(
            &|p: crate::Point| p[0].max(p[2]) - 0.2,
            None,
            Options {
                periodic: [false, false, true],
                adapt: Some(Default::default()),
                ..Options::default()
            },
            Some(&[[1., 0., 0., -0.2], [0., 0., 1., -0.2]]),
        );
        assert!(
            matches!(result, Err(crate::MeshingError::InvalidOptions(message)) if message.contains("planes invariant on periodic axes"))
        );
    }

    #[test]
    fn joint_projection_preserves_strata_locks_and_domain() {
        let planes = normalized(&[
            [1., 0., 0., -0.5],
            [0., 1., 0., -0.5],
            [0., 0., 1., -0.5],
            [1., 1., 1., -1.5],
        ])
        .unwrap();
        let face = project(
            &planes,
            [0.7, 0.1, 0.2],
            1,
            [false, false, true],
            1e-12,
            None,
        )
        .unwrap();
        assert_eq!(face, [0.5, 0.1, 0.2]);
        assert!(project(&planes, [0.7, 0.7, 0.2], 1, [false; 3], 1e-12, None).is_none());
        let corner = project(&planes, [0.7; 3], 15, [false; 3], 1e-12, None).unwrap();
        assert!(corner.into_iter().all(|x| (x - 0.5).abs() < 1e-12));
        assert_eq!(active(&planes, corner, 1e-12), 15);
    }

    #[test]
    fn plane_normalization_handles_extreme_positive_scalings() {
        let planes = normalized(&[[1e300, 0., 0., -5e299], [0., 1e-300, 0., -5e-301]]).unwrap();
        let p = project(&planes, [0.7, 0.7, 0.2], 3, [false; 3], 1e-12, None).unwrap();
        assert_eq!(p, [0.5, 0.5, 0.2]);
        assert!(normalized(&[[1e-300, 0., 0., 1e300]]).is_none());
    }

    #[test]
    fn nonconvex_projection_admits_either_branch_but_rejects_the_notch() {
        let planes = [[1., 0., 0., -0.05], [0., 1., 0., -0.05]];
        let regions = [1, 2];
        let p = [0.1, -0.1, 0.];
        assert_eq!(
            project(&planes, p, 0, [false; 3], 1e-12, Some(&regions)),
            Some(p)
        );
        assert!(project(&planes, p, 0, [false; 3], 1e-12, None).is_none());
        assert!(
            project(
                &planes,
                [0.1, 0.1, 0.],
                0,
                [false; 3],
                1e-12,
                Some(&regions)
            )
            .is_none()
        );
        let p = project(
            &planes,
            [0.1, 0.1, 0.],
            1,
            [false; 3],
            1e-12,
            Some(&regions),
        )
        .unwrap();
        assert!((p[0] - 0.05).abs() < 1e-12 && p[1] == 0.1);
    }
}
