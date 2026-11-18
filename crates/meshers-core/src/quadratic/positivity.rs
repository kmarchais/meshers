//! Bound the determinant on a reference tet using the affine variation of its Jacobian.
use super::jacobian_matrix;
use crate::{Point, cross, dot};

type Matrix = [Point; 3];

/// The input coordinates are relative to the first node, including when fitting.
pub(super) fn positive(x: &[Point; 10]) -> bool {
    let corners = std::array::from_fn(|i| {
        let mut l = [0.; 4];
        l[i] = 1.;
        jacobian_matrix(x, l)
    });
    // Scale per physical coordinate, allowing for cancellation in the quadratic
    // derivatives and in the triple products below. Keep this guard on subdivision.
    let scales: Point =
        std::array::from_fn(|axis| x.iter().map(|p| p[axis].abs()).fold(0., f64::max));
    let guard = 4096. * f64::EPSILON * scales.iter().product::<f64>();
    let mut budget = 4096;
    bound(&corners, guard, 0, &mut budget)
}

fn bound(corners: &[Matrix; 4], guard: f64, depth: usize, budget: &mut usize) -> bool {
    if *budget == 0 {
        return false;
    }
    *budget -= 1;
    if corners.iter().any(|m| {
        let d = dot(m[0], cross(m[1], m[2]));
        !d.is_finite() || d <= guard
    }) {
        return false;
    }
    // Each row is a convex combination of its four corner rows. By
    // multilinearity the determinant is a convex combination of these 64
    // triple products (even when the barycentric weights are shared).
    let mut mixed = [[[0.; 4]; 4]; 4];
    for (b, jb) in corners.iter().enumerate() {
        for (c, jc) in corners.iter().enumerate() {
            let normal = cross(jb[1], jc[2]);
            for (a, ja) in corners.iter().enumerate() {
                mixed[a][b][c] = dot(ja[0], normal);
            }
        }
    }
    if mixed.iter().flatten().flatten().all(|&d| d > guard) {
        return true;
    }
    if depth >= 20 {
        return false;
    }
    // Split where the determinant bound varies. Shear can change the matrix
    // substantially along an edge while leaving every triple product unchanged.
    let mut edge = (0, 1);
    let mut largest = 0.;
    for i in 0..4 {
        for j in 0..i {
            let mut difference: f64 = 0.;
            for a in 0..4 {
                for b in 0..4 {
                    difference = difference
                        .max((mixed[i][a][b] - mixed[j][a][b]).abs())
                        .max((mixed[a][i][b] - mixed[a][j][b]).abs())
                        .max((mixed[a][b][i] - mixed[a][b][j]).abs());
                }
            }
            if difference > largest {
                largest = difference;
                edge = (i, j);
            }
        }
    }
    let (i, j) = edge;
    let midpoint = std::array::from_fn(|r| {
        std::array::from_fn(|c| 0.5 * corners[i][r][c] + 0.5 * corners[j][r][c])
    });
    [i, j].into_iter().all(|k| {
        let mut child = *corners;
        child[k] = midpoint;
        bound(&child, guard, depth + 1, budget)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quadratic::{EDGES, jacobian, samples};

    #[test]
    fn fold_between_positive_samples_and_valid_twist() {
        let vertices: [Point; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
        let nodes: [Point; 10] = std::array::from_fn(|i| {
            if i < 4 {
                return vertices[i];
            }
            let [a, b] = EDGES[i - 4];
            std::array::from_fn(|k| (vertices[a][k] + vertices[b][k]) / 2.)
        });
        let folded = nodes.map(|[r, s, t]| [r, (r - 0.38) * s, (r - 0.42) * t]);
        assert!(samples().iter().all(|&l| jacobian(&folded, l) > 0.));
        assert!(!positive(&folded));
        let twist =
            nodes.map(|[r, s, t]| [r, (r - 0.4) * s - 0.025 * t, (r - 0.4) * t + 0.025 * s]);
        assert!(positive(&twist));
        assert!(positive(&nodes));
    }
}
