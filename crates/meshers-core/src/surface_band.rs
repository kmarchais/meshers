//! Scalar band and box constraints shared by the direct triangle extractor.
use crate::{MeshingError, Point, implicit::ScalarField, norm, sub};

type Result<T> = std::result::Result<T, MeshingError>;

fn failed(message: &str) -> MeshingError {
    MeshingError::GenerationFailed(message.into())
}

fn mul(p: Point, scale: f64) -> Point {
    p.map(|v| v * scale)
}

fn unit(p: Point) -> Point {
    mul(p, 1. / norm(p).max(1e-300))
}

/// An implicit band intersected with an axis-aligned box.
pub struct Band<'a, F: ScalarField> {
    pub field: &'a F,
    pub bounds: [Point; 2],
    pub levels: [f64; 2],
}

impl<F: ScalarField> Band<'_, F> {
    pub(crate) fn validate(&self) -> Result<()> {
        if !(self.levels[0].is_finite()
            && self.levels[1].is_finite()
            && self.levels[0] < self.levels[1])
            || (0..3).any(|i| {
                !self.bounds[0][i].is_finite()
                    || !self.bounds[1][i].is_finite()
                    || self.bounds[0][i] >= self.bounds[1][i]
            })
        {
            return Err(MeshingError::InvalidOptions(
                "Finite ordered band levels and bounds are required".into(),
            ));
        }
        if !self.length().is_finite() || self.length() <= 0. {
            return Err(MeshingError::InvalidOptions(
                "Box extent must have a finite positive length".into(),
            ));
        }
        Ok(())
    }

    fn length(&self) -> f64 {
        norm(sub(self.bounds[1], self.bounds[0]))
    }

    fn gradient(&self, p: Point) -> Point {
        self.field.gradient(p).unwrap_or_else(|| {
            let h = self.length() * 1e-6;
            std::array::from_fn(|i| {
                let mut a = p;
                let mut b = p;
                a[i] += h;
                b[i] -= h;
                (self.field.value(a) - self.field.value(b)) / (2. * h)
            })
        })
    }

    pub(crate) fn normal(&self, p: Point, label: u8) -> Point {
        if label < 2 {
            mul(unit(self.gradient(p)), if label == 0 { 1. } else { -1. })
        } else {
            let mut n = [0.; 3];
            n[(label as usize - 2) / 2] = if label.is_multiple_of(2) { -1. } else { 1. };
            n
        }
    }

    pub(crate) fn project(&self, mut p: Point, mask: u8) -> Result<Point> {
        let free: Point = std::array::from_fn(|i| {
            if (mask >> (2 + 2 * i)) & 3 == 0 {
                1.
            } else {
                0.
            }
        });
        for (i, coordinate) in p.iter_mut().enumerate() {
            for side in 0..2 {
                if mask & (1 << (2 + 2 * i + side)) != 0 {
                    *coordinate = self.bounds[side][i];
                }
            }
        }
        if mask & 3 != 0 {
            let target = self.levels[usize::from(mask & 1 != 0)];
            let mut converged = false;
            for _ in 0..20 {
                let residual = self.field.value(p) - target;
                let g = self.gradient(p);
                let g = std::array::from_fn(|i| g[i] * free[i]);
                let gn = norm(g);
                if !residual.is_finite() || !gn.is_finite() || gn <= f64::MIN_POSITIVE {
                    return Err(failed("Singular/non-finite field during projection"));
                }
                let distance = residual / gn;
                if distance.abs() <= self.length() * 1e-12 {
                    converged = true;
                    break;
                }
                p = sub(
                    p,
                    mul(
                        unit(g),
                        distance.clamp(-self.length() * 0.05, self.length() * 0.05),
                    ),
                );
            }
            if !converged {
                return Err(failed("Surface projection did not converge"));
            }
        }
        if (0..3).any(|i| {
            !p[i].is_finite()
                || p[i] < self.bounds[0][i] - self.length() * 1e-10
                || p[i] > self.bounds[1][i] + self.length() * 1e-10
        }) {
            return Err(failed("Projection left the clipping box"));
        }
        Ok(p)
    }
}
