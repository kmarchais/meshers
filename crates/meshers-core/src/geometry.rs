//! Internal geometry context shared by specialized and generic algorithms.
use crate::Point;
pub(crate) trait Geometry: Sync {
    fn constraints(&self) -> Option<&[&dyn crate::implicit::ScalarField]> {
        None
    }
    fn halfspaces(&self) -> Option<&[[f64; 4]]> {
        None
    }
    fn halfspace_regions(&self) -> Option<&[crate::implicit::intersection::ConstraintMask]> {
        None
    }
    /// Physical reference length for dimensionless optimizer penalties.
    fn optimizer_length(&self) -> f64 {
        1.
    }
    fn refinement_quality_floor(&self) -> f64 {
        0.25
    }
    fn minimum_quality(&self) -> f64 {
        0.
    }
    fn value(&self, p: Point) -> f64;
    fn batched(&self) -> bool {
        false
    }
    fn values(&self, points: &[Point]) -> Vec<f64> {
        points.iter().map(|&p| self.value(p)).collect()
    }
    fn target_lengths(&self, _points: &[Point]) -> Option<Vec<f64>> {
        None
    }
    fn gradients(&self, points: &[Point]) -> Vec<Point> {
        points.iter().map(|&p| self.gradient(p)).collect()
    }
    fn check(&self) -> Result<(), String> {
        Ok(())
    }
    fn gradient(&self, p: Point) -> Point;
    fn rank(&self, p: Point) -> [i64; 3];
    fn bounds(&self) -> [Point; 2];
    fn band(&self) -> bool {
        true
    }
    /// Implicit pipelines preserve this limit while improving element quality.
    fn geometry_tolerance(&self) -> Option<f64> {
        None
    }
    fn locked(&self, p: Point) -> [bool; 3] {
        let [lo, hi] = self.bounds();
        std::array::from_fn(|a| p[a] == lo[a] || p[a] == hi[a])
    }
    fn contains(&self, p: Point) -> bool {
        let [lo, hi] = self.bounds();
        (0..3).all(|a| p[a] >= lo[a] && p[a] <= hi[a])
    }
    fn residual(&self, p: Point, t: f64) -> f64 {
        let v = self.value(p);
        if self.band() { v.abs() - t } else { v }
    }
    fn level(&self, value: f64, t: f64) -> f64 {
        if self.band() { value.signum() * t } else { 0. }
    }
    /// Distance from each sample to the wall point that Newton steps along its
    /// free gradient reach.
    ///
    /// The first-order estimate `|f| / |grad f|` is a distance only for a distance
    /// function: it misjudges any other field by the curvature of its values. The
    /// length of the path walked to the wall bounds the distance from above
    /// whatever the scaling of the field. Samples that do not converge keep the
    /// first-order estimate, infinite at a vanishing gradient.
    fn wall_distances(
        &self,
        samples: &[(Point, [bool; 3])],
        threshold: f64,
    ) -> Result<Vec<f64>, String> {
        Ok(self
            .wall_walks(samples, threshold)?
            .into_iter()
            .map(|walk| walk.0)
            .collect())
    }
    /// `wall_distances` with the point each walk ended at.
    fn wall_walks(
        &self,
        samples: &[(Point, [bool; 3])],
        threshold: f64,
    ) -> Result<Vec<(f64, Point)>, String> {
        let g = self;
        const STEPS: usize = 8;
        let [lo, hi] = g.bounds();
        let scale = (0..3).map(|a| hi[a] - lo[a]).fold(0., f64::max);
        let resolved = 1e-9 * scale;
        let free = |mut gradient: Point, locked: [bool; 3]| {
            for a in 0..3 {
                if locked[a] {
                    gradient[a] = 0.;
                }
            }
            gradient
        };
        let start: Vec<Point> = samples.iter().map(|s| s.0).collect();
        let values = g.values(&start);
        let gradients = g.gradients(&start);
        let levels: Vec<f64> = values.iter().map(|&v| g.level(v, threshold)).collect();
        let mut distance: Vec<f64> = (0..samples.len())
            .map(|i| wall_distance(values[i] - levels[i], gradients[i], samples[i].1))
            .collect();
        // On the wall already, held at a box corner, or without a usable gradient.
        let mut active: Vec<usize> = (0..samples.len())
            .filter(|&i| {
                distance[i] > resolved && distance[i].is_finite() && samples[i].1 != [true; 3]
            })
            .collect();
        let mut at = start.clone();
        let mut state: Vec<(f64, Point)> = values.into_iter().zip(gradients).collect();
        for _ in 0..STEPS {
            if active.is_empty() {
                break;
            }
            g.check()?;
            for &i in &active {
                let (value, gradient) = state[i];
                let gradient = free(gradient, samples[i].1);
                let slope = crate::norm(gradient);
                let step = ((value - levels[i]) / slope).clamp(-0.05 * scale, 0.05 * scale) / slope;
                at[i] = std::array::from_fn(|a| at[i][a] - step * gradient[a]);
            }
            let moved: Vec<Point> = active.iter().map(|&i| at[i]).collect();
            let values = g.values(&moved);
            let gradients = g.gradients(&moved);
            let mut next = Vec::with_capacity(active.len());
            for ((&i, value), gradient) in active.iter().zip(values).zip(gradients) {
                state[i] = (value, gradient);
                let residual = wall_distance(value - levels[i], gradient, samples[i].1);
                if residual <= resolved {
                    distance[i] = crate::norm(crate::sub(at[i], start[i])) + residual;
                } else if residual.is_finite() {
                    next.push(i);
                }
            }
            active = next;
        }
        Ok(distance.into_iter().zip(at).collect())
    }
    fn cap(&self, points: [Point; 3]) -> bool {
        let [lo, hi] = self.bounds();
        (0..3).any(|a| points.iter().all(|p| p[a] == lo[a]) || points.iter().all(|p| p[a] == hi[a]))
    }
}
pub(crate) fn wall_distance(residual: f64, mut gradient: Point, locked: [bool; 3]) -> f64 {
    let residual = residual.abs();
    if residual == 0. {
        return 0.;
    }
    // Values are normalized by the preflight field amplitude. At a fixed box
    // corner, roundoff at a true wall must not become an infinite 1/0 distance.
    // Retain the full-gradient estimate only for roundoff-sized corner residuals;
    // ordinary misses and genuine zero gradients still reject.
    if locked != [true; 3] || residual > 64. * f64::EPSILON {
        for a in 0..3 {
            if locked[a] {
                gradient[a] = 0.;
            }
        }
    }
    let n = crate::norm(gradient);
    if n > 0. && n.is_finite() {
        residual / n
    } else {
        f64::INFINITY
    }
}
pub(crate) struct Legacy;
impl Geometry for Legacy {
    fn value(&self, p: Point) -> f64 {
        crate::field(p)
    }
    fn gradient(&self, p: Point) -> Point {
        crate::gradient(p)
    }
    fn rank(&self, p: Point) -> [i64; 3] {
        crate::rank(p)
    }
    fn bounds(&self) -> [Point; 2] {
        [[0.; 3], [1.; 3]]
    }
}
