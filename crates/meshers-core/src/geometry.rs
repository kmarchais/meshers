//! Internal geometry context shared by specialized and generic algorithms.
use crate::Point;
pub(crate) trait Geometry: Sync {
    fn value(&self, p: Point) -> f64;
    fn batched(&self) -> bool {
        false
    }
    fn values(&self, points: &[Point]) -> Vec<f64> {
        points.iter().map(|&p| self.value(p)).collect()
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
    fn cap(&self, points: [Point; 3]) -> bool {
        let [lo, hi] = self.bounds();
        (0..3).any(|a| points.iter().all(|p| p[a] == lo[a]) || points.iter().all(|p| p[a] == hi[a]))
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
