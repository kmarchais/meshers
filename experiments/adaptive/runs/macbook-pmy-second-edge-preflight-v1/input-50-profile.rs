//! Optional coarse phase timings; no clocks or allocations when disabled.
use std::time::Instant;
pub(crate) struct Profile {
    scope: &'static str,
    last: Option<Instant>,
    phases: Vec<(&'static str, f64)>,
    counts: Vec<(&'static str, usize)>,
}
impl Profile {
    pub(crate) fn new(scope: &'static str) -> Self {
        Self {
            scope,
            last: std::env::var_os("MESHER_CPU_PROFILE").map(|_| Instant::now()),
            phases: Vec::new(),
            counts: Vec::new(),
        }
    }
    pub(crate) fn mark(&mut self, name: &'static str) {
        if let Some(last) = self.last {
            let now = Instant::now();
            self.phases
                .push((name, now.duration_since(last).as_secs_f64()));
            self.last = Some(now);
        }
    }
    pub(crate) fn count(&mut self, name: &'static str, value: usize) {
        if self.last.is_some() {
            self.counts.push((name, value));
        }
    }
}
impl Drop for Profile {
    fn drop(&mut self) {
        if self.last.is_some() {
            self.mark("cleanup_and_other");
            eprintln!(
                "CPU_PROFILE {}",
                serde_json::json!({"scope":self.scope,"phases":self.phases,"counts":self.counts})
            );
        }
    }
}
