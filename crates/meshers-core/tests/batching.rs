use meshers_core::{
    Point,
    implicit::{self, Options, Region, ScalarField},
};
use std::sync::atomic::{AtomicUsize, Ordering};

struct Field {
    batch: bool,
    largest: AtomicUsize,
}
impl ScalarField for Field {
    fn value(&self, p: Point) -> f64 {
        meshers_core::field(p)
    }
    fn prefers_batches(&self) -> bool {
        self.batch
    }
    fn values(&self, p: &[Point]) -> Vec<f64> {
        self.largest.fetch_max(p.len(), Ordering::Relaxed);
        p.iter().map(|&p| self.value(p)).collect()
    }
}

#[test]
fn batched_optimizer_preserves_scalar_mesh() {
    for threads in [0, 1, 4] {
        let options = Options {
            cells: [8; 3],
            region: Region::Band {
                lower: -0.5,
                upper: 0.5,
            },
            periodic: [true; 3],
            geometry_tolerance: 0.1,
            optimize_passes: 1,
            threads,
            ..Options::default()
        };
        let a = Field {
            batch: false,
            largest: AtomicUsize::new(0),
        };
        let b = Field {
            batch: true,
            largest: AtomicUsize::new(0),
        };
        let old = implicit::generate(&a, options).unwrap();
        let new = implicit::generate(&b, options).unwrap();
        assert_eq!(
            old.mesh.points, new.mesh.points,
            "coordinates threads={threads}"
        );
        assert_eq!(old.mesh.tets, new.mesh.tets);
        assert_eq!(old.mesh.surface, new.mesh.surface);
        assert_eq!(old.periodic_pairs, new.periodic_pairs);
        assert_eq!(
            old.diagnostics.maximum_sampled_surface_error,
            new.diagnostics.maximum_sampled_surface_error
        );
        assert!(b.largest.load(Ordering::Relaxed) > 8);
    }
}

struct Sphere(bool);
impl ScalarField for Sphere {
    fn value(&self, p: Point) -> f64 {
        ((p[0] - 0.47).powi(2) + (p[1] - 0.51).powi(2) + (p[2] - 0.49).powi(2)).sqrt() - 0.29
    }
    fn prefers_batches(&self) -> bool {
        self.0
    }
}
#[test]
fn nonperiodic_sphere_batches_preserve_mesh() {
    for threads in [0, 1] {
        let options = Options {
            cells: [10; 3],
            geometry_tolerance: 0.03,
            optimize_passes: 2,
            threads,
            ..Options::default()
        };
        let a = implicit::generate(&Sphere(false), options).unwrap();
        let b = implicit::generate(&Sphere(true), options).unwrap();
        assert_eq!(a.mesh.points, b.mesh.points);
        assert_eq!(a.mesh.tets, b.mesh.tets);
        assert_eq!(a.mesh.surface, b.mesh.surface);
        assert_eq!(
            a.diagnostics.maximum_sampled_surface_error,
            b.diagnostics.maximum_sampled_surface_error
        );
    }
}
