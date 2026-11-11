//! Experimental optimizer adapter contract. No device dependencies in the core.
//! Factories create an engine per immutable-topology optimization stage. Engines
//! return one displacement per group in color order. Groups in a color have no
//! shared tetrahedra. Implementations with device acceptance must perform the
//! same float64 geometry, inversion and quality checks as CPU acceptance before
//! returning displacements. They are trusted numerical extensions, not a safe
//! way to accept arbitrary unvalidated vertex moves.
use crate::{Mesh, Point};

/// Member node indices, incident tetrahedra, and incident implicit-boundary faces.
pub type Group = (Vec<usize>, Vec<usize>, Vec<[usize; 3]>);

pub trait Factory {
    fn create(
        &self,
        mesh: &Mesh,
        groups: &[Group],
        incident: &[Vec<usize>],
        node_group: &[usize],
    ) -> Result<Box<dyn Optimizer>, String>;
}

pub trait Optimizer {
    fn device_acceptance(&self) -> bool {
        false
    }
    fn proposals(
        &mut self,
        mesh: &Mesh,
        color: &[usize],
        threshold: f64,
        h: f64,
        weight: f64,
        pass: usize,
    ) -> Result<Vec<Point>, String>;
    fn synchronize_updates(&mut self, mesh: &Mesh, nodes: &[usize]) -> Result<(), String>;
}
