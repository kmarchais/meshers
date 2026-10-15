//! Isolated direct-surface research API; not included in published wheels.
use crate::{CancellationToken, MeshingFailure, make_field};
use meshers_core::{MeshingError, Point, quads::Band, remesh, triangles};
use numpy::{IntoPyArray, ndarray::Array2};
use pyo3::{exceptions::PyValueError, prelude::*, types::PyDict};

#[allow(clippy::too_many_arguments)]
#[pyfunction(signature = (field, gradient, bounds, cells, band, smoothing_iterations, improvement_rounds, polish_passes, batch_size, cancel, check_geometry=true, adaptive=None))]
pub(super) fn generate_surface<'py>(
    py: Python<'py>,
    field: Py<PyAny>,
    gradient: Option<Py<PyAny>>,
    bounds: [Point; 2],
    cells: usize,
    band: [f64; 2],
    smoothing_iterations: usize,
    improvement_rounds: usize,
    polish_passes: usize,
    batch_size: usize,
    cancel: Option<PyRef<'py, CancellationToken>>,
    check_geometry: bool,
    adaptive: Option<(f64, f64, f64, f64, usize, bool, bool)>,
) -> PyResult<Bound<'py, PyDict>> {
    if smoothing_iterations > 100 || improvement_rounds > 100 || polish_passes > 100 {
        return Err(PyValueError::new_err("Surface passes must be at most 100"));
    }
    let f = make_field(
        py,
        field,
        gradient,
        bounds,
        cancel.as_ref().map(|v| v.flag.clone()).unwrap_or_default(),
        batch_size,
    )?;
    let geometry = Band {
        field: &f,
        bounds,
        levels: band,
    };
    let result = py.detach(|| {
        let start = std::time::Instant::now();
        let mut mesh = triangles::extract(&geometry, cells)?;
        let extraction = start.elapsed().as_secs_f64();
        let mut stats = None;
        if let Some((tolerance, min_length, max_length, gradation, iterations, balanced, fast)) =
            adaptive
        {
            let options = remesh::Options {
                tolerance,
                min_length,
                max_length,
                gradation,
                iterations,
                balanced,
                fast,
            };
            stats = Some(remesh::adapt(&mut mesh, &geometry, &options)?);
        } else {
            triangles::smooth(&mut mesh, &geometry, smoothing_iterations)?;
            triangles::improve(&mut mesh, &geometry, improvement_rounds)?;
            triangles::polish(&mut mesh, &geometry, polish_passes)?;
        }
        let (minimum_shape, volume) = triangles::validate(&mesh)?;
        let error = if check_geometry {
            Some(triangles::sampled_error(&mesh, &geometry)?)
        } else {
            None
        };
        Ok::<_, MeshingError>((
            mesh,
            extraction,
            start.elapsed().as_secs_f64(),
            minimum_shape,
            volume,
            error,
            stats,
        ))
    });
    if let Some(error) = f.error.lock().unwrap().take() {
        return Err(error);
    }
    py.check_signals()?;
    let (mesh, extraction, total, minimum_shape, volume, error, stats) =
        result.map_err(|error| match error {
            MeshingError::InvalidOptions(message) => PyValueError::new_err(message),
            MeshingError::GenerationFailed(message) => MeshingFailure::new_err(message),
        })?;
    let output = PyDict::new(py);
    output.set_item(
        "points",
        Array2::from_shape_vec(
            (mesh.points.len(), 3),
            mesh.points.into_iter().flatten().collect(),
        )
        .unwrap()
        .into_pyarray(py),
    )?;
    output.set_item(
        "triangles",
        Array2::from_shape_vec(
            (mesh.faces.len(), 3),
            mesh.faces.into_iter().flatten().map(|v| v as u64).collect(),
        )
        .unwrap()
        .into_pyarray(py),
    )?;
    output.set_item("labels", mesh.labels.into_pyarray(py))?;
    output.set_item("extraction_seconds", extraction)?;
    output.set_item("compute_seconds", total)?;
    output.set_item("minimum_triangle_shape", minimum_shape)?;
    output.set_item("volume", volume)?;
    output.set_item("maximum_sampled_surface_error", error)?;
    output.set_item("geometry_checked", check_geometry)?;
    if let Some(stats) = stats {
        let adaptive = PyDict::new(py);
        adaptive.set_item("splits", stats.splits)?;
        adaptive.set_item("collapses", stats.collapses)?;
        adaptive.set_item("flips", stats.flips)?;
        adaptive.set_item("moves", stats.moves)?;
        adaptive.set_item("error_splits", stats.error_splits)?;
        adaptive.set_item("maximum_error", stats.maximum_error)?;
        let phases = ["sizing", "split", "collapse", "flip", "relax", "enforce"];
        let seconds = PyDict::new(py);
        for (phase, s) in phases.iter().zip(stats.seconds) {
            seconds.set_item(phase, s)?;
        }
        adaptive.set_item("seconds", seconds)?;
        output.set_item("adaptive", adaptive)?;
    }
    Ok(output)
}
