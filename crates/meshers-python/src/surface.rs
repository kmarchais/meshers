use meshers_core::{MeshingError, Point, surface_band::Band, triangles};
use numpy::{IntoPyArray, ndarray::Array2};
use pyo3::{exceptions::PyValueError, prelude::*, types::PyDict};

use crate::{CancellationToken, MeshingFailure, make_field};

#[allow(clippy::too_many_arguments)]
#[pyfunction]
pub(super) fn generate_surface<'py>(
    py: Python<'py>,
    field: Py<PyAny>,
    gradient: Option<Py<PyAny>>,
    bounds: [Point; 2],
    cells: usize,
    band: [f64; 2],
    periodic: [bool; 3],
    smoothing_iterations: usize,
    improvement_rounds: usize,
    polish_passes: usize,
    batch_size: usize,
    cancel: Option<PyRef<'py, CancellationToken>>,
) -> PyResult<Bound<'py, PyDict>> {
    if smoothing_iterations > 100 || improvement_rounds > 100 || polish_passes > 100 {
        return Err(PyValueError::new_err("surface passes must be at most 100"));
    }
    if periodic.iter().any(|&v| v) && (smoothing_iterations != 0 || improvement_rounds != 0) {
        return Err(PyValueError::new_err(
            "periodic surfaces currently require zero smoothing and improvement rounds",
        ));
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
    let start = std::time::Instant::now();
    let result = py.detach(|| {
        let mut mesh = triangles::extract(&geometry, cells)?;
        if periodic.iter().any(|&v| v) {
            triangles::polish_periodic(&mut mesh, &geometry, periodic, polish_passes)?;
        } else {
            triangles::smooth(&mut mesh, &geometry, smoothing_iterations)?;
            triangles::improve(&mut mesh, &geometry, improvement_rounds)?;
            triangles::polish(&mut mesh, &geometry, polish_passes)?;
        }
        Ok::<_, MeshingError>(mesh)
    });
    if let Some(error) = f.error.lock().unwrap().take() {
        return Err(error);
    }
    py.check_signals()?;
    let result = result.map_err(|e| match e {
        MeshingError::InvalidOptions(message) => PyValueError::new_err(message),
        MeshingError::GenerationFailed(message) => MeshingFailure::new_err(message),
    })?;
    let output = PyDict::new(py);
    output.set_item(
        "points",
        Array2::from_shape_vec(
            (result.points.len(), 3),
            result.points.into_iter().flatten().collect(),
        )
        .unwrap()
        .into_pyarray(py),
    )?;
    output.set_item(
        "triangles",
        Array2::from_shape_vec(
            (result.faces.len(), 3),
            result
                .faces
                .into_iter()
                .flatten()
                .map(|v| v as i64)
                .collect(),
        )
        .unwrap()
        .into_pyarray(py),
    )?;
    output.set_item("labels", result.labels.into_pyarray(py))?;
    output.set_item("seconds", start.elapsed().as_secs_f64())?;
    Ok(output)
}
