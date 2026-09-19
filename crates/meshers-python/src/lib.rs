mod compiled;
mod intersection;
use meshers_core::{
    MeshingError, Point,
    implicit::{self, ScalarField},
};
use numpy::{IntoPyArray, PyReadonlyArray1, PyReadonlyArray2, ndarray::Array2};
use pyo3::{
    exceptions::{PyRuntimeError, PyTypeError, PyValueError},
    prelude::*,
    types::PyDict,
};
use std::f64::consts::TAU;
use std::sync::{
    Arc, Mutex,
    atomic::{AtomicBool, AtomicUsize, Ordering},
};

pyo3::create_exception!(_meshers, MeshingFailure, PyRuntimeError);
pyo3::create_exception!(_meshers, CancelledError, PyRuntimeError);

#[pyclass(skip_from_py_object)]
#[derive(Clone, Default)]
struct CancellationToken {
    flag: Arc<AtomicBool>,
}
#[pymethods]
impl CancellationToken {
    #[new]
    fn new() -> Self {
        Self::default()
    }
    fn cancel(&self) {
        self.flag.store(true, Ordering::Relaxed);
    }
    #[getter]
    fn cancelled(&self) -> bool {
        self.flag.load(Ordering::Relaxed)
    }
}

struct Field {
    compiled: Option<(compiled::Kernel, bool)>,
    native: Option<String>,
    callback: Option<Py<PyAny>>,
    derivative: Option<Py<PyAny>>,
    error: Mutex<Option<PyErr>>,
    cancellation: Arc<AtomicBool>,
    step: Point,
    calls: AtomicUsize,
    checks: AtomicUsize,
    batch_size: usize,
}
fn native_value(name: &str, p: Point) -> f64 {
    let [x, y, z] = p.map(|v| TAU * v);
    match name {
        "gyroid" => x.sin() * y.cos() + y.sin() * z.cos() + z.sin() * x.cos(),
        "schwarz_p" => x.cos() + y.cos() + z.cos(),
        _ => {
            x.sin() * y.sin() * z.sin()
                + x.sin() * y.cos() * z.cos()
                + x.cos() * y.sin() * z.cos()
                + x.cos() * y.cos() * z.sin()
        }
    }
}
fn native_gradient(name: &str, p: Point) -> Point {
    let [x, y, z] = p.map(|v| TAU * v);
    let [sx, sy, sz] = [x.sin(), y.sin(), z.sin()];
    let [cx, cy, cz] = [x.cos(), y.cos(), z.cos()];
    match name {
        "gyroid" => [cx * cy - sz * sx, cy * cz - sx * sy, cz * cx - sy * sz],
        "schwarz_p" => [-sx, -sy, -sz],
        _ => [
            cx * sy * sz + cx * cy * cz - sx * sy * cz - sx * cy * sz,
            sx * cy * sz - sx * sy * cz + cx * cy * cz - cx * sy * sz,
            sx * sy * cz - sx * cy * sz - cx * sy * sz + cx * cy * cz,
        ],
    }
    .map(|v| TAU * v)
}
impl Field {
    fn save(&self, error: PyErr) {
        let mut slot = self.error.lock().unwrap();
        if slot.is_none() {
            *slot = Some(error);
        }
    }
    fn evaluate(&self, points: &[Point]) -> Vec<f64> {
        if self.check().is_err() {
            return vec![f64::NAN; points.len()];
        }
        if let Some(name) = &self.native {
            return points.iter().map(|&p| native_value(name, p)).collect();
        }
        let mut values = Vec::with_capacity(points.len());
        for chunk in points.chunks(self.batch_size) {
            self.calls.fetch_add(1, Ordering::Relaxed);
            let result = Python::attach(|py| -> PyResult<Vec<f64>> {
                let x = chunk
                    .iter()
                    .map(|p| p[0])
                    .collect::<Vec<_>>()
                    .into_pyarray(py);
                let y = chunk
                    .iter()
                    .map(|p| p[1])
                    .collect::<Vec<_>>()
                    .into_pyarray(py);
                let z = chunk
                    .iter()
                    .map(|p| p[2])
                    .collect::<Vec<_>>()
                    .into_pyarray(py);
                let value = self.callback.as_ref().unwrap().call1(py, (x, y, z))?;
                let array = value
                    .extract::<PyReadonlyArray1<'_, f64>>(py)
                    .map_err(|_| {
                        PyTypeError::new_err("field must return a 1D float64 NumPy array")
                    })?;
                let data = array.as_array();
                if data.len() != chunk.len() {
                    return Err(PyValueError::new_err(
                        "field returned the wrong number of values",
                    ));
                }
                if data.iter().any(|v| !v.is_finite()) {
                    return Err(PyValueError::new_err("field returned nonfinite values"));
                }
                Ok(data.iter().copied().collect())
            });
            match result {
                Ok(v) => values.extend(v),
                Err(e) => {
                    self.save(e);
                    return vec![f64::NAN; points.len()];
                }
            }
        }
        values
    }
}
impl ScalarField for Field {
    fn value(&self, p: Point) -> f64 {
        if let Some((kernel, _)) = self.compiled {
            if self.check().is_err() {
                return f64::NAN;
            }
            return match compiled::evaluate(kernel, p) {
                Ok(out) => out[0],
                Err(()) => {
                    self.save(PyValueError::new_err("compiled field evaluation failed"));
                    f64::NAN
                }
            };
        }
        if let Some(name) = &self.native {
            if self.check().is_err() {
                return f64::NAN;
            }
            return native_value(name, p);
        }
        self.evaluate(&[p])[0]
    }
    fn values(&self, p: &[Point]) -> Vec<f64> {
        if self.compiled.is_some() {
            return p.iter().map(|&p| self.value(p)).collect();
        }
        self.evaluate(p)
    }
    fn prefers_batches(&self) -> bool {
        self.compiled.is_none() && (self.native.is_none() || self.derivative.is_some())
    }
    fn gradient(&self, p: Point) -> Option<Point> {
        if let Some((kernel, has_gradient)) = self.compiled {
            if !has_gradient {
                return None;
            }
            if self.check().is_err() {
                return Some([f64::NAN; 3]);
            }
            return Some(match compiled::evaluate(kernel, p) {
                Ok(out) => [out[1], out[2], out[3]],
                Err(()) => {
                    self.save(PyValueError::new_err("compiled gradient evaluation failed"));
                    [f64::NAN; 3]
                }
            });
        }
        if let Some(name) = self.native.as_ref().filter(|_| self.derivative.is_none()) {
            return Some(native_gradient(name, p));
        }
        self.gradients(&[p])[0]
    }
    fn gradients(&self, points: &[Point]) -> Vec<Option<Point>> {
        if let Some(name) = self.native.as_ref().filter(|_| self.derivative.is_none()) {
            return points
                .iter()
                .map(|&p| Some(native_gradient(name, p)))
                .collect();
        }
        let mut output = Vec::with_capacity(points.len());
        for chunk in points.chunks(self.batch_size.min(4096)) {
            if self.check().is_err() {
                return vec![Some([f64::NAN; 3]); points.len()];
            }
            if let Some(derivative) = &self.derivative {
                self.calls.fetch_add(1, Ordering::Relaxed);
                let result = Python::attach(|py| -> PyResult<Vec<Option<Point>>> {
                    let arrays: Vec<_> = (0..3)
                        .map(|a| {
                            chunk
                                .iter()
                                .map(|p| p[a])
                                .collect::<Vec<_>>()
                                .into_pyarray(py)
                        })
                        .collect();
                    let value = derivative.call1(py, (&arrays[0], &arrays[1], &arrays[2]))?;
                    let array = value
                        .extract::<PyReadonlyArray2<'_, f64>>(py)
                        .map_err(|_| {
                            PyTypeError::new_err(
                                "gradient must return a (N, 3) float64 NumPy array",
                            )
                        })?;
                    let a = array.as_array();
                    if a.shape() != [chunk.len(), 3] || a.iter().any(|v| !v.is_finite()) {
                        return Err(PyValueError::new_err("invalid gradient shape or values"));
                    }
                    Ok((0..chunk.len())
                        .map(|i| Some([a[[i, 0]], a[[i, 1]], a[[i, 2]]]))
                        .collect())
                });
                match result {
                    Ok(v) => output.extend(v),
                    Err(e) => {
                        self.save(e);
                        return vec![Some([f64::NAN; 3]); points.len()];
                    }
                }
            } else {
                let mut samples = Vec::with_capacity(chunk.len() * 6);
                for &p in chunk {
                    for a in 0..3 {
                        let mut lo = p;
                        let mut hi = p;
                        lo[a] -= self.step[a];
                        hi[a] += self.step[a];
                        samples.extend([lo, hi]);
                    }
                }
                let values = self.evaluate(&samples);
                output.extend(values.as_chunks::<6>().0.iter().map(|v| {
                    Some(std::array::from_fn(|a| {
                        (v[2 * a + 1] - v[2 * a]) / (2. * self.step[a])
                    }))
                }));
            }
        }
        output
    }
    fn check(&self) -> Result<(), String> {
        if self.cancellation.load(Ordering::Relaxed) {
            self.save(CancelledError::new_err("generation cancelled"));
            return Err("generation cancelled".into());
        }
        if self.error.lock().unwrap().is_some() {
            return Err("Python field failed".into());
        }
        if self
            .checks
            .fetch_add(1, Ordering::Relaxed)
            .is_multiple_of(256)
            && let Err(error) = Python::attach(|py| py.check_signals())
        {
            self.save(error);
            return Err("Python interrupted".into());
        }
        // Native field calls avoid taking the GIL per point.
        Ok(())
    }
}

fn make_field(
    py: Python<'_>,
    field: Py<PyAny>,
    gradient: Option<Py<PyAny>>,
    bounds: [Point; 2],
    cancellation: Arc<AtomicBool>,
    batch_size: usize,
) -> PyResult<Field> {
    if batch_size == 0 || batch_size > 65536 {
        return Err(PyValueError::new_err("batch_size must be 1..65536"));
    }
    let compiled = field
        .extract::<PyRef<'_, compiled::CompiledField>>(py)
        .ok()
        .map(|c| (c.kernel, c.has_gradient));
    if compiled.is_some() && gradient.is_some() {
        return Err(PyValueError::new_err(
            "supply the gradient to compile_field, not generate",
        ));
    }
    let native = field.extract::<String>(py).ok();
    if let Some(name) = &native {
        if !["gyroid", "schwarz_p", "schwarz_d"].contains(&name.as_str()) {
            return Err(PyValueError::new_err("unknown native field"));
        }
    } else if !field.bind(py).is_callable() {
        return Err(PyTypeError::new_err(
            "field must be callable or a native field name",
        ));
    }
    if gradient.as_ref().is_some_and(|g| !g.bind(py).is_callable()) {
        return Err(PyTypeError::new_err("gradient must be callable"));
    }
    Ok(Field {
        compiled,
        native,
        callback: Some(field),
        derivative: gradient,
        error: Mutex::new(None),
        cancellation,
        step: std::array::from_fn(|a| (bounds[1][a] - bounds[0][a]) * 1e-5),
        calls: AtomicUsize::new(0),
        checks: AtomicUsize::new(0),
        batch_size,
    })
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
fn generate<'py>(
    py: Python<'py>,
    field: Py<PyAny>,
    gradient: Option<Py<PyAny>>,
    bounds: [Point; 2],
    cells: [usize; 3],
    periodic: [bool; 3],
    band: Option<[f64; 2]>,
    tolerance: f64,
    minimum_quality: f64,
    budget: usize,
    passes: usize,
    snap: f64,
    threads: usize,
    batch_size: usize,
    cancel: Option<PyRef<'py, CancellationToken>>,
) -> PyResult<Bound<'py, PyDict>> {
    let f = make_field(
        py,
        field,
        gradient,
        bounds,
        cancel.as_ref().map(|v| v.flag.clone()).unwrap_or_default(),
        batch_size,
    )?;
    let options = implicit::Options {
        bounds,
        cells,
        periodic,
        region: band.map_or(implicit::Region::Negative, |[lower, upper]| {
            implicit::Region::Band { lower, upper }
        }),
        geometry_tolerance: tolerance,
        minimum_quality,
        max_tetrahedra: budget,
        optimize_passes: passes,
        snap,
        threads,
    };
    let start = std::time::Instant::now();
    let result = py.detach(|| implicit::generate(&f, options));
    if let Some(error) = f.error.lock().unwrap().take() {
        return Err(error);
    }
    py.check_signals()?;
    let result = result.map_err(|e| match e {
        MeshingError::InvalidOptions(message) => PyValueError::new_err(message),
        MeshingError::GenerationFailed(message) => MeshingFailure::new_err(message),
    })?;
    let output = PyDict::new(py);
    let points = Array2::from_shape_vec(
        (result.mesh.points.len(), 3),
        result.mesh.points.into_iter().flatten().collect(),
    )
    .unwrap();
    let cells = Array2::from_shape_vec(
        (result.mesh.tets.len(), 4),
        result
            .mesh
            .tets
            .into_iter()
            .flatten()
            .map(|i| i as i64)
            .collect(),
    )
    .unwrap();
    let surface = Array2::from_shape_vec(
        (result.mesh.surface.len(), 3),
        result
            .mesh
            .surface
            .into_iter()
            .flatten()
            .map(|i| i as i64)
            .collect(),
    )
    .unwrap();
    output.set_item("points", points.into_pyarray(py))?;
    output.set_item("tetrahedra", cells.into_pyarray(py))?;
    output.set_item("surface", surface.into_pyarray(py))?;
    output.set_item("boundary_tags", result.boundary_tags.into_pyarray(py))?;
    let pairs = pyo3::types::PyList::empty(py);
    for axis in result.periodic_pairs {
        pairs.append(
            Array2::from_shape_vec(
                (axis.len(), 2),
                axis.into_iter().flatten().map(|i| i as i64).collect(),
            )
            .unwrap()
            .into_pyarray(py),
        )?;
    }
    output.set_item("periodic_pairs", pairs)?;
    let diagnostics = PyDict::new(py);
    diagnostics.set_item("empty", result.diagnostics.empty)?;
    diagnostics.set_item("volume", result.diagnostics.volume)?;
    diagnostics.set_item(
        "minimum_mmg_quality",
        result.diagnostics.minimum_mmg_quality,
    )?;
    diagnostics.set_item(
        "sampled_surface_error",
        result.diagnostics.maximum_sampled_surface_error,
    )?;
    diagnostics.set_item(
        "elements_below_quality_01",
        result.diagnostics.elements_below_quality_01,
    )?;
    diagnostics.set_item(
        "reverted_snap_vertices",
        result.diagnostics.reverted_snap_vertices,
    )?;
    diagnostics.set_item("callback_calls", f.calls.load(Ordering::Relaxed))?;
    diagnostics.set_item("seconds", start.elapsed().as_secs_f64())?;
    output.set_item("diagnostics", diagnostics)?;
    Ok(output)
}

/// Return the available native CPU worker count.
#[pyfunction]
fn available_threads() -> usize {
    meshers_core::available_threads()
}

#[pymodule(gil_used = true)]
fn _meshers(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(available_threads, m)?)?;
    m.add_class::<compiled::CompiledField>()?;
    m.add_function(wrap_pyfunction!(compiled::_compile_expression, m)?)?;
    m.add_function(wrap_pyfunction!(generate, m)?)?;
    m.add_function(wrap_pyfunction!(intersection::generate_intersection, m)?)?;
    m.add_class::<CancellationToken>()?;
    m.add("MeshingError", m.py().get_type::<MeshingFailure>())?;
    m.add("CancelledError", m.py().get_type::<CancelledError>())?;
    Ok(())
}
