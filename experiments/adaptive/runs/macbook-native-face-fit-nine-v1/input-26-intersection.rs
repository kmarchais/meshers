//! Python handoff to the constrained native mesher; Python assembles named metadata.
use super::*;
use meshers_core::implicit::intersection;

#[allow(clippy::too_many_arguments)]
#[pyfunction]
pub(super) fn generate_intersection<'py>(
    py: Python<'py>,
    fields: Vec<Py<PyAny>>,
    mapping: Option<Vec<Py<PyAny>>>,
    bounds: [Point; 2],
    cells: [usize; 3],
    periodic: [bool; 3],
    minimum_quality: f64,
    budget: usize,
    passes: usize,
    snap: f64,
    batch_size: usize,
    cancel: Option<PyRef<'py, CancellationToken>>,
) -> PyResult<Bound<'py, PyDict>> {
    let flag = cancel.as_ref().map(|v| v.flag.clone()).unwrap_or_default();
    let fields = fields
        .into_iter()
        .map(|f| make_field(py, f, None, bounds, flag.clone(), batch_size))
        .collect::<PyResult<Vec<_>>>()?;
    let mapping = mapping
        .map(|fs| {
            fs.into_iter()
                .map(|f| make_field(py, f, None, bounds, flag.clone(), batch_size))
                .collect::<PyResult<Vec<_>>>()
        })
        .transpose()?;
    if mapping.as_ref().is_some_and(|fs| fs.len() != 3) {
        return Err(PyValueError::new_err(
            "coordinate map needs three components",
        ));
    }
    let refs: Vec<&dyn ScalarField> = fields.iter().map(|f| f as &dyn ScalarField).collect();
    let map = |p: Point| {
        mapping
            .as_ref()
            .map_or(p, |fs| std::array::from_fn(|a| fs[a].value(p)))
    };
    let start = std::time::Instant::now();
    let result = py.detach(|| {
        intersection::generate(
            &refs,
            map,
            implicit::Options {
                bounds,
                cells,
                periodic,
                snap,
                minimum_quality,
                max_tetrahedra: budget,
                optimize_passes: passes,
                geometry_tolerance: 0.,
                ..Default::default()
            },
        )
    });
    for f in fields.iter().chain(mapping.iter().flatten()) {
        if let Some(error) = f.error.lock().unwrap().take() {
            return Err(error);
        }
    }
    py.check_signals()?;
    let result = result.map_err(MeshingFailure::new_err)?;
    let out = PyDict::new(py);
    for (name, values) in [
        ("points", result.mesh.points),
        ("parameters", result.parameters),
    ] {
        out.set_item(
            name,
            Array2::from_shape_vec((values.len(), 3), values.into_iter().flatten().collect())
                .unwrap()
                .into_pyarray(py),
        )?;
    }
    out.set_item(
        "tetrahedra",
        Array2::from_shape_vec(
            (result.mesh.tets.len(), 4),
            result
                .mesh
                .tets
                .into_iter()
                .flatten()
                .map(|v| v as i64)
                .collect(),
        )
        .unwrap()
        .into_pyarray(py),
    )?;
    out.set_item(
        "surface",
        Array2::from_shape_vec(
            (result.mesh.surface.len(), 3),
            result
                .mesh
                .surface
                .into_iter()
                .flatten()
                .map(|v| v as i64)
                .collect(),
        )
        .unwrap()
        .into_pyarray(py),
    )?;
    out.set_item("constraint_masks", result.constraints.into_pyarray(py))?;
    let d = PyDict::new(py);
    d.set_item(
        "reverted_snap_vertices",
        result.quality.reverted_snap_vertices,
    )?;
    d.set_item("minimum_mmg_quality", result.quality.minimum_quality)?;
    d.set_item(
        "elements_below_quality_01",
        result.quality.elements_below_01,
    )?;
    d.set_item(
        "initial_minimum_mmg_quality",
        result.quality.initial_minimum_quality,
    )?;
    d.set_item(
        "initial_elements_below_quality_01",
        result.quality.initial_elements_below_01,
    )?;
    d.set_item("collapsed_vertices", result.quality.collapsed_vertices)?;
    d.set_item(
        "accepted_vertex_moves",
        result.quality.accepted_vertex_moves,
    )?;
    d.set_item(
        "accepted_reconnections",
        result.quality.accepted_reconnections,
    )?;
    d.set_item("seconds", start.elapsed().as_secs_f64())?;
    d.set_item("threads", 1)?;
    d.set_item(
        "callback_calls",
        fields
            .iter()
            .chain(mapping.iter().flatten())
            .map(|f| f.calls.load(Ordering::Relaxed))
            .sum::<usize>(),
    )?;
    out.set_item("diagnostics", d)?;
    Ok(out)
}
