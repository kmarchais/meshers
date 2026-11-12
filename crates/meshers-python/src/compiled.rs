//! Owned native expression compiler embedded in the wheel. No external compiler.
use cranelift_codegen::ir::{
    AbiParam, InstBuilder, MemFlagsData, Value, condcodes::FloatCC, types,
};
use cranelift_frontend::{FunctionBuilder, FunctionBuilderContext};
use cranelift_jit::{JITBuilder, JITModule};
use cranelift_module::{Linkage, Module, default_libcall_names};
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray2};
use pyo3::{
    exceptions::{PyRuntimeError, PyValueError},
    prelude::*,
};
use std::collections::HashMap;
use std::sync::Mutex;
mod polyhedral;

pub(crate) type Kernel = unsafe extern "C" fn(f64, f64, f64, *mut f64) -> i32;

// The finalized module is never modified or accessed concurrently. Only its machine
// code is executed, through CompiledField, whose Python owner outlives all worker calls.
// The module's non-Send builder closures are unused after finalization. Drop unmaps
// memory only after the last owner is gone, on whichever thread releases that owner.
struct CodeMemory(Option<JITModule>);
unsafe impl Send for CodeMemory {}
unsafe impl Sync for CodeMemory {}
impl Drop for CodeMemory {
    fn drop(&mut self) {
        if let Some(module) = self.0.take() {
            // SAFETY: no kernel pointer is publicly exposed; every call retains owner.
            unsafe { module.free_memory() };
        }
    }
}

#[pyclass(frozen)]
pub(crate) struct CompiledField {
    pub(crate) kernel: Kernel,
    pub(crate) value_kernel: Kernel,
    #[pyo3(get)]
    pub(crate) has_gradient: bool,
    #[pyo3(get)]
    node_count: usize,
    #[pyo3(get)]
    _halfspaces: Option<Vec<[f64; 4]>>,
    #[pyo3(get)]
    _polyhedral_regions: Option<polyhedral::Regions>,
    #[pyo3(get)]
    _sdf_intersection: bool,
    constraint_data: Option<(Vec<Node>, Vec<Vec<usize>>)>,
    constraint_cache: Mutex<Option<Vec<Py<CompiledField>>>>,
    _memory: CodeMemory,
}

pub(crate) fn evaluate(kernel: Kernel, p: [f64; 3]) -> Result<[f64; 4], ()> {
    let mut out = [f64::NAN; 4];
    // SAFETY: validated graph compilation emits this exact ABI, writing <=4 doubles;
    // CompiledField remains owned by the caller throughout detached meshing.
    let status = unsafe { kernel(p[0], p[1], p[2], out.as_mut_ptr()) };
    if status == 0 { Ok(out) } else { Err(()) }
}
#[pymethods]
impl CompiledField {
    /// Return cached evaluators for validated zero-sublevel intersection branches.
    fn _intersection_fields(&self, py: Python<'_>) -> PyResult<Option<Vec<Py<CompiledField>>>> {
        let Some((nodes, outputs)) = &self.constraint_data else {
            return Ok(None);
        };
        let mut cached = self
            .constraint_cache
            .lock()
            .map_err(|_| PyRuntimeError::new_err("retained constraint cache is unavailable"))?;
        if cached.is_none() {
            let fields = outputs
                .iter()
                .map(|outputs| {
                    Py::new(
                        py,
                        _compile_expression(nodes.clone(), outputs.clone(), None, false)?,
                    )
                })
                .collect::<PyResult<Vec<_>>>()?;
            *cached = Some(fields);
        }
        Ok(cached
            .as_ref()
            .map(|fields| fields.iter().map(|f| f.clone_ref(py)).collect()))
    }
    /// Evaluate owned point batches without a Python call per point.
    fn values<'py>(
        &self,
        py: Python<'py>,
        points: PyReadonlyArray2<'py, f64>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let array = points.as_array();
        if array.ncols() != 3 {
            return Err(PyValueError::new_err("points must have shape (n, 3)"));
        }
        // Copy before releasing the GIL: another Python thread may modify the
        // source array, and strided arrays need not have contiguous storage.
        let owned: Vec<[f64; 3]> = array
            .rows()
            .into_iter()
            .map(|row| [row[0], row[1], row[2]])
            .collect();
        let mut result = Vec::with_capacity(owned.len());
        for chunk in owned.chunks(4096) {
            py.check_signals()?;
            let values = py.detach(|| {
                chunk
                    .iter()
                    .map(|&point| {
                        let out = evaluate(self.value_kernel, point)
                            .map_err(|_| PyRuntimeError::new_err("compiled field failed"))?;
                        if !out[0].is_finite() {
                            return Err(PyValueError::new_err("nonfinite compiled field"));
                        }
                        Ok(out[0])
                    })
                    .collect::<PyResult<Vec<f64>>>()
            })?;
            result.extend(values);
        }
        py.check_signals()?;
        Ok(result.into_pyarray(py))
    }
    fn __call__(&self, x: f64, y: f64, z: f64) -> PyResult<f64> {
        let out = evaluate(self.value_kernel, [x, y, z])
            .map_err(|_| PyRuntimeError::new_err("compiled field failed"))?;
        if !out[0].is_finite() {
            return Err(PyValueError::new_err("nonfinite compiled field"));
        }
        Ok(out[0])
    }
    fn gradient(&self, x: f64, y: f64, z: f64) -> PyResult<[f64; 3]> {
        if !self.has_gradient {
            return Err(PyValueError::new_err("field has no compiled gradient"));
        }
        let out = evaluate(self.kernel, [x, y, z])
            .map_err(|_| PyRuntimeError::new_err("compiled field failed"))?;
        let g = [out[1], out[2], out[3]];
        if g.iter().any(|v| !v.is_finite()) {
            return Err(PyValueError::new_err("nonfinite compiled gradient"));
        }
        Ok(g)
    }
}

macro_rules! unary_math {
    ($($name:ident => $method:ident),*) => {$(extern "C" fn $name(x:f64)->f64{x.$method()})*};
}
unary_math!(asin=>asin,acos=>acos,atan=>atan,sin=>sin,cos=>cos,tan=>tan,exp=>exp,ln=>ln,sinh=>sinh,cosh=>cosh,tanh=>tanh);
extern "C" fn atan2(a: f64, b: f64) -> f64 {
    a.atan2(b)
}
extern "C" fn pow(a: f64, b: f64) -> f64 {
    a.powf(b)
}
extern "C" fn minimum(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.min(b)
    }
}
extern "C" fn maximum(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

type Node = (String, Vec<usize>, f64);
fn emit_kernel(
    module: &mut JITModule,
    imports: &HashMap<&str, cranelift_module::FuncId>,
    nodes: &[Node],
    outputs: &[usize],
    name: &str,
) -> PyResult<cranelift_module::FuncId> {
    let err =
        |e: &dyn std::fmt::Display| PyRuntimeError::new_err(format!("embedded compiler: {e}"));
    // External math calls are not eliminated as dead code by the JIT. Prune
    // nodes unreachable from this kernel's outputs before lowering the graph.
    let mut required = vec![false; nodes.len()];
    for &output in outputs {
        required[output] = true;
    }
    for i in (0..nodes.len()).rev() {
        if required[i] && nodes[i].0 != "input" {
            for &argument in &nodes[i].1 {
                required[argument] = true;
            }
        }
    }
    let mut indices = vec![usize::MAX; nodes.len()];
    let mut pruned = Vec::new();
    for (i, (op, args, constant)) in nodes.iter().enumerate() {
        if required[i] {
            indices[i] = pruned.len();
            let arguments = if op == "input" {
                args.clone()
            } else {
                args.iter().map(|&arg| indices[arg]).collect()
            };
            pruned.push((op.clone(), arguments, *constant));
        }
    }
    let outputs: Vec<_> = outputs.iter().map(|&output| indices[output]).collect();
    let nodes = pruned;
    let mut ctx = module.make_context();
    ctx.func
        .signature
        .params
        .extend((0..3).map(|_| AbiParam::new(types::F64)));
    ctx.func
        .signature
        .params
        .push(AbiParam::new(module.target_config().pointer_type()));
    ctx.func.signature.returns.push(AbiParam::new(types::I32));
    let id = module
        .declare_function(name, Linkage::Local, &ctx.func.signature)
        .map_err(|e| err(&e))?;
    let mut fc = FunctionBuilderContext::new();
    {
        let mut b = FunctionBuilder::new(&mut ctx.func, &mut fc);
        let block = b.create_block();
        b.append_block_params_for_function_params(block);
        b.switch_to_block(block);
        let params = b.block_params(block).to_vec();
        let mut values: Vec<Value> = Vec::with_capacity(nodes.len());
        for (op, args, c) in &nodes {
            let v = match op.as_str() {
                "input" => params[args[0]],
                "constant" => b.ins().f64const(*c),
                "add" => b.ins().fadd(values[args[0]], values[args[1]]),
                "sub" => b.ins().fsub(values[args[0]], values[args[1]]),
                "mul" => b.ins().fmul(values[args[0]], values[args[1]]),
                "div" => b.ins().fdiv(values[args[0]], values[args[1]]),
                "neg" => b.ins().fneg(values[args[0]]),
                "abs" => b.ins().fabs(values[args[0]]),
                "sqrt" => b.ins().sqrt(values[args[0]]),
                // Keep cheap scalar operations in the native instruction stream.
                // Nonfinite min/max results retain NumPy's NaN propagation.
                "pow" if nodes[args[1]].0 == "constant" && nodes[args[1]].2 == 2. => {
                    b.ins().fmul(values[args[0]], values[args[0]])
                }
                "minimum" => b.ins().fmin(values[args[0]], values[args[1]]),
                "maximum" => b.ins().fmax(values[args[0]], values[args[1]]),
                "select" => {
                    let zero = b.ins().f64const(0.);
                    let cond = b.ins().fcmp(FloatCC::NotEqual, values[args[0]], zero);
                    b.ins().select(cond, values[args[1]], values[args[2]])
                }
                "lt" | "le" | "gt" | "ge" | "eq" | "ne" => {
                    let cc = match op.as_str() {
                        "lt" => FloatCC::LessThan,
                        "le" => FloatCC::LessThanOrEqual,
                        "gt" => FloatCC::GreaterThan,
                        "ge" => FloatCC::GreaterThanOrEqual,
                        "eq" => FloatCC::Equal,
                        _ => FloatCC::NotEqual,
                    };
                    let cond = b.ins().fcmp(cc, values[args[0]], values[args[1]]);
                    let one = b.ins().f64const(1.);
                    let zero = b.ins().f64const(0.);
                    b.ins().select(cond, one, zero)
                }
                name => {
                    let callee = module.declare_func_in_func(imports[name], b.func);
                    let argv: Vec<_> = args.iter().map(|&a| values[a]).collect();
                    let call = b.ins().call(callee, &argv);
                    b.inst_results(call)[0]
                }
            };
            values.push(v);
        }
        for (offset, &node) in outputs.iter().enumerate() {
            b.ins().store(
                MemFlagsData::new(),
                values[node],
                params[3],
                (offset * 8) as i32,
            );
        }
        let status = b.ins().iconst(types::I32, 0);
        b.ins().return_(&[status]);
        b.seal_all_blocks();
        b.finalize(module.target_config());
    }
    module.define_function(id, &mut ctx).map_err(|e| err(&e))?;
    Ok(id)
}

/// Recognize exactly sqrt(sum(max(f_i, 0)^2)) + min(max(f_i), 0).
/// Restrict branches to unit-distance convex affine or rigid radial forms;
/// arbitrary scaling/nonlinear branches keep scalar behavior. In real arithmetic
/// the zero sublevel set is f_i <= 0. Scalar/gradient kernels remain unchanged.
fn sdf_intersection_roots(nodes: &[Node], output: usize) -> Option<Vec<usize>> {
    fn leaves(nodes: &[Node], root: usize, operation: &str) -> Vec<usize> {
        let mut pending = vec![root];
        let mut roots = Vec::new();
        let mut seen = vec![false; nodes.len()];
        while let Some(index) = pending.pop() {
            if seen[index] {
                continue;
            }
            seen[index] = true;
            if nodes[index].0 == operation {
                pending.extend(nodes[index].1.iter().rev());
            } else {
                roots.push(index);
            }
        }
        roots
    }
    fn squared_base(nodes: &[Node], index: usize) -> Option<usize> {
        let (op, args, _) = &nodes[index];
        match op.as_str() {
            "pow" if nodes[args[1]].0 == "constant" && nodes[args[1]].2 == 2. => Some(args[0]),
            "mul" if args[0] == args[1] => Some(args[0]),
            _ => None,
        }
    }
    fn unit_branch(nodes: &[Node], root: usize) -> bool {
        let unit = |plane: &[f64; 4]| {
            let squared: f64 = plane[..3].iter().map(|v| v * v).sum();
            (squared - 1.).abs() <= 1e-12
        };
        if let Some(forms) = polyhedral::convex_forms(nodes, root) {
            return !forms.is_empty() && forms.iter().all(unit);
        }
        let (op, args, _) = &nodes[root];
        if op != "sub"
            || nodes[args[1]].0 != "constant"
            || nodes[args[1]].2 <= 0.
            || nodes[args[0]].0 != "sqrt"
        {
            return false;
        }
        // A radial distance in a rigid affine frame. Count every squared
        // term: deduplicating a repeated basis vector would hide a scale.
        let mut pending = vec![nodes[args[0]].1[0]];
        let mut forms = Vec::<[f64; 4]>::new();
        while let Some(index) = pending.pop() {
            if nodes[index].0 == "add" {
                pending.extend(&nodes[index].1);
            } else {
                let Some(base) = squared_base(nodes, index) else {
                    return false;
                };
                let Some(affine) = polyhedral::convex_forms(nodes, base) else {
                    return false;
                };
                if affine.len() != 1 || !unit(&affine[0]) {
                    return false;
                }
                if forms
                    .iter()
                    .any(|p| (0..3).map(|a| p[a] * affine[0][a]).sum::<f64>().abs() > 1e-12)
                {
                    return false;
                }
                forms.push(affine[0]);
            }
            if pending.len() + forms.len() > 3 {
                return false;
            }
        }
        (2..=3).contains(&forms.len())
    }
    let zero = |index: usize| nodes[index].0 == "constant" && nodes[index].2 == 0.;
    let unclamp = |index: usize, operation: &str| {
        let (op, args, _) = &nodes[index];
        if op != operation {
            return None;
        }
        if zero(args[0]) {
            Some(args[1])
        } else if zero(args[1]) {
            Some(args[0])
        } else {
            None
        }
    };
    let (op, args, _) = &nodes[output];
    if op != "add" {
        return None;
    }
    for (outside, inside) in [(args[0], args[1]), (args[1], args[0])] {
        if nodes[outside].0 != "sqrt" {
            continue;
        }
        let Some(inner) = unclamp(inside, "minimum") else {
            continue;
        };
        let roots = leaves(nodes, inner, "maximum");
        if !(2..=meshers_core::implicit::intersection::MAX_CONSTRAINTS).contains(&roots.len())
            || roots.iter().any(|&i| !unit_branch(nodes, i))
        {
            continue;
        }
        let mut squared = Vec::new();
        let mut complete = true;
        for term in leaves(nodes, nodes[outside].1[0], "add") {
            let Some(base) = squared_base(nodes, term) else {
                complete = false;
                break;
            };
            let Some(root) = unclamp(base, "maximum") else {
                complete = false;
                break;
            };
            squared.push(root);
        }
        squared.sort_unstable();
        squared.dedup();
        let mut expected = roots.clone();
        expected.sort_unstable();
        if complete && squared == expected {
            return Some(roots);
        }
    }
    None
}

#[pyfunction(signature = (nodes, outputs, constraint_outputs=None, *, retain_sdf=false))]
pub(crate) fn _compile_expression(
    nodes: Vec<Node>,
    outputs: Vec<usize>,
    constraint_outputs: Option<Vec<Vec<usize>>>,
    retain_sdf: bool,
) -> PyResult<CompiledField> {
    if nodes.is_empty()
        || nodes.len() > 20000
        || ![1, 4].contains(&outputs.len())
        || outputs.iter().any(|&i| i >= nodes.len())
    {
        return Err(PyValueError::new_err(
            "invalid expression dimensions (maximum 20000 nodes)",
        ));
    }
    for (i, (op, args, c)) in nodes.iter().enumerate() {
        let arity = match op.as_str() {
            "constant" => 0,
            "input" | "neg" | "abs" | "sqrt" | "sin" | "cos" | "tan" | "exp" | "ln" | "sinh"
            | "cosh" | "tanh" | "asin" | "acos" | "atan" => 1,
            "add" | "sub" | "mul" | "div" | "pow" | "minimum" | "maximum" | "lt" | "le" | "gt"
            | "ge" | "eq" | "ne" | "atan2" => 2,
            "select" => 3,
            _ => return Err(PyValueError::new_err("unknown expression operation")),
        };
        if args.len() != arity
            || !c.is_finite()
            || (op == "input" && args[0] > 2)
            || (op != "input" && args.iter().any(|&a| a >= i))
        {
            return Err(PyValueError::new_err("invalid expression operands"));
        }
    }
    if let Some(branches) = &constraint_outputs {
        let mut roots = Vec::new();
        let mut pending = vec![outputs[0]];
        let mut seen = vec![false; nodes.len()];
        while let Some(i) = pending.pop() {
            if seen[i] {
                continue;
            }
            seen[i] = true;
            if nodes[i].0 == "maximum" {
                pending.push(nodes[i].1[1]);
                pending.push(nodes[i].1[0]);
            } else {
                roots.push(i);
            }
        }
        if !(2..=meshers_core::implicit::intersection::MAX_CONSTRAINTS).contains(&roots.len())
            || branches.len() != roots.len()
            || branches.iter().zip(&roots).any(|(branch, &root)| {
                ![1, 4].contains(&branch.len())
                    || branch.first() != Some(&root)
                    || branch.iter().any(|&i| i >= nodes.len())
                    || nodes[root].0 == "constant"
            })
        {
            return Err(PyValueError::new_err(
                "retained constraints must match the outer maximum branches",
            ));
        }
    }
    let sdf_constraints = if retain_sdf && constraint_outputs.is_none() {
        sdf_intersection_roots(&nodes, outputs[0])
            .map(|roots| roots.into_iter().map(|root| vec![root]).collect())
    } else {
        None
    };
    let sdf_intersection = sdf_constraints.is_some();
    let constraint_outputs = constraint_outputs.or(sdf_constraints);
    let err =
        |e: &dyn std::fmt::Display| PyRuntimeError::new_err(format!("embedded compiler: {e}"));
    let mut builder = JITBuilder::with_flags(&[("opt_level", "speed")], default_libcall_names())
        .map_err(|e| err(&e))?;
    let math: [(&str, *const u8, usize); 15] = [
        ("atan2", atan2 as *const u8, 2),
        ("asin", asin as *const u8, 1),
        ("acos", acos as *const u8, 1),
        ("atan", atan as *const u8, 1),
        ("sin", sin as *const u8, 1),
        ("cos", cos as *const u8, 1),
        ("tan", tan as *const u8, 1),
        ("exp", exp as *const u8, 1),
        ("ln", ln as *const u8, 1),
        ("sinh", sinh as *const u8, 1),
        ("cosh", cosh as *const u8, 1),
        ("tanh", tanh as *const u8, 1),
        ("pow", pow as *const u8, 2),
        ("minimum", minimum as *const u8, 2),
        ("maximum", maximum as *const u8, 2),
    ];
    for &(name, address, _) in &math {
        builder.symbol(name, address);
    }
    let mut memory = CodeMemory(Some(JITModule::new(builder)));
    let module = memory.0.as_mut().expect("module initialized");
    let mut imports = HashMap::new();
    for &(name, _, arity) in &math {
        let mut sig = module.make_signature();
        sig.params
            .extend((0..arity).map(|_| AbiParam::new(types::F64)));
        sig.returns.push(AbiParam::new(types::F64));
        let id = module
            .declare_function(name, Linkage::Import, &sig)
            .map_err(|e| err(&e))?;
        imports.insert(name, id);
    }
    let id = emit_kernel(module, &imports, &nodes, &outputs, "field")?;
    let value_id = if outputs.len() == 4 {
        emit_kernel(module, &imports, &nodes, &[outputs[0]], "field_value")?
    } else {
        id
    };
    module.finalize_definitions().map_err(|e| err(&e))?;
    let address = module.get_finalized_function(id);
    // SAFETY: the function signature and stores above exactly match Kernel.
    let kernel = unsafe { std::mem::transmute::<*const u8, Kernel>(address) };
    let value_address = module.get_finalized_function(value_id);
    // SAFETY: emit_kernel uses the same ABI for both output sets.
    let value_kernel = unsafe { std::mem::transmute::<*const u8, Kernel>(value_address) };
    let halfspaces = polyhedral::halfspaces(&nodes, outputs[0]);
    let regions = if halfspaces.is_none() {
        polyhedral::regions(&nodes, outputs[0])
    } else {
        None
    };
    // Existing exact affine routes take precedence; do not retain an unused
    // second representation or compile the parent evaluator a second time.
    let constraint_outputs =
        constraint_outputs.filter(|_| halfspaces.is_none() && regions.is_none());
    Ok(CompiledField {
        kernel,
        value_kernel,
        has_gradient: outputs.len() == 4,
        node_count: nodes.len(),
        _halfspaces: halfspaces,
        _polyhedral_regions: regions,
        _sdf_intersection: sdf_intersection && constraint_outputs.is_some(),
        constraint_data: constraint_outputs.map(|outputs| (nodes, outputs)),
        constraint_cache: Mutex::new(None),
        _memory: memory,
    })
}
