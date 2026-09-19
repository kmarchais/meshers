//! Owned native expression compiler embedded in the wheel. No external compiler.
use cranelift_codegen::ir::{
    AbiParam, InstBuilder, MemFlagsData, Value, condcodes::FloatCC, types,
};
use cranelift_frontend::{FunctionBuilder, FunctionBuilderContext};
use cranelift_jit::{JITBuilder, JITModule};
use cranelift_module::{Linkage, Module, default_libcall_names};
use pyo3::{
    exceptions::{PyRuntimeError, PyValueError},
    prelude::*,
};
use std::collections::HashMap;

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
    #[pyo3(get)]
    pub(crate) has_gradient: bool,
    #[pyo3(get)]
    node_count: usize,
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
    fn __call__(&self, x: f64, y: f64, z: f64) -> PyResult<f64> {
        let out = evaluate(self.kernel, [x, y, z])
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
#[pyfunction]
pub(crate) fn _compile_expression(
    nodes: Vec<Node>,
    outputs: Vec<usize>,
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
        .declare_function("field", Linkage::Local, &ctx.func.signature)
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
    module.finalize_definitions().map_err(|e| err(&e))?;
    let address = module.get_finalized_function(id);
    // SAFETY: the function signature and stores above exactly match Kernel.
    let kernel = unsafe { std::mem::transmute::<*const u8, Kernel>(address) };
    Ok(CompiledField {
        kernel,
        has_gradient: outputs.len() == 4,
        node_count: nodes.len(),
        _memory: memory,
    })
}
