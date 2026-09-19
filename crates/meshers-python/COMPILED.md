# Automatic native fields

Pass the usual NumPy-style field to `meshers.generate`. The wheel includes the
Cranelift compiler; users do not install Rust, a C compiler, Numba or SymPy.

```python
import numpy as np
import meshers


def field(x, y, z):
    return np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y) + np.sin(2 * np.pi * z)


mesh = meshers.generate(field, band=(-0.5, 0.5), cells=24)
```

Preparation traces scalar arithmetic and supported NumPy operations into a shared
expression graph, differentiates it and emits float64 machine code in memory.
No generated source, subprocess or external compiler runs. The mesher releases
the GIL and evaluates the compiled graph directly on its CPU workers.

## Supported fields

- Arithmetic, powers, sin/cos/tan, exp/log, sqrt, sinh/cosh/tanh, abs, min/max,
  comparisons and `np.where`.
- Captured scalar parameters, callable thickness and blends of TPMS functions.
- Automatic first derivatives, or a supplied gradient returning three expressions.
  The current algorithm does not need a Hessian.

Functions must be deterministic and pointwise. Tracing executes the function once;
do not use side effects or quantities that depend on array length or batch contents.
Unsupported array operations use NumPy callbacks with a warning. `compile=False`
selects callbacks explicitly. Arbitrary Python programs are not compilable.
At nonsmooth points, abs uses derivative zero; min/max use the first branch at a
tie; `where` uses the selected branch. Discontinuous fields or undefined derivatives
still need appropriate geometry and resolution checks.

The 32-entry in-memory cache keys the expression graph, including captured scalar
values. Each call retraces to observe changed parameters; mesh resolution does not
invalidate the compiled graph. `compile_field` remains available for explicit,
immutable snapshots but is unnecessary in normal use. `gradient=False` on that
optional API retains finite differences for comparison.

## Distribution

Install the wheel matching the operating system and architecture. Python 3.10+
uses the stable Python ABI. Wheels must be built by maintainers, with Rust, before
distribution; installing from a source archive still requires a build toolchain.
Use `pip install --only-binary=meshers <wheel-or-release>` to avoid a source build.
Install with `python -m pip install meshers`.

CI builds and tests Windows, Linux and macOS wheels, including Apple Silicon
and Intel Mac artifacts.
Embedded JIT compilation requires permission to allocate executable memory.
The supported NumPy fallback remains available where that is disallowed.

See the [Python examples](../../docs/site/examples.md) for compiled fields and callbacks.
