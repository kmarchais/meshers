# Fields, thickness and grading

Define a deterministic pointwise function `f(x, y, z)`. Without `band`, meshers
keeps the negative side. With `band=(lower, upper)`, it keeps the interval in raw
field units. Levelset width is not automatically a constant physical thickness.

```python
import numpy as np
import meshers


def graded_sheet(x, y, z):
    a, b, c = 2 * np.pi * x, 2 * np.pi * y, 2 * np.pi * z
    gyroid = np.sin(a) * np.cos(b) + np.sin(b) * np.cos(c) + np.sin(c) * np.cos(a)
    schwarz = np.cos(a) + np.cos(b) + np.cos(c)
    weight = 0.2 + 0.3 * x
    half_width = 0.4 + 0.1 * z
    return ((1 - weight) * gyroid + weight * schwarz) / half_width


mesh = meshers.generate(graded_sheet, band=(-1, 1), cells=24, geometry_tolerance=0.05)
```

## Automatic preparation

The wheel traces arithmetic and supported NumPy operations into an expression
graph, differentiates the full graph and compiles it in memory with Cranelift.
Grading derivatives are included. No Hessian is needed by the current algorithm.

Supported operations include powers, sin/cos/tan, exp/log, sqrt, sinh/cosh/tanh,
abs, minimum/maximum, comparisons and `np.where`. Python branches on coordinates,
arbitrary array operations and arbitrary Python programs are not traceable.

Unsupported operations produce a warning and use batched NumPy callbacks.
`compile=False` selects that path explicitly. Callbacks must return float64 `(N,)`
arrays; supplied gradient callbacks return `(N,3)`. Avoid side effects and
operations that depend on batch length or values elsewhere in a batch.

The 32-entry cache includes captured constants, not mesh resolution. Each
preparation retraces to observe changed parameters. An explicit `compile_field`
call is optional and returns an immutable snapshot.

At abs(0), the derivative is zero; min/max ties use the first branch; `where`
uses the chosen branch. These conventions do not make nonsmooth or discontinuous
geometries smooth. Keep thickness positive throughout the domain.


## Coordinate maps and nonsmooth points

The embedded compiler supports `arctan2`, `arccos`, `arcsin`, `arctan` and
positional `clip` in addition to scalar arithmetic and trigonometric functions.
Geometry providers own their coordinate maps, units, envelopes and grading.
Meshers evaluates the resulting Cartesian field and meshes its negative region.

Automatic derivatives choose a zero subgradient for a vector norm at its origin.
This avoids `0/0` inside common distance fields; it does not make a singular
coordinate map smooth. Invalid domains and nonfinite derivatives still fail.
Full spherical maps need a unique inside/outside limit at their poles. A swept field must be
continuous along its curve and match across its angular seam.


An angular field can describe incompatible inside/outside regions arbitrarily
close to the same axis point. Increasing resolution cannot repair that geometry.
The geometry provider must define a regular field or let the caller explicitly
exclude the singular region. An excluded region may use a regular positive field
extension. The caller's negative region remains the meshing contract.
