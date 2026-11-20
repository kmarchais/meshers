# Automatic native fields

Meshers traces a pointwise function, differentiates its expression graph and emits float64 machine code with Cranelift. CPU workers evaluate the compiled field with the GIL released.

```python
import numpy as np
import meshers


def field(x, y, z):
    return np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y) + np.sin(2 * np.pi * z)


mesh = meshers.generate(field, band=(-0.5, 0.5), cells=24)
```

- **Operations** — scalar arithmetic, comparisons, supported trigonometric and exponential functions, min/max and `np.where`. At a min/max tie, the derivative averages both gradients.
- **Reuse** — automatic preparation retraces captured parameters. A 32-entry cache reuses matching compiled expressions across resolutions. `compile_field` returns an immutable snapshot.
- **Callbacks** — use `compile_expression=False` when executable memory is unavailable. Unsupported expressions also fall back to batched NumPy with a warning. Return one value per point and keep values independent of batch size.

[Field conventions](https://kmarchais.github.io/meshers/fields/) · [Evaluator examples](https://github.com/kmarchais/meshers/blob/main/crates/meshers-python/examples/evaluators.py)
