# Python examples

These scripts use the installed public API. Run them from any writable directory;
exports are written there. Install h5py for the examples that write VTKHDF.
From a development checkout, use `uv run python` in place of `python`.

```sh
python crates/meshers-python/examples/basic.py
```

| Script | Features |
| --- | --- |
| [basic.py](https://github.com/kmarchais/meshers/blob/main/crates/meshers-python/examples/basic.py) | All three built-in TPMS fields, a custom negative-side solid, bands, bounds, per-axis resolution, arrays, volume and surface export |
| [graded.py](https://github.com/kmarchais/meshers/blob/main/crates/meshers-python/examples/graded.py) | Spatially varying thickness and automatic derivatives |
| [evaluators.py](https://github.com/kmarchais/meshers/blob/main/crates/meshers-python/examples/evaluators.py) | Automatic compilation, explicit reusable compilation, supplied gradients, batched callbacks, finite differences |
| [periodic.py](https://github.com/kmarchais/meshers/blob/main/crates/meshers-python/examples/periodic.py) | Three-axis node correspondence, quality floor, element budget, optimization, snapping and worker count |
| [intersection.py](https://github.com/kmarchais/meshers/blob/main/crates/meshers-python/examples/intersection.py) | Named surfaces, feature edges, parameter coordinates, straight periodicity, curved maps and rigid periodic transforms |
| [cancellation.py](https://github.com/kmarchais/meshers/blob/main/crates/meshers-python/examples/cancellation.py) | Cancellation tokens, invalid input and meshing failure handling |

CI runs all six scripts unchanged in temporary directories and checks exported
coordinates, connectivity and positive volumes. API regression tests cover
additional option combinations and failure cases. Every public option is
listed in the [API reference](api/python.md).

For a larger mapped TPMS example:

```sh
python crates/meshers-python/examples/curved_tpms.py torus torus.vtkhdf
python crates/meshers-python/examples/curved_tpms.py cylinder cylinder.vtkhdf
```

The simple scripts are intended as starting points. Increase resolution and
choose tolerances for your geometry; their settings are not universal defaults.
