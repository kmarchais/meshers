# Meshers Python 0.1.0

In-process CPU float64 meshing of bounded implicit solids. This is the v0.1.0 release candidate. No subprocess or intermediate mesh file is needed.

## Install and use

Install a prebuilt wheel matching your platform. Users need Python, not Rust,
Numba or SymPy. The candidate has not yet been published on PyPI.

```sh
python -m pip install /path/to/meshers.whl h5py
```

```python
import meshers

mesh = meshers.generate("gyroid", cells=24, band=(-0.5, 0.5), periodic=(True, True, True))
points, tetrahedra = mesh.points, mesh.tetrahedra
mesh.write_vtkhdf("gyroid.vtkhdf")  # requires h5py
```

Native fields: `gyroid`, `schwarz_p`, `schwarz_d`; their period is one physical
coordinate unit. Bounds default to the unit cube and use
`(xmin, xmax, ymin, ymax, zmin, zmax)`. `cells` accepts an integer or three integers.
Without `band`, the solid is `field <= 0`. Periodicity defaults to disabled.

## Custom and graded functions

Pass a pointwise NumPy-style callable `field(x, y, z)`. Supported expressions
compile automatically inside the wheel, with automatic gradients and cached
machine code. CPU workers evaluate float64 directly, with the GIL released.
This supports callable grading, thickness and blends without an explicit
compilation step. See [supported operations and limits](https://github.com/kmarchais/meshers/blob/main/crates/meshers-python/COMPILED.md).

Unsupported expressions fall back to batched NumPy callbacks with a warning.
Select that path explicitly with `compile=False`; callbacks must return `(N,)`
float64 arrays, and optional gradients `(N,3)`. Functions must be deterministic
and independent of batch contents. Callback evaluation runs on the calling
Python thread.

See the [graded example](https://github.com/kmarchais/meshers/blob/main/crates/meshers-python/examples/graded.py).

## Results and failures

- Owned NumPy arrays: float64 points, int64 tetrahedra and boundary triangles.
- Boundary tags: 0 implicit, 1/2 x-/x+, 3/4 y-/y+, 5/6 z-/z+.
- Per-axis periodic node pairs plus geometry and quality diagnostics.
- VTKHDF volume or `surface=True`, periodic master IDs and physical shift vectors.
  Volume files include `Volume` and `MMGQuality`, computed using the MMG formula
  in NumPy, not by calling mmgpy.
- Invalid options raise `ValueError`; failed meshing raises `MeshingError`.
  Callback exceptions are preserved. `KeyboardInterrupt` is propagated.
- Pass `cancel=CancellationToken()` and call `token.cancel()` from another thread
  for cooperative cancellation. A callback must return before it can be interrupted
  by that token. Tokens remain cancelled after use.

Output arrays own their storage and are mutable. Mutating them can invalidate
connectivity, metadata or diagnostics. Surface output is the closed solid boundary,
not an independent zero-level surface mesher. Small features may be missed at
insufficient resolution; successful sampling is not geometry certification.

## Verification and remaining work

The test suite covers automatic derivatives, cache invalidation, cancellation,
periodicity, bounded callbacks, named intersections, physical coordinate scales
and VTKHDF round trips. Linux and Windows source and wheel checks have passed
locally. Final macOS and portable artifact checks remain required before
publication. See the [release validation](https://github.com/kmarchais/meshers/blob/main/docs/site/release.md)
for the tested source, platform results and remaining gates. Thin, nonsmooth or
poorly resolved geometry still needs application-specific convergence checks.

See the [runnable examples](https://github.com/kmarchais/meshers/blob/main/docs/site/examples.md) and [performance regression checks](https://github.com/kmarchais/meshers/blob/main/docs/site/performance.md).

## Maintainer build

Maintainers need a compatible Rust toolchain and Python to build wheels:

```sh
python -m pip install maturin
maturin build --release --interpreter python --out dist
```

End users install the resulting wheels. Source installations require Rust.

## Retained intersections and mapped periodic ends

`meshers.generate_intersection` accepts named negative sublevel constraints,
an optional coordinate map and explicit rigid periodic transforms. It returns
the same Mesh type, with named surface access, feature edges and parameter
coordinates. See the [intersection guide](https://github.com/kmarchais/meshers/blob/main/docs/site/intersections.md)
and `examples/curved_tpms.py`. Source builds require Rust 1.95 or newer.
