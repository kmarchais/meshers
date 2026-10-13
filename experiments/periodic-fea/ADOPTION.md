# Prototype direct periodic generation in Microgen

Start with the existing one-pass native generator. The tested Gyroid and Schwarz
P profiles generate admitted periodic meshes 7.49 and 6.35 times faster than the
standard native-seed MMG controls. Use [SUMMARY.json](SUMMARY.json) for the
measurements and [RESULTS.md](RESULTS.md) for reference sensitivity and failed
algorithm controls. Production Meshers has not changed during this investigation.

Use these resolutions for the studied unit-cell fixtures, not as defaults for
arbitrary geometries. Worst strain-energy differences stay below 1% against
grid 80 and two quadratic extrapolations. The extrapolations assume a leading
error proportional to 1/cells^2; they do not certify continuum accuracy.

## Generate the studied periodic Gyroid

Use raw band [-0.55, 0.55] in the unit cube. Normalize the field for the public
Meshers band [-1, 1], then construct Microgen's public `BoxMesh` from the arrays.
This uses no VTK geometry algorithm or MMG remeshing step.

```python
import numpy as np
import pyvista as pv
import meshers
from microgen import BoxMesh, is_periodic

def field(x, y, z):
    a, b, c = 2 * np.pi * x, 2 * np.pi * y, 2 * np.pi * z
    return (np.sin(a) * np.cos(b) + np.sin(b) * np.cos(c)
            + np.sin(c) * np.cos(a)) / 0.55

native = meshers.generate(
    field,
    bounds=(-0.5, 0.5) * 3,
    band=(-1, 1),
    cells=48,
    periodic=(True, True, True),
    optimize_passes=1,
    minimum_quality=0.1,
    geometry_tolerance=0.002,
    max_tetrahedra=2_000_000,
    threads=1,
)
box = BoxMesh(
    np.asarray(native.points),
    {pv.CellType.TETRA: np.asarray(native.tetrahedra, dtype=np.int64)},
)
assert is_periodic(box.nodes_coords, tol=1e-10)
```

For the studied Schwarz P fixture, use
`(cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z)) / 0.6` and `cells=52`.
Keep the other options fixed. These profiles use the verified packaged PR #5
build, source commit `66b800cf80e791788b980c01a8438f0e7f8ec8c6`; the installed
released package was not the candidate. Reuse the frozen native package as
described in [README.md](README.md).

## Verify the application boundary

`verify_microgen_bridge.py` checks actual admitted native artifacts through
Microgen 2.0.0b1's public `BoxMesh.from_pyvista` and `to_pyvista` methods. Both
roundtrips preserve every coordinate and tetrahedron exactly, and both pass
Microgen's public `is_periodic` at tolerance 1e-10. The original independent
cap-triangle and periodic-equivalence checks remain necessary; the public node
check alone does not check triangle conformity or FEA accuracy.

Run the bridge check on each adopted fixture:

```powershell
benchmarks/microgen/.venv/Scripts/python.exe experiments/periodic-fea/verify_microgen_bridge.py <native-artifact/result.json>
```

Keep the same geometry, volume, positive/manifold/nonoverlap, shape and periodic
admission checks when adding a production backend. Reject failed meshes rather
than reporting them as qualified FEA results. Use a fixture-specific convergence
study for new bands, thin features and other geometries. Printing and visualization
need a separate benchmark; these FEA profiles do not establish the fastest
geometry-only surface route.

Do not port the tested scalar residual map or uniform fixed-surface HXT control
on these results. The scalar map gave little benefit at similar element count
and required prior FEA. HXT was slower and failed quality or accuracy at
comparable complexity. Anisotropic adaptation and periodic restricted Delaunay
remain unmeasured alternatives, not demonstrated improvements.
