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
These targets concern homogenized stiffness. The [local stress checks](ACCURACY.md)
show larger field disagreement; these profiles do not promise 1% local stress error.

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

For the studied Gyroid, an isolated native accuracy prototype now supports
`interior_refine_nodes`. Use the unchanged field above with `cells=48`,
`geometry_tolerance=0.0002`, `minimum_quality=0.1`, `optimize_passes=1`,
`threads=1` and `interior_refine_nodes=15196`. This option requires the verified
experimental wheel from `runs/accuracy-native-body-build-v1/package`, source
commit `83be74cb`; the released package and older PR #5 wheel do not expose it.
Import that isolated package before generating the mesh.

The profile produces 135,444 nodes and 621,403 tetrahedra, matching MMG's total
node count with 2.99% fewer tetrahedra. Its stiffness and mixed-load local-stress
scores are lower against both native80 and MMG64. Three measured live generation
runs give median 5.54 s versus standard MMG's 39.51 s, a 7.13-times ratio under
the same one-thread calibration. Minimum shape quality and wall-normal RMS
remain slightly worse than MMG, while all common gates pass. The independent
geometry admission tolerance remains 0.002; the native 0.0002 request concerns
its sampled forward boundary estimate. See [ACCURACY.md](ACCURACY.md) for both
distance directions, field scores, reference sensitivity and timing scope.
Local stress disagreement remains about 10%, rather than a certified error.

`runs/accuracy-native-matched-generation-v1/microgen-bridge.json` verifies the
actual matched-node and 20,000-node Gyroid meshes through public `BoxMesh`
roundtrips. Coordinates and tetrahedra remain exact, and both pass the public
periodic-node check at 1e-10. These checks do not change the geometry or run MMG.

The same isolated accuracy wheel also has an admitted Schwarz P profile at
matched total nodes. Use the fixed raw band [-0.6, 0.6] above with `cells=50`,
`geometry_tolerance=0.001` and `interior_refine_nodes=8168`, keeping the other
options fixed. The native request remains separate from the common geometry
admission gate of 0.002. This produces 123,172 nodes and 604,833 tetrahedra,
3.93% fewer tetrahedra than the standard MMG control at the same node count.
Stiffness and mixed-load local stress disagreement are lower against both
native80 and MMG64, including every stress scramble. Wall-normal RMS is worse
than MMG and local stress disagreement remains about 11%.

Three measured live Schwarz P generations give median 4.03 s versus 30.56 s
for standard MMG adaptation of a freshly generated native seed, a 7.57-times
ratio on one thread. All eight runs reproduce their admitted meshes exactly.
The public BoxMesh roundtrip and periodic-node checks also pass for the new
Schwarz P body and seed. See [ACCURACY.md](ACCURACY.md) for the fixed field,
reference sensitivity, unchanged gates and timing scope. These experimental
fixture profiles do not establish defaults for other TPMS or spinodoids.

Schwarz D also has an admitted matched-node baseline using Microgen's fixed
catalogue width 0.8502652628143951, `cells=40`, `geometry_tolerance=0.001` and
`interior_refine_nodes=11968` in the same experimental wheel. It matches MMG's
85,347 nodes with 12.40% fewer tetrahedra, and both stiffness and local stress
disagreement are lower against native80 and MMG64. Three measured live runs
give median 3.19 s versus 32.45 s for the standard MMG40seed48 recipe, a
10.17-times ratio on one thread. All eight meshes reproduce their admitted
arrays exactly; public Microgen roundtrips and periodic-node checks pass.
Stiffness disagreement remains 1.46-1.78%, and local stress disagreement is
15-16%. This profile remains an admitted baseline while the investigation
works toward the existing 1% stiffness selection target. See
[ACCURACY.md](ACCURACY.md) for the stronger MMG control and fixed gates.

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
