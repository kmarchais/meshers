# Intersections and curved periodic boundaries

Use `generate` for one implicit field or a band. Use `generate_intersection`
when distinct surfaces meet and their intersection curves must remain mesh edges.
Both return the same `Mesh` type. Custom grading, blends, lattices and random
fields use the existing callable/compiled-field interface.

```python
import meshers
import numpy as np

mesh = meshers.generate_intersection(
    {
        "sphere": lambda x, y, z: x * x + y * y + z * z - 0.7**2,
        "cut": lambda x, y, z: z - 0.25,
    },
    bounds=(-0.8, 0.8, -0.8, 0.8, -0.8, 0.8),
    cells=24,
    optimize_passes=8,
    geometry_tolerance=0.02,
)
cut_triangles = mesh.surface_for("cut")
intersection_edges = mesh.feature_edges
mesh.write_vtkhdf("cut-sphere.vtkhdf")
```

Each named constraint keeps its negative side. A sheet uses two constraints,
`field - upper` and `lower - field`. Constraints are intersected inside each
background cell before tetrahedralization. Quality-aware construction and
constrained optimization preserve surface identities and feature topology.
If joint snapping pinches the solid, it is backed out and the same geometry/grid
is retried; `reverted_snap_vertices` records the rollback.

## Parameter maps and periodicity

The optional `coordinate_map` returns three expressions in parameter coordinates.
It maps the **background before cutting**, avoiding inversion from bending thin
finished tetrahedra. It must preserve orientation and be one-to-one away from paired periodic ends;
local orientation checks do not certify the absence of global self-intersections.
Fields, `bounds`, `cells` and
`periodic` describe the parameter domain; returned `points` are physical and
`parameters` retain their original coordinate meaning.

For a torus sector with major radius R, axial parameter z becomes angle z/R:

```python
R = 3.0
length = 1.0
angle = length / R
rotation = np.eye(4)
rotation[:2, :2] = [
    [np.cos(angle), -np.sin(angle)],
    [np.sin(angle), np.cos(angle)],
]

mesh = meshers.generate_intersection(
    {"radial": lambda x, y, z: x - 0.5, "width": lambda x, y, z: y - 0.6},
    bounds=(0, 1, 0, 1, 0, length),
    cells=12,
    coordinate_map=lambda x, y, z: ((R + x) * np.cos(z / R), (R + x) * np.sin(z / R), -y),
    periodic=(False, False, True),
    periodic_transforms={2: rotation},
    geometry_tolerance=0.02,
)
```

`periodic_transforms` maps each enabled axis (0, 1 or 2) to a rigid 4x4 matrix
from the low end to the high end. For unmapped geometry, box translations are
automatic. Mapped periodic axes require explicit transforms: generation checks
both node correspondence and matching end triangles after optimization.
For vector problems, rotate vector degrees of freedom using the matrix's 3x3
block; coordinate shifts alone do not express rotational periodicity.

A complete cylinder/60-degree torus TPMS example is in
[`curved_tpms.py`](https://github.com/kmarchais/meshers/blob/main/crates/meshers-python/examples/curved_tpms.py).

## Inspection and export

- `constraint_names`: names in insertion order.
- `constraint_masks`: per-vertex bit j identifies membership in constraint j.
- `surface_for(name)`: triangles assigned to a named surface.
- `feature_edges`: boundary edges shared by two or more named constraints.
- `parameters`, `periodic_pairs`, `periodic_transforms`: coordinate and pairing metadata.
- `diagnostics`: quality before/after optimization, accepted operations, volume,
  timing, callbacks and sampled geometry errors when requested.

VTKHDF retains parameter coordinates and constraint masks as point arrays,
paired node IDs and transform matrices as field arrays, and feature edges.
`ConstraintNamesUTF8` is a uint8 field array containing UTF-8 JSON, readable by
VTK without HDF5 string conversion assumptions.

## Acceptance and limits

`minimum_quality` is a strict acceptance gate; unmet requests raise
`MeshingError`. It does not round sharp junctions or silently relax the target.
Very acute input angles can make a requested minimum impossible.

`geometry_tolerance=None` skips the optional physical surface check. A positive
value samples facet vertices, edge midpoints and centroids against the assigned
analytic constraints. Mapped checks invert the coordinate map locally and
transform gradients to physical coordinates. Singular/inconsistent inversion
is rejected. This is a first-order distance estimate and an acceptance gate,
not adaptive refinement or a Hausdorff bound.

The constrained path supports 2–8 fields, 4–128 intervals per axis, 0–20 quality
passes and one native worker. Straight-sided elements approximate curved walls;
it does not produce high-order curved finite elements. Callbacks work but
scalar evaluation can be expensive; supported expressions compile automatically.
Sampling cannot certify sub-grid features. Solver qualification remains
application-specific; see the [elasticity verification](solver-verification.md).
