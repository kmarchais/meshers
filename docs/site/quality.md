# Periodicity and quality

## Request periodicity explicitly

```python
mesh = meshers.generate("gyroid", band=(-0.5, 0.5), cells=24, periodic=(True, True, True))
```

The field must match under translations across each requested pair of box faces.
Linear grading across an axis generally breaks that periodicity. Without a
constraint, matching boundary meshes are not guaranteed, even for a periodic field.

`periodic_pairs` contains matching node IDs for each axis. VTKHDF stores
`PeriodicMasterId` and physical `PeriodicShift` point arrays. Applications that
transform a mesh must transform its periodic shift vectors consistently.

## Inspect the result

| Measure | Interpretation |
| --- | --- |
| Positive tetrahedron volume | Necessary for a valid element orientation |
| MMG quality | 1 for a regular tetrahedron; approaching 0 indicates degeneracy |
| Sampled boundary error | Estimated geometric deviation in physical units |
| Volume distribution | Very small elements can increase FEM conditioning and timestep costs |

VTKHDF includes `Volume` and `MMGQuality`. Quality uses the MMG formula in NumPy;
it does not call mmgpy. Boundary tags are 0 for implicit faces, then 1 and 2 for x-/x+,
3 and 4 for y-/y+, and 5 and 6 for z-/z+.

Increase resolution and check solution convergence for the intended FEM model.
Minimum quality alone does not establish simulation accuracy. Small components
can be missed by sampling, and the geometry check is not a global proof.

Invalid options raise `ValueError`; unsuccessful geometry/topology checks raise
`MeshingError`. A `CancellationToken` allows cooperative cancellation, including
from another Python thread. Callback exceptions are preserved.

## Require a minimum element quality

```python
mesh = meshers.generate(
    "gyroid",
    band=(-0.5, 0.5),
    cells=24,
    periodic=(True, True, True),
    minimum_quality=0.1,
)
```

`minimum_quality` is a required lower bound on the MMG tetrahedron-quality
metric, in `[0, 1]`. The default `0` keeps the earlier geometry-only acceptance
policy. A nonempty mesh below the requested bound raises `MeshingError` and
returns no partial mesh. Empty geometries remain valid empty results. The same
option is available as `implicit::Options::minimum_quality` in Rust.

The example value 0.1 is an explicit acceptance choice, not a universal FEM
criterion. This option does not automatically remesh, change the field, or certify
solver accuracy. Appropriate shape quality, element sizes, conditioning and
solution convergence depend on the problem.

Diagnostics include `minimum_mmg_quality`, `elements_below_quality_01` (a count),
and `reverted_snap_vertices`. Positive volumes and finite positive quality are
checked regardless of the requested floor.

## Safe snapping and boundary checks

Snapping is retained where it passes the boundary checks. If snapped lattice
vertices introduce non-manifold boundary edges or vertex links, Meshers rolls
back those snaps and recomputes the cuts. Periodic copies are reverted together.
This repair is limited to eight clipping attempts; unresolved topology still
raises an error. It does not shift the cell phase, delete connected components,
or relax the boundary-error limit.

Both edge incidence and vertex links are checked before optimization and on the
finished mesh. Two shells touching at a single vertex are rejected even when
every edge has exactly two incident triangles. These combinatorial checks do
not prove that the mesh is free of geometric self-intersections or unresolved
features.
