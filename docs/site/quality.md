# Check the mesh before a solve

| Inspect | What to check |
| --- | --- |
| Orientation | Positive signed Tet4 volumes. For Tet10, exported `Volume` integrates the curved mapping. |
| Shape | MMG quality ranges from 0 (degenerate) to 1 (regular). `minimum_quality` rejects a nonempty mesh below your floor. Tet10 `CornerMMGQuality` describes its linear parent. |
| Boundaries | Tag 0 marks implicit faces. Tags 1/2, 3/4 and 5/6 mark low/high x, y and z caps. `sampled_surface_error` describes the linear parent. Curved faces need their own distance audit. |
| Periodicity | Request matching axes explicitly. Check `periodic_pairs` against the intended translation or [rigid transform](intersections.md). Transform the shift vectors if you transform the mesh. |
| Solver accuracy | Refine the mesh and measure the change in the physical quantity you need. [Elasticity example](solver-verification.md) |

## Missing material

Set `expected_components` on the development branch to catch a lost island or severed ligament, even when every surviving element has positive volume.

- Meshers counts groups connected across triangular faces **before** joining periodic images. Vertex or edge contact does not join groups.
- Two defects can preserve the count. Use per-component volume and probes in thin features to check identity.
- A sampled-empty Tet10 result has connectivity shapes `(0, 10)` and `(0, 6)`. Require a positive component count when empty output would be an error.

## Topology and curved orientation

Meshers checks boundary edge incidence and vertex links, and retries cuts if snapping pinches the boundary. Periodic copies roll back together.

During Tet10 elevation and fitting, the development branch bisects regions where determinant bounds cannot rule out a fold. A fold or unresolved bound reduces the correction or rejects the fit. `minimum_jacobian_ratio` reports only sampled points.

## Geometry measures

`sampled_surface_error` is the largest distance from a facet sample to the wall, measured by walking from the sample to the wall rather than by `|f| / |grad f|`. A walk that reaches the wall bounds the distance from above whether or not the field is a distance function. A sample whose first-order estimate is under a quarter of `geometry_tolerance` keeps that estimate, as does a walk that does not converge. `geometry_tolerance` gates the result. A sample in a box plane of an adaptive mesh is measured within that plane, unless the solid passes the plane by less than that. Across a periodic plane the wall of the next cell counts as the wall inside does, so both images of a sample measure alike.

`sampled_coverage_error` measures the other direction: from the wall of the solid, sampled where lines twice finer than the lattice cross it, three times with adaptation, to the meshed wall. Only this measure sees a feature the mesh lost or a crease it cut across. On the lattice path, an edge whose ends have sharply different wall normals is checked for a fold between them: where the normal jumps over a short step, the refinement places its node in the fold, or part of the way down when the elements around the edge would flatten. Adaptive meshing does not do this yet. It is reported, capped at four tolerances, and does not reject a mesh.

Periodic images have bitwise identical coordinates in the plane of their faces, in linear and quadratic output. A lattice node snapped onto the wall is projected once for all its images, and every later operation moves images by the same steps.

## Volume of linear meshes

With every vertex on the wall, flat facets lie inside a convex wall and the mesh loses volume. On linear meshes each wall vertex therefore moves along the wall normal by the mean offset of its facets, by at most half of `geometry_tolerance`, and stays in its box planes. Wall vertices of a Tet4 mesh are then close to the wall, not on it. `balanced_wall_vertices` counts them. Quadratic meshes are unchanged.

## Quality recovery and narrow gaps

A mesh below `minimum_quality` is remeshed locally: short edges collapse, long ones split, interior edges and faces swap. A wall edge swaps for the other diagonal of its two facets where the solid is convex along it and the new facets stay within `geometry_tolerance` of the wall. With a quality floor set, the same operations also work on the elements below an MMG quality of 0.2 (`improved_elements`), and are kept when the final checks accept the result. Without a floor, the mesh is returned as cut.

If the floor is still missed, a last resort trades geometry for shape around the poor elements. Two crossings of the wall with one box edge merge, a crossing next to a box corner slides along its edge, and other wall vertices may leave the wall. Vertices stay in their box planes, periodic images move alike, and the result stays within `geometry_tolerance`. Each move is accepted on the samples and the measure of the final geometry check. `defeaturing_moves` counts these moves and `generate` emits a `UserWarning` when it is nonzero.

Before cutting, lattice edges that hide a wall crossing deeper than a quarter of `geometry_tolerance` are split, so a gap or film thinner than a cell gets its own nodes. `refined_lattice_edges` counts the splits. The search samples lattice edges: it sees gaps and films, not a ligament or channel thinner than a cell that passes between edges. It covers plain fields; meshes with a supplied background, retained constraints or half-spaces skip it, as they skip volume balancing and the coverage measure.

Where the search does not run, an element left outside the solid by more than `geometry_tolerance` is removed when the boundary stays valid (`removed_bridging_elements`), otherwise the mesh is rejected. An adapted mesh with such an element is refused.

When the final checks refuse an adapted mesh of a plain field, `generate` returns the uniform lattice mesh of the same field, without the requested sizes or layers, sets `adaptation_abandoned` and emits a `UserWarning`. Meshes with retained constraints or half-spaces raise `MeshingError` instead.

Each of these two degraded results can be refused instead of returned with a warning. `fallback=False` raises the refusal of the adapted mesh, which names the gate it missed. `defeaturing=False` keeps quality recovery to the operations that leave the wall where it is, and rejects a mesh they cannot lift to `minimum_quality`. `strict=True` turns both off; either can be turned back on beside it. A mesh that needs neither is the same whatever these are set to.

```python
mesh = meshers.generate(
    "gyroid", band=(-0.5, 0.5), cells=16, adapt=True, minimum_quality=0.1, strict=True
)
```

Invalid options raise `ValueError`. A failed meshing gate raises `MeshingError`. [Python API](api/python.md)
