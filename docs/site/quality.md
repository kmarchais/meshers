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

`sampled_surface_error` is the largest distance from a facet sample to the wall, measured by walking from the sample to the wall rather than by `|f| / |grad f|`. It bounds the distance from above whether or not the field is a distance function, and `geometry_tolerance` gates it. A sample in a box plane of an adaptive mesh is measured within that plane, unless the solid passes the plane by less than that.

`sampled_coverage_error` measures the other direction: from the wall of the solid, sampled where lines twice finer than the lattice cross it, three times with adaptation, to the meshed wall. Only this measure sees a feature the mesh lost or a crease it cut across. It is reported, capped at four tolerances, and does not reject a mesh.

Periodic images have bitwise identical coordinates in the plane of their faces, in linear and quadratic output.

## Volume of linear meshes

With every vertex on the wall, flat facets lie inside a convex wall and the mesh loses volume. On linear meshes each wall vertex therefore moves along the wall normal by the mean offset of its facets, by at most half of `geometry_tolerance`, and stays in its box planes. Wall vertices of a Tet4 mesh are then close to the wall, not on it. `balanced_wall_vertices` counts them. Quadratic meshes are unchanged.

## Quality recovery and narrow gaps

A mesh below `minimum_quality` is remeshed locally: short edges collapse, long ones split, interior edges and faces swap. With a quality floor set, the same operations also work on the elements below an MMG quality of 0.2 (`improved_elements`), and are kept when the final checks accept the result. Without a floor, the mesh is returned as cut.

If the floor is still missed, a last resort trades geometry for shape around the poor elements. Two crossings of the wall with one box edge merge, a crossing next to a box corner slides along its edge, and other wall vertices may leave the wall. Vertices stay in their box planes, periodic images move alike, and the result stays within `geometry_tolerance`. `defeaturing_moves` counts these moves and `generate` emits a `UserWarning` when it is nonzero.

Before cutting, lattice edges that hide a wall crossing deeper than a quarter of `geometry_tolerance` are split, so a gap or film thinner than a cell gets its own nodes. `refined_lattice_edges` counts the splits. An element left outside the solid by more than `geometry_tolerance` is removed when the boundary stays valid (`removed_bridging_elements`), otherwise the mesh is rejected. Adaptive meshing takes both steps too.

When the final checks refuse an adapted mesh, `generate` returns the uniform lattice mesh of the same field, sets `adaptation_abandoned` and emits a `UserWarning`.

Invalid options raise `ValueError`. A failed meshing gate raises `MeshingError`. [Python API](api/python.md)
