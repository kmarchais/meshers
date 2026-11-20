# Check the mesh before a solve

| Inspect | What to check |
| --- | --- |
| Orientation | Positive signed Tet4 volumes. For Tet10, exported `Volume` integrates the curved mapping. |
| Shape | MMG quality ranges from 0 (degenerate) to 1 (regular). `minimum_quality` rejects a nonempty mesh below your floor. Tet10 `CornerMMGQuality` describes its linear parent. |
| Boundaries | Tag 0 marks implicit faces. Tags 1/2, 3/4 and 5/6 mark low/high x, y and z caps. `sampled_surface_error` describes the linear parent. Curved faces need their own distance audit. |
| Periodicity | Request matching axes explicitly. Check `periodic_pairs` against the intended translation or [rigid transform](intersections.md). Transform the shift vectors if you transform the mesh. |
| Solver accuracy | Refine the mesh and measure the change in the physical quantity you need. [Elasticity example](solver-verification.md) |

## Missing material

On the development branch, set `expected_components` when you know the number of disconnected solids. A lost island or severed ligament can otherwise leave valid, positive elements.

- Meshers counts groups connected across triangular faces **before** joining periodic images. Vertex or edge contact does not join groups.
- Two defects can preserve the count. Use per-component volume and probes in thin features to check identity.
- A sampled-empty Tet10 result has connectivity shapes `(0, 10)` and `(0, 6)`. Require a positive component count when empty output would be an error.

## Topology and curved orientation

Meshers checks boundary edge incidence and vertex links, and retries cuts if snapping pinches the boundary. Periodic copies roll back together.

During Tet10 elevation and fitting on the development branch, determinant bounds cover the space between Jacobian samples. Inconclusive regions are bisected. A fold or unresolved bound limits the edge correction or rejects the fit. `minimum_jacobian_ratio` is the **sampled** ratio.

Invalid options raise `ValueError`. A failed meshing gate raises `MeshingError`. [Python API](api/python.md)
