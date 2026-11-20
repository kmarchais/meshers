# Python API

The default wheel exposes volume meshing, field compilation and cancellation. `Field` is the callable type accepted by the mesher.

::: meshers.generate

`meshers.MAX_CONSTRAINTS` gives the named-constraint limit of the installed wheel.

::: meshers.generate_intersection

::: meshers.adapt_intersection

### Intersection keyword options

Pass these directly to `generate_intersection` or `adapt_intersection`, for example
`adapt_intersection(seed, fields, bounds=box, max_length=0.05)`.

| Keyword | Default | Meaning |
| --- | --- | --- |
| `seed_tolerance` | `0.01` | Maximum sampled distance for the starting mesh during adaptation. |
| `min_length`, `max_length` | Final geometry tolerance, shortest physical box side / 4 | Smallest and largest target lengths during adaptation. |
| `gradation`, `adapt_iterations` | `0.3`, `4` | Maximum target-length growth per unit distance and number of size-driven rounds. |
| `quality_guard`, `adapt_tolerance` | `0.3`, final geometry tolerance | Local quality floor and distance target for adaptive operations. |
| `optimize_passes` | `4` | Constrained optimization rounds (0–20). |
| `max_tetrahedra` | Generate 2,000,000 · adapt 250,000 | Output element budget. |
| `compile_expression`, `batch_size` | Automatic compilation enabled, `4096` | Trace supported fields and maps. Evaluate callbacks in batches of at most this size (1–65,536). |
| `cancel` | `None` | Cooperative cancellation token. |

| Operation | Additional keyword | Default | Meaning |
| --- | --- | --- | --- |
| Generate | `coordinate_map` | Identity | Map parameter coordinates to physical coordinates before cutting. |
| Generate | `periodic`, `closed_axes` | Three false flags each | Pair opposite caps on periodic axes. Close coincident seams before repair on closed axes. |
| Generate | `periodic_transforms` | Box translations for unmapped grids | Map each periodic axis's low face to its high face with a rigid 4×4 matrix. Mapped periodic axes require an explicit transform. Closed axes require identity. |
| Generate | `adapt`, `inverse_map`, `physical_bounds` | `False`, `None`, `None` | Adapt a mapped seed in physical space. Supply both maps and physical bounds. External periodic faces are unsupported on this path. |
| Generate | `snap` | `0.2` | Joint background snap fraction (0–0.2). |
| Adapt | `constraint_masks` | Seed masks | Override per-vertex labels. Bit `j` identifies field `j` in the supplied mapping. |

::: meshers.Mesh

::: meshers.compile_field

::: meshers.CompiledField
    options:
      members: ["__call__", "gradient", "values", "has_gradient", "node_count", "halfspaces", "polyhedral_regions", "intersection_fields", "sdf_intersection"]

::: meshers.CancellationToken

::: meshers.MeshingError

::: meshers.CancelledError

## Experimental surface build

Check `meshers.SURFACE_MESHING_AVAILABLE` before calling the triangle-only mesher. It requires a wheel built with the `experimental-surfaces` Cargo feature.

::: meshers.surface.generate_surface

::: meshers.surface.SurfaceMesh
