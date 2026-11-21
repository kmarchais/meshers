# Python API

<ol class="mesh-flow">
  <li><strong>Field</strong>Python callable or compiled expression</li>
  <li><strong>generate</strong>Set bounds, cells and geometry options</li>
  <li><strong>Mesh</strong>Read points, tetrahedra and surface arrays</li>
</ol>

The default wheel exposes volume meshing, field compilation and cancellation. `Field` is the callable type accepted by the mesher.

::: meshers.generate

`meshers.MAX_CONSTRAINTS` gives the named-constraint limit of the installed wheel.

::: meshers.generate_intersection

::: meshers.adapt_intersection

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
