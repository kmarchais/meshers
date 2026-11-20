# Python API

The default wheel exposes volume meshing, field compilation and cancellation. `Field` is the callable type accepted by the mesher.

::: meshers.generate

::: meshers.generate_intersection

::: meshers.adapt_intersection

::: meshers.Mesh

::: meshers.compile_field

::: meshers.CompiledField
    options:
      members: ["__call__", "gradient", "values", "has_gradient", "node_count"]

::: meshers.CancellationToken

::: meshers.MeshingError

::: meshers.CancelledError

## Experimental surface build

With the `experimental-surfaces` Cargo feature, the Python package also exposes a triangle-only mesher for printing and preview.

::: meshers.surface.generate_surface

::: meshers.surface.SurfaceMesh
