# Tetrahedral meshes from implicit solids

![Triangulated gyroid sheet. Teal triangles follow the two implicit walls and copper triangles close the sheet at the box faces.](assets/gyroid.png)

Teal faces follow the gyroid walls. Copper faces cap the box. Use either set to assign boundary conditions.

```python
import meshers

mesh = meshers.generate("gyroid", band=(-0.5, 0.5), cells=12, periodic=(True, True, True))
```

- [Generate and export your first mesh](getting-started.md)
- [Keep a cut face or pair periodic boundaries](intersections.md)
- [Check the mesh](quality.md) and [see how displacement error changes with refinement](solver-verification.md)
