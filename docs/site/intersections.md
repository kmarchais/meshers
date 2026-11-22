---
icon: lucide/layers
---

# Named surfaces and curved domains

<div class="guide-pair" markdown="block">

<figure class="guide-figure" id="cut-layers">
  <img src="assets/cut-sphere.svg" alt="Cut sphere with a triangulated wall and planar cap">
  <fieldset><legend>Visible boundary</legend>
  <label><input type="radio" name="boundary" value="" checked> Both</label>
  <label><input type="radio" name="boundary" value="wall"> Spherical wall</label>
  <label><input type="radio" name="boundary" value="cap"> Planar cut</label>
  </fieldset>
  <figcaption>The triangles come from the two named surfaces of a 12-cell-per-axis mesh.</figcaption>
</figure>

<div markdown="block">

Each field names a boundary and keeps its negative side. Here `cut` labels the copper face at `z = 0.25`.

```python
import meshers

mesh = meshers.generate_intersection(
    {
        "sphere": lambda x, y, z: x**2 + y**2 + z**2 - 0.7**2,
        "cut": lambda x, y, z: z - 0.25,
    },
    bounds=(-0.8, 0.8) * 3,
    cells=12,
    optimize_passes=2,
    geometry_tolerance=0.04,
)
cut_triangles = mesh.surface_for("cut")
rim_edges = mesh.feature_edges
```

VTKHDF preserves the named surfaces, per-vertex `constraint_masks`, feature edges, chart coordinates and periodic metadata.

</div>
</div>

<script type="module" src="assets/guide-controls.js"></script>

## Map and pair the ends

| Geometry | What to supply | What the solver receives |
| --- | --- | --- |
| Straight periodic cell | `periodic` on matching axes | Low/high node pairs and box translations |
| Curved sector | An orientation-preserving, one-to-one `coordinate_map` and a rigid 4×4 low-to-high transform per periodic axis | Physical points, chart coordinates and transforms. Use the 3×3 rotation block for vector degrees of freedom. |
| Full turn | `closed_axes` on the periodic seam and the identity transform | Matching end triangles weld into an interior seam. The welded axis has no external pairs. |

The map acts on the background **before cutting**. Fields, bounds and cell counts use chart coordinates. [Cylinder and torus example](https://github.com/kmarchais/meshers/blob/main/crates/meshers-python/examples/curved_tpms.py)

## Adapt in physical coordinates

- **Mapped seed** — pass `adapt=True`, `inverse_map` and `physical_bounds`. Resolve each feature with positively oriented cells and meet `seed_tolerance` (default 0.01).
- **Existing Tet4 seed** — pass physical fields and matching masks to `adapt_intersection`. Physical adaptation uses one worker on nonperiodic solids. Returned `parameters` are physical coordinates.
- **Final mesh** — require `minimum_quality` and `geometry_tolerance`. Meshers samples walls and chart caps, including triangle interiors. Inspect `mapped_cap_errors` per cap. Check between-sample error and global self-intersections separately.

[Python options and limits](api/python.md) · [Quality checks](quality.md)
