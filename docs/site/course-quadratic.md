---
icon: lucide/spline
---

# Fit a curved tetrahedron without folding it

<div class="guide-pair" markdown="block">

<figure class="guide-figure" data-guide="fold">
  <div class="guide-diagram"><img src="assets/playground/fold.svg" alt="Graph of physical position against reference position, with endpoints and the movable midpoint marked"></div>
  <div class="guide-control mesh-slider"><label for="setting">Midside displacement (% of edge)</label><output for="setting"></output><input id="setting" type="range" min="0" max="50" value="10"></div>
  <figcaption>One-dimensional analogue. Dashed line: straight map. Solid curve: quadratic map through the three marked nodes.</figcaption>
</figure>

<div markdown="block">

- **4 corners** fix the parent tetrahedron.
- **6 shared midside nodes** control edge curvature.
- **One edge in the widget:** `x(t) = t + 4dt(1 − t)`. Its midpoint slope is always 1.

| Set displacement d | Inspect the right end |
| --- | --- |
| 10% | The curve rises everywhere. |
| 25% | The endpoint slope reaches zero. |
| 35% | The curve turns back. Two reference points can map to the same position. |

</div>
</div>

## Check the whole reference cell

| Whole-cell test | Action |
| --- | --- |
| Bound the determinant with triple products of the affine Jacobian rows | Accept the region if its lower bound clears the guard. |
| Bounds inconclusive | Bisect and check both children. |
| Fold detected or subdivision budget exhausted | Reduce the correction or reject the fit. |

**Read the diagnostic:** `minimum_jacobian_ratio` samples the mapping. The whole-cell bound checks orientation. Measure boundary fit separately. [Wall fitting](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/quadratic/fit.rs) · [Subdivision and bound](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/quadratic/positivity.rs#L23)

**Next:** [Extract a surface mesh](course-surfaces.md).

<script type="module" src="assets/guide-controls.js"></script>
