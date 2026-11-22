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

Tet10 adds one shared midside node to each of the six edges of a Tet4. Moving wall nodes can fit curvature while keeping the corner nodes fixed. Straight corner tetrahedra can have positive volume even when the curved map turns back inside an element.

**Try it.** The edge map is `x(t) = t + 4dt(1 − t)`. Its midpoint slope is always 1. At a displacement of 25% the slope at the right endpoint reaches zero. Beyond that point two different reference positions can map to the same physical position.

</div>
</div>

## Check the whole reference cell

The 3D test asks whether the Jacobian determinant stays positive. For a quadratic element each Jacobian row varies linearly in reference coordinates. Meshers bounds the determinant using triple products of the corner rows. A strictly positive lower bound admits that region. An inconclusive bound triggers subdivision. Exhausting the subdivision budget leaves the fit uncertified.

The fitting step can reduce a correction or reject it. The reported `minimum_jacobian_ratio` is a sampled diagnostic. Use the positivity check to establish orientation, and a separate boundary check to measure fit. [Wall fitting](https://github.com/kmarchais/meshers/blob/main/crates/meshers-core/src/quadratic/fit.rs) · [Determinant bound](https://github.com/kmarchais/meshers/blob/main/crates/meshers-core/src/quadratic/positivity.rs)

**Next:** [Extract a surface mesh](course-surfaces.md).

<script type="module" src="assets/guide-controls.js"></script>
