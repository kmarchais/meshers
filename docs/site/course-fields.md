---
icon: lucide/activity
---

# From an implicit field to a solid

<div class="guide-pair" markdown="block">

<figure class="guide-figure" data-guide="band">
  <div class="guide-diagram"><img src="assets/playground/band.svg" alt="A band narrows from left to right around a horizontal zero contour"></div>
  <div class="guide-control mesh-slider"><label for="setting">Band half-width</label><output for="setting"></output><input id="setting" type="range" min="5" max="50" value="30"></div>
  <figcaption>Section through f(x, y) = (1 + x)y. The two walls lie at y = ±b/(1 + x).</figcaption>
</figure>

<div markdown="block">

| Region | Keep the point when |
| --- | --- |
| Solid | `f ≤ 0` |
| Sheet | `−b ≤ f ≤ b` |

1. Set the half-width **b = 0.45**.
2. Compare the left and right ends. The zero line stays fixed.
3. Read the physical thickness at **x = 0** (0.90) and **x = 2** (0.30).

</div>
</div>

## How Meshers evaluates a field

| Trace this operation | At x = 2, y = 0.1 |
| --- | ---: |
| Add: `u = 1 + x` | 3 |
| Multiply: `f = u × y` | 0.3 |
| Differentiate: `∂f/∂x = y`, `∂f/∂y = u` | 0.1, 3 |

- **Compiled field:** reuse the traced operations and derivatives at each query point.
- **NumPy callback:** supply a gradient or use finite differences. A near-zero gradient makes wall projection unreliable.

[Trace the field](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-python/python/meshers/expression.py#L587) · [Differentiate the graph](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-python/python/meshers/expression.py#L406) · [Field options](fields.md)

## Explore the algorithms

| Follow the volume mesh | What you can investigate |
| --- | --- |
| [Cut the background](course-cutting.md) | Which edge crossings become mesh nodes? |
| [Adapt the elements](course-adaptation.md) | Where does a shorter edge pay off? |
| [Constrain and pair boundaries](course-constraints.md) | Which nodes may move together? |
| [Fit curved elements](course-quadratic.md) | How can moving one midside node fold a mapping? |
| [Extract surfaces](course-surfaces.md) | When do triangles or quads suffice? |
| [Verify the FEA result](course-verification.md) | Does refinement reduce physical error? |

<script type="module" src="assets/guide-controls.js"></script>
