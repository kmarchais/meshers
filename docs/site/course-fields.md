---
icon: lucide/activity
---

# From an implicit field to a solid

Follow a gyroid from a field evaluation to a stiffness calculation. The same steps apply to fields you write yourself.

<div class="guide-pair" markdown="block">

<figure class="guide-figure" data-guide="band">
  <div class="guide-diagram"><img src="assets/playground/band.svg" alt="A band narrows from left to right around a horizontal zero contour"></div>
  <div class="guide-control mesh-slider"><label for="setting">Band half-width</label><output for="setting"></output><input id="setting" type="range" min="5" max="50" value="30"></div>
  <figcaption>Section through f(x, y) = (1 + x)y. The two walls lie at y = ±b/(1 + x).</figcaption>
</figure>

<div markdown="block">

A field assigns one number to each point. With a sublevel set, keep points where `f ≤ 0`. With a sheet, keep points between two levels. The zero contour can stay fixed while the sheet thickness changes.

**Try it.** Set the half-width to 0.45. On the left the sheet is 0.90 units thick. On the right it is 0.30. A field-value tolerance only measures a distance directly for a signed distance field.

</div>
</div>

## How Meshers evaluates a field

For a traceable Python function, Meshers records arithmetic as a graph. In the example above it records an addition and a multiplication. Differentiating the same graph gives the gradient `(y, 1 + x, 0)`. The compiled evaluator can reuse these operations at many points.

A NumPy callback handles fields the tracer cannot express. Supply its gradient, or Meshers estimates it with finite differences. A near-zero gradient makes projecting onto a wall unreliable. [Python field options](fields.md) · [Expression tracer](https://github.com/kmarchais/meshers/blob/main/crates/meshers-python/python/meshers/expression.py)

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
