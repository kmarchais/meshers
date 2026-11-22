# Fields, thickness and grading

For a solid, Meshers keeps `f ≤ 0`. For a band, it keeps `lower ≤ f ≤ upper`. A fixed band gets thinner where the field changes faster.

<figure class="guide-figure" data-guide="band">
  <div class="guide-diagram"><img src="assets/playground/band.svg" alt="Analytic section of a graded sheet"></div>
  <div class="guide-control mesh-slider"><label for="band-setting">Band half-width</label><output for="band-setting"></output><input id="band-setting" type="range" min="5" max="50" value="30"></div>
  <figcaption>Vertical section through f = (1 + x)y, from x = 0 to x = 2.</figcaption>
</figure>

<script type="module" src="assets/guide-controls.js"></script>

To grade a sheet, divide the field by a positive half-width and use `band=(-1, 1)`. Match widths across periodic ends. [Graded sheet example](https://github.com/kmarchais/meshers/blob/main/crates/meshers-python/examples/graded.py)

| Field implementation | Evaluation and derivatives |
| --- | --- |
| Arithmetic, supported NumPy functions, comparisons, `np.where` | Meshers traces and compiles the expression and its gradient, including derivatives of grading terms. |
| Unsupported operation or `compile_expression=False` | Batched NumPy callback returning float64 `(N,)`. Supply gradients of shape `(N, 3)` or use finite differences. |
| `compile_field(f)` | Reusable, immutable snapshot. Automatic preparation retraces captured values and caches up to 32 compiled expressions. |

**Nonsmooth fields.** `abs` has zero derivative at zero. Min/max ties average both gradients. `where` uses the selected branch.

For angular fields, keep the inside/outside classification consistent at the axis and across the seam. [Evaluator examples](https://github.com/kmarchais/meshers/blob/main/crates/meshers-python/examples/evaluators.py)

The original `compile=` keyword is accepted as a deprecated alias for `compile_expression=`.
