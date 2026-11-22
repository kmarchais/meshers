---
icon: lucide/pyramid
---

# Cut a background tetrahedron

<div class="guide-pair" markdown="block">

<figure class="guide-figure" data-guide="cut">
  <div class="guide-diagram"><img src="assets/playground/cut.svg" alt="A triangular cut face joins three roots on the edges from the origin"></div>
  <div class="guide-control mesh-slider"><label for="setting">Cut level</label><output for="setting"></output><input id="setting" type="range" min="5" max="95" value="50"></div>
  <figcaption>Reference cell with corners O and the three unit axis points. Copper dots are the roots of x + y + z = t.</figcaption>
</figure>

<div markdown="block">

The volume mesher starts with a body-centred cubic lattice. Inside the box, two neighbouring cube centres and the endpoints of an edge on their shared face form a tetrahedron. At the box boundary a centre joins a face triangle.

**Try it.** Keep `x + y + z ≤ t`. At `t = 0.5` the kept tetrahedron has half the edge length along each axis and one eighth of the original volume. A clipped cell with two inside vertices has a different shape and needs several tetrahedra.

</div>
</div>

## Make adjacent cuts agree

<figure class="guide-figure"><div class="guide-diagram">
<svg viewBox="0 0 500 190" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="A tetrahedron between two cube centres">
  <title>A tetrahedron between two cube centres</title>
  <g fill="none" stroke="currentColor" opacity=".42">
    <path d="M65 56h300l55-27H120Z M65 143h300l55-27H120Z M65 56v87m150-87v87m150-87v87m55-114v87m-150-87v87m-150-87v87"/>
  </g>
  <path d="m215 56 55-27v87l-55 27Z" fill="var(--mint)" fill-opacity=".1" stroke="var(--mint)" stroke-dasharray="4 3"/>
  <path d="m168 86 47-30v87Z m150 0-103-30v87Z" fill="var(--selected)" stroke="var(--mint)" stroke-width="2"/>
  <path d="M168 86h150" stroke="var(--mint)" stroke-dasharray="4 3"/>
  <path d="M215 56v87" stroke="var(--chart-copper)" stroke-width="4"/>
  <g fill="var(--paper)"><circle cx="168" cy="86" r="4"/><circle cx="318" cy="86" r="4"/></g>
  <text x="250" y="175" text-anchor="middle">Two cell centres + one shared face edge</text>
</svg></div><figcaption>Four choices of edge on the shared square face give four tetrahedra between these centres.</figcaption></figure>

1. Sample the field at lattice vertices. A constrained snap can move a nearby vertex onto the wall while keeping its box coordinates on their planes.
2. Locate crossings on bracketed edges. A secant step uses the endpoint values. Bisection keeps the search inside the bracket when that step comes too close to an end. Neighbours reuse the same edge root for the same level.
3. Split each clipped region into tetrahedra with compatible face diagonals. Inspect boundary edges and vertex links. A pinched boundary can trigger a retry with affected snaps rolled back.

A small inclusion can sit between lattice samples. A conforming mesh alone cannot establish that it contains that inclusion. [Thin-feature controls](benchmark-gaps.md) · [Volume construction](https://github.com/kmarchais/meshers/blob/main/crates/meshers-core/src/implicit.rs)

**Next:** [Size and improve the elements](course-adaptation.md).

<script type="module" src="assets/guide-controls.js"></script>
