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

| In the body-centred cubic grid | Make a tetrahedron from |
| --- | --- |
| Interior | Two cube centres + one shared-face edge |
| Box boundary | One centre + one triangle on the box |

1. Move the level to **t = 0.5**.
2. Trace the three edges from O to the copper dots.
3. Halving each edge leaves **½ × ½ × ½ = ⅛** of the volume.

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

1. Sample each lattice vertex. Snap nearby vertices to the wall with box-plane coordinates locked.
2. Bracket the root on each crossing edge. Use a secant step or bisect. Share the root with neighbouring tetrahedra.
3. Split each clipped cell using matching face diagonals. If snapping pinches the boundary, roll back the affected snaps and retry.

**Failure to try:** place a small void between lattice samples. Check the void centre as well as the mesh connectivity. [Cavity replay](benchmark-gaps.md#reproduce-the-lost-cavity) · [Build the background](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/implicit.rs#L1579) · [Clip a tetrahedron](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/implicit.rs#L1737)

**Next:** [Size and improve the elements](course-adaptation.md).

<script type="module" src="assets/guide-controls.js"></script>
