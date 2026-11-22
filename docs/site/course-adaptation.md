---
icon: lucide/scaling
---

# Spend elements where curvature requires them

<div class="guide-pair" markdown="block">

<figure class="guide-figure" data-guide="sag">
  <div class="guide-diagram"><img src="assets/playground/sag.svg" alt="An arc, its chord and the maximum separation between them"></div>
  <div class="guide-control mesh-slider"><label for="setting">Chord length</label><output for="setting"></output><input id="setting" type="range" min="8" max="95" value="60"></div>
  <figcaption>Exact circle construction. The dashed segment measures the largest gap between the arc and its straight chord.</figcaption>
</figure>

<div markdown="block">

A flat triangle spans a curved patch. For a unit circle, a chord of length h has sag `1 − √(1 − h²/4)`. **Try it.** Compare h = 1.20 and 0.60. Halving the edge reduces the gap from 0.200 to about 0.046.

On a mesh, the local radius is unknown. Meshers estimates triangle sag from edge length ℓ and the change θ in vertex normals, approximately ℓθ/6. A tighter tolerance calls for shorter boundary edges. Target sizes grade into the interior.

</div>
</div>

## Change a neighbourhood

<figure class="guide-figure"><div class="guide-diagram">
<svg viewBox="0 0 480 178" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="A two-to-three swap in a five-vertex cavity">
  <title>A two-to-three swap in a five-vertex cavity</title>
  <defs><g id="cavity" fill="none" stroke="currentColor"><path d="M35 82 102 16 165 93 95 140Z M102 16 112 62 95 140 M35 82 112 62 165 93"/></g></defs>
  <use href="#cavity"/><use href="#cavity" x="245"/>
  <path d="m35 82 77-20 53 31Z" fill="var(--selected)" stroke="var(--mint)" stroke-width="2"/>
  <path d="m280 82 130 11" stroke="currentColor" stroke-dasharray="4 3"/>
  <path d="m347 16-7 124" stroke="var(--chart-copper)" stroke-width="3"/>
  <text x="100" y="166" text-anchor="middle">Two tets share a face</text><text x="345" y="166" text-anchor="middle">Three tets share an edge</text>
</svg></div><figcaption>A 2–3 swap keeps all five vertices and the cavity boundary. It changes their connectivity.</figcaption></figure>

| Local move | Effect |
| --- | --- |
| Split a long edge | Insert a node and retriangulate the incident cells to resolve detail. |
| Collapse a short edge | Remove a node when the remaining cells and constraints allow it. |
| Swap connectivity | Change which tetrahedra fill a cavity to improve shape. |
| Relocate a vertex | Improve neighbouring shapes. Project constrained motion back onto the wall or feature. |

Size and shape are separate objectives. The vertex optimizer scores the tetrahedra incident on a proposed move. If the smooth objective stalls at a poor interior vertex, a recovery search targets the worst incident shape directly. Inverting a cell invalidates the move. Translated copies of a periodic vertex move as a group. [Vertex optimization](https://github.com/kmarchais/meshers/blob/main/crates/meshers-core/src/optimize.rs) · [Adaptation loop](https://github.com/kmarchais/meshers/blob/main/crates/meshers-core/src/tet_adapt.rs)

**Inspect a sliver.** Filter the [gyroid cutaway](quality.md) to the low end of the shape scale. The score compares volume with that of a regular tetrahedron at the same root-mean-square edge length. Refining to smaller elements does not by itself raise this score.

## Refine to a budget

Surface-edge bisection splits the entire periodic orbit of a selected edge so that both sides remain conforming. Interior refinement can then insert nodes within a budget while keeping the boundary triangles fixed. This serves a different purpose from moving a wall closer to its level set. During vertex optimization, groups with disjoint incident cells can be evaluated in parallel without proposing two changes to the same tetrahedron. [Refinement methods](https://github.com/kmarchais/meshers/blob/main/crates/meshers-core/src/adaptive.rs)

**Next:** [Constrain and pair the boundaries](course-constraints.md).

<script type="module" src="assets/guide-controls.js"></script>
