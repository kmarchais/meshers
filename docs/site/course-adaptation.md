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

| Set chord length h | Read the sag on the unit circle |
| --- | ---: |
| 1.20 | 0.200 |
| 0.60 | 0.046 |

**Try it:** halve the chord. The gap shrinks by more than half.

- **Exact circle:** sag = `1 − √(1 − h²/4)`.
- **Mesh estimate:** edge length × angle between vertex normals / 6. Shorten boundary edges where sag exceeds tolerance.

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

1. Score the incident tetrahedra before and after a proposed move.
2. Reject an inversion. Move periodic copies together.
3. If interior optimization stalls, search for a move that improves the worst incident shape.

[2–3 swap](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/tet_adapt.rs#L1220) · [Recovery search](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/optimize.rs#L464)

**Test shape independently:** filter the [gyroid cutaway](quality.md) to its lowest scores. A sliver encloses little volume for its edge lengths.

## Refine to a budget

| Keep fixed | Where refinement acts |
| --- | --- |
| Periodic correspondence | Bisect all periodic copies of a selected surface edge together. |
| Boundary triangles | Add interior nodes up to the node budget. |
| Neighbour coordinates during scoring | Evaluate vertex moves with disjoint incident cells in parallel. |

[Surface and interior refinement](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/adaptive.rs)

**Next:** [Constrain and pair the boundaries](course-constraints.md).

<script type="module" src="assets/guide-controls.js"></script>
