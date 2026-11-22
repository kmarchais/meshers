---
icon: lucide/grid-2x2
---

# Where does a surface face come from?

<div id="surface-extraction">
<fieldset class="mesh-fields"><legend>Step through both extractors</legend>
  <div class="mesh-segments">
    <label><input type="radio" name="stage" value="0" checked><span>Sample signs</span></label>
    <label><input type="radio" name="stage" value="1"><span>Find crossings</span></label>
    <label><input type="radio" name="stage" value="2"><span>Place vertices</span></label>
    <label><input type="radio" name="stage" value="3"><span>Connect faces</span></label>
  </div>
</fieldset>
<div class="guide-pair">
<figure class="guide-figure" id="triangle">
  <h2>Triangle · cut one tetrahedron</h2>
  <div class="guide-diagram"><img src="assets/playground/extract-triangle.svg" alt="Three roots on edges from A form a triangle inside a tetrahedron."></div>
  <figcaption aria-live="polite">One corner below the plane gives three crossing edges and one triangle.</figcaption>
  <fieldset class="mesh-fields"><legend>Plane x + y + z = t</legend>
  <div class="mesh-segments">
    <label><input type="radio" name="level" value="0.5" checked><span>0.5</span></label>
    <label><input type="radio" name="level" value="1.5"><span>1.5</span></label>
    <label><input type="radio" name="level" value="3.5"><span>3.5</span></label>
  </div>
  </fieldset>
</figure>
<figure class="guide-figure" id="quad">
  <h2>Quad · circle one grid edge</h2>
  <div class="guide-diagram"><img src="assets/playground/extract-quad.svg" alt="Four dual vertices, one per cube, connect around the crossed central edge."></div>
  <figcaption aria-live="polite">Each cube supplies a vertex at the mean of its loop’s roots. Joining the four vertices gives one quad.</figcaption>
</figure>
</div>
</div>

## Predict the triangle count

1. Read the field at A, B, C, D: **0, 1, 2, 3**.
2. Set **t = 1.5**. Which edges bracket 1.5? On AD, the crossing is halfway from 0 to 3.
3. Count your crossings, then select **Connect faces**.

<details markdown="block"><summary>Check the crossing edges and face counts</summary>

| Plane | Edges with opposite signs | Emitted faces |
| --- | --- | --- |
| t = 0.5 | AB, AC, AD | One triangle |
| t = 1.5 | AC, AD, BC, BD | Two triangles sharing a diagonal |
| t = 3.5 | None | None |

</details>

## Follow the implementation

| Triangle extraction | Uniform dual contouring |
| --- | --- |
| [Split a cube into six tetrahedra](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/triangles.rs#L381-L420) | [Group a cell’s crossings into connected loops](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/quads.rs#L423-L465) |
| [Reuse each edge crossing](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/triangles.rs#L171-L239) | [Average the roots of each loop](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/quads.rs#L466-L481) |
| [Choose the triangle fan with the best worst score](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/triangles.rs#L310-L339) | [Join the four vertices around a crossed edge](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/quads.rs#L488-L525) |

## Remesh a curved patch

<details class="edge-split"><summary>Split a chord and project its new midpoint onto the wall</summary>
<svg viewBox="0 0 500 145" role="img" aria-label="A straight chord misses the arc. Its midpoint moves onto the arc, leaving two shorter chords.">
<title>One constrained edge split on a curved wall</title>
<g fill="none" stroke-width="2">
  <path d="M30 96 Q115 -28 200 96 M285 96 Q370 -28 455 96" stroke="var(--mint)"/>
  <path d="M30 96H200 M285 96H455" stroke="currentColor" stroke-dasharray="4 4"/>
  <path class="new-chords" pathLength="1" d="M285 96 370 34 455 96" stroke="var(--chart-copper)"/>
  <path d="M370 88V43m-5 7 5-7 5 7" stroke="currentColor"/>
</g>
<circle cx="115" cy="96" r="4" fill="var(--chart-copper)"/><circle cx="370" cy="34" r="4" fill="var(--chart-copper)"/>
<g fill="currentColor" text-anchor="middle"><text x="115" y="126">Midpoint on chord</text><text x="370" y="126">New node on wall</text></g>
</svg>
<ul><li>Project onto the active wall, cap or feature curve.</li><li>Recheck distances at edge midpoints and triangle centroids.</li></ul>
</details>

[Split](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/remesh.rs#L351) · [Collapse](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/remesh.rs#L387) · [Flip](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/remesh.rs#L457) · [Relax](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/src/remesh.rs#L496)

- **Python triangles:** check `meshers.SURFACE_MESHING_AVAILABLE`.
- **Rust quads:** experimental uniform extraction. Feature-aligned layouts and subcell topology remain open.
- **Solid FEA:** these methods emit surface faces. [Generate volume elements](course-cutting.md).

<script type="module" src="assets/extraction.js"></script>
