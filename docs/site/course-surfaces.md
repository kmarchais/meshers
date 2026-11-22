---
icon: lucide/grid-2x2
---

# Extract triangles or quads directly

<div class="guide-pair" markdown="block">

<figure class="guide-figure"><div class="guide-diagram">
<svg viewBox="0 0 464 202" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Triangle through three edge roots, and a quad joining four dual vertices around one crossed grid edge">
  <title>Primal crossings and dual connectivity</title>
  <g stroke="currentColor" opacity=".55" fill="none">
    <path d="M32 133 182 148 126 23 32 133 180 70 182 148 M126 23 180 70"/>
    <path d="M275 28H443M275 88H443M275 148H443M280 23V153M359 23V153M438 23V153"/>
  </g>
  <g stroke="var(--mint)" fill="var(--selected)" stroke-width="2.5">
    <path d="m79 78 27 62 1-39Z"/>
    <path d="m315 57 77 6 9 56-91-2Z"/>
  </g>
  <g fill="var(--chart-copper)">
    <circle cx="79" cy="78" r="4"/><circle cx="106" cy="140" r="4"/><circle cx="107" cy="101" r="4"/>
    <circle cx="315" cy="57" r="4"/><circle cx="392" cy="63" r="4"/><circle cx="401" cy="119" r="4"/><circle cx="310" cy="117" r="4"/>
  </g>
  <circle cx="359" cy="88" r="5" fill="var(--paper)"/>
  <text x="108" y="185" text-anchor="middle">Roots on edges</text><text x="359" y="185" text-anchor="middle">Vertices in cells</text>
</svg></div>
<figcaption>Left: one triangle cut from a tetrahedron. Right: end-on view of a crossed grid edge (dark dot). One local contour from each surrounding cell supplies a quad corner.</figcaption></figure>

<div markdown="block">

When only the skin is needed, experimental surface builds can skip the volume mesh. The triangle extractor divides grid cubes into tetrahedra with a consistent Freudenthal subdivision. It retains crossing edges and boundary polygons.

The quad extractor uses uniform dual contouring. It groups crossings into connected contour loops inside each cell and places a vertex at the mean crossing position of each group. A sign-changing grid edge connects vertices from the surrounding cells.

**Trace one face.** Count the four cells around the dot on the right. The quad crosses their boundaries. Its four vertices belong to local contour groups, so counting grid corners would give the wrong connectivity.

</div>
</div>

## Remesh the triangular surface

The surface remesher estimates a target edge length from curvature, then splits, collapses, flips and relaxes edges. Wall, cap and feature-curve vertices have different projection constraints. A final refinement pass checks sampled centroid and midpoint distances.

Check the installed build with `meshers.SURFACE_MESHING_AVAILABLE`. These experimental surface routes produce faces. A quad surface still needs a separate volume discretization for a solid FEA problem. Feature-aligned quad layouts and guaranteed recovery of subcell topology remain open. [Triangle extraction](https://github.com/kmarchais/meshers/blob/main/crates/meshers-core/src/triangles.rs) · [Surface adaptation](https://github.com/kmarchais/meshers/blob/main/crates/meshers-core/src/remesh.rs) · [Quad extraction](https://github.com/kmarchais/meshers/blob/main/crates/meshers-core/src/quads.rs)

**Next:** [Verify the finite element result](course-verification.md).
