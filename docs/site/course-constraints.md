---
icon: lucide/combine
---

# Keep cuts, seams and periodic pairs distinct

<div class="guide-pair" markdown="block">

<figure class="guide-figure"><img src="assets/cut-sphere.svg" alt="Copper planar cap meets the teal spherical wall along a circular rim"><figcaption>The wall has one active constraint, the cap has another, and the rim has both.</figcaption></figure>

<div markdown="block">

Intersect two sublevel sets to cut a sphere with a plane. Within a background tetrahedron, Meshers linearly interpolates each field and clips the convex cell against each constraint. It tetrahedralizes the surviving region after all cuts. The boundary labels preserve which field made each face.

**Predict a move.** A cap vertex can slide in the plane. At the rim the same move must also stay on the sphere. Projecting the rim onto only one wall could round the crease or pull the cap off its plane. [Isolate each surface in the viewer](intersections.md).

</div>
</div>

## Put a moved node back on its feature

Meshers reads the active fields from the node’s constraint mask. It linearizes their equations at the proposed point, solves the resulting tangent planes together, and repeats while holding locked box coordinates fixed. For flat polyhedral features the plane equations are already exact. For example, projecting onto both `x = 0.5` and `y = 0.5` sends `(0.7, 0.7, 0.2)` to `(0.5, 0.5, 0.2)`. The remaining coordinate is free along the edge. [Joint Newton projection](https://github.com/kmarchais/meshers/blob/main/crates/meshers-core/src/projection.rs) · [Affine projection](https://github.com/kmarchais/meshers/blob/main/crates/meshers-core/src/polyhedral.rs)

The same labels let Meshers adapt a supplied Tet4 seed in physical coordinates. For a mapped seed, the inverse map relates physical points to the chart where the constraint fields are defined. [Seed and map requirements](intersections.md#adapt-in-physical-coordinates)

## Opposite faces

<div class="guide-pair" markdown="block">

<figure class="guide-figure"><img src="assets/playground/pairs.svg" alt="Dashed translation lines join corresponding nodes on opposite box faces"><figcaption>A translated pair has two node IDs. Its tangential coordinates and the triangulation on the two caps must match.</figcaption></figure>

<div markdown="block">

For periodic volume adaptation, identify opposite copies in a cell that wraps around. Compute edge lengths using the nearest periodic image. Unfold once adaptation finishes to recover matching face meshes.

For a curved sector, the coordinate map acts on the background before cutting. A rigid transform pairs the ends in physical space. Apply its rotation to vector degrees of freedom. With a welded full turn, the end triangles form an interior seam and there is no external pair on that axis. [Maps and exported metadata](intersections.md#map-and-pair-the-ends) · [Intersection construction](https://github.com/kmarchais/meshers/blob/main/crates/meshers-core/src/implicit/intersection.rs)

</div>
</div>

**Next:** [Fit curved elements](course-quadratic.md).
