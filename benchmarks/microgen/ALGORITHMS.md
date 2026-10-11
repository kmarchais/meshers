# Algorithm investigations for Microgen adoption

The benchmark should decide which implementation to improve. A different language
does not by itself make the same algorithm faster. Compare measured end-to-end
cost at the same geometry and solution-error target before choosing a replacement.

## Fast surfaces

Printing and visualization need a closed surface with accurate volume and
boundaries. Generating and optimizing a complete tetrahedral mesh to obtain that
surface wastes work. The surface mode measures Microgen's public surface path
and the full volume-to-boundary cost of Meshers 0.1.0 separately.

VTK's [Flying Edges implementation](https://vtk.org/doc/nightly/html/classvtkFlyingEdges3D.html)
uses four passes, edge ownership, computational trimming and preallocated output.
It is a useful baseline for a native structured-grid surface generator. VTK warns
that the filter can emit degenerate triangles, so extraction speed alone is not
an admission criterion. Clipped sheet walls also need consistent box caps and
shared feature curves. This benchmark does not substitute Flying Edges for
Microgen's public generator and then claim that it measured the current workflow.

## Accurate volume meshes

The original Meshers issue reproductions remain useful. In particular, extra
optimization passes can lose a previously satisfied geometry tolerance. A quality
objective should constrain surface error while improving element shape; moving
wall vertices back to the field does not constrain the edges and face interiors.
The benchmark samples those interiors and retains nonmonotonic outcomes.

[Isosurface stuffing](https://people.eecs.berkeley.edu/~jrs/papers/stuffing.pdf)
provides a grid-based reference with bounded dihedral angles under its algorithm's
assumptions. Its cut stencils and snapping rules deserve comparison with Meshers'
current extraction. Those guarantees do not transfer automatically to a different
stencil set, clipped sheet bands, arbitrary sharp unions or later optimization.

[HXT](https://arxiv.org/abs/2008.08508) provides a fast parallel Delaunay route to
quality tetrahedral meshes. It is a candidate for an implementation or additional
control after the surface is valid. Thread scaling, surface recovery and actual
FEA error should all enter the comparison.

[fTetWild](https://arxiv.org/abs/1908.03581) targets robust tetrahedralization of
imperfect triangle soups. It is relevant when imported or combined geometry
requires robustness. That input problem differs from extracting a known smooth
implicit function, so its published results do not establish that it will be
the fastest route for Microgen's procedural fields.

## MMG controls

MMG's [Hausdorff option](https://mmgtools.org/tutorials/tutorials_detailed_options.html)
constrains approximation to its reconstructed ideal boundary. It is not a promise
of that distance to the original analytic Microgen field. The independent volume
and two-direction boundary probes remain necessary after MMG remeshing.

The [level-set workflow](https://mmgtools.org/tutorials/tutorials_mmg3d_implicit.html)
starts from a background volume mesh carrying a sampled level set. It is a useful
future control alongside the current VTK-to-MMG pipelines. It needs its own name,
input sampling cost and retained parameters. It must not silently replace a
failed Microgen remeshing result.

The first adoption milestone is measured coverage of the low-error FEA targets
and the fast geometry targets across the fixed corpus. Published algorithm claims
suggest experiments; the retained benchmark results decide whether those claims
hold for these geometries, limits and implementations.
