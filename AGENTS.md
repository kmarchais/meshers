# Working on Meshers

Meshers generates and adapts meshes for finite element analysis. The primary workload is Microgen's implicit solids, including periodic TPMS and graded fields. Judge changes by the cost of resolving a physical problem to a stated accuracy.

## What to measure

- Check that the returned mesh represents the intended solid. Measure boundary distance in both directions, signed volume by component or material, retained features, and the boundary labels used by the solver. Include mapped caps and periodic seams. Record which checks are sampled.
- Check cell orientation, conformity and shape distributions on the final coordinates. For higher-order cells, inspect the curved mapping. Count elements, nodes and free degrees of freedom.
- For an FEA claim, state the load, material, boundary conditions, element formulation, reference solution and error tolerance. Compare assembly and solve costs as well as generation time and peak memory. Keep failed and over-budget cases in the results.
- Compare the same physical geometry on both routes. Identify the exact Microgen backend, seed, repair steps and solver settings. A fixed-count comparison answers a different question from a fixed-error comparison. Label each result accordingly.

Run the release fixtures and comparative FEA studies on every available, compatible configuration. Record any pending or unsupported runs alongside the completed ones. Capture CPU model, logical CPU count, physical RAM, OS and architecture when a benchmark runs. Publish dated configuration snapshots and numerical settings. Keep hostnames, local account paths and queue identifiers private. Record a new snapshot after a configuration change.

## How to make changes

Start with a reproducible discrepancy or profile on a pinned build. Preserve the failing input, baseline and replay command. Check the proposed change against an independent measurement and a case outside the tuning set before changing defaults. Record regressions in any of the competing goals. Reuse existing fixtures, solver studies and release benchmarks.

Preserve node ordering, material and face labels, boundary conditions and periodic pairs through adaptation and export. Test what the downstream solver reads. Keep diagnostics consistent with the returned mesh after fitting or repair.

CAD/BRep, triangle and quad surfaces, tetrahedra and hex cells are part of the intended scope. Identify which input and element combinations work in the installed build and which are research prototypes. For printing or visualization, state the separate acceptance criteria.

Use the development commands in [README.md](README.md) and the [performance protocol](docs/site/performance.md). Send heavy builds, sweeps and solver comparisons through the shared compute queue described in Kevin's machine instructions.
