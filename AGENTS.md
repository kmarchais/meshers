# Working on Meshers

Meshers is almost entirely vibe coded. Use reproducible benchmarks to guide agents as they extend support for implicit geometries (FRep), especially Microgen’s periodic and graded solids. For finite element analysis (FEA), compare the cost of reaching a stated physical accuracy. For CAD and new cell types, distinguish released support from prototypes.

## Measure algorithm changes

Compare pinned base and candidate builds on the same geometry, hardware and settings. Preserve the input, seed, backend, repair steps and replay command. Choose a suite from the [benchmark index](docs/site/benchmark-index.md) and check its [open targets](docs/site/benchmark-gaps.md). Reuse the release fixtures, add a regression case for the feature or defect, and check a case outside the tuning set before changing defaults.

| Metric | What to record |
| --- | --- |
| Computation | Generation, assembly and solve time separately, plus peak memory. |
| Boundary and volume accuracy | Sample distances in both directions, including caps and seams. Check signed volume by component or material and probe thin features. |
| Mesh quality and size | Orientation, conformity and shape distributions on final coordinates. Check curved mappings. Count elements, nodes and free degrees of freedom. |
| FEA convergence | Error against a reference under refinement. State the load, material, boundary conditions, element formulation and tolerance. |

Verify node ordering, labels, boundary conditions and periodic pairs in the exported mesh the solver reads. Run release fixtures and comparative FEA studies on available compatible configurations. Keep failures, fallbacks and pending runs visible. Record dated snapshots of CPU model, logical core count, physical RAM, OS and architecture and keep machine names private.

## Checks and pull requests

Every PR that adds a feature or changes how an algorithm works must update the documentation in the same PR. Update the relevant user guide or [features and algorithms](docs/site/course-fields.md) pages, including any examples, figures and limits affected by the change.

Use the pinned formatter, linter, type checker and tests to guide edits. Run Ruff and ty for Python, rustfmt and Clippy for Rust, and the relevant regression tests. Follow the [development commands](README.md#develop) and [benchmark protocol](docs/site/performance.md). Queue heavy builds and studies on the shared compute queue.

Keep PR descriptions to a purpose sentence and the evidence. An algorithm PR should add a capability or improve a measured result. Include a compact before/after table for the affected metrics, a replay command or results link, and regression results for existing features. Quantify any tradeoff so reviewers can decide whether it is worthwhile. For documentation-only changes, show the changed page or example.
