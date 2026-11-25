# Periodic cap traces, 24 November 2026

Node-to-node periodic constraints tie degrees of freedom on opposite caps. Between the nodes, each cap interpolates those values over its own triangles. The old worker checked the displacement of every listed pair and counted the incident boundary triangles on each corner edge. Both checks can pass while a degree of freedom is unconstrained or the two finite-element traces differ.

## Controlled damage

Each control starts with the 24-cell gyroid from the release worker. Node and face IDs below are zero-based examples from the baseline Windows wheel.

| Change to returned mesh | What the old worker accepts | What the new worker checks |
| --- | --- | --- |
| Tet4, remove x pair (10, 105) | The other 333 x pairs translate by one. | Every node used by a low-side volume face needs an image. |
| Tet4, flip the high-side diagonal (101, 105) to (102, 103) | Two new positive tetrahedra fill the same pyramid. The surface still closes and every pair is intact. | The tagged triangles agree with their own volume on each side. The mapped low trace must also equal the high trace. |
| Tet10, remove midpoint pair (17889, 18237) | Corner connectivity and all remaining displacements are intact. | A quadratic trace needs the three edge degrees of freedom as well as its corners. |
| Tet10, exchange slots 3 and 4 of high-side face 158 | The triangle still has the same three corners and the same six node IDs. | The midpoint for edge (101, 102) must occupy that edge's slot in the volume trace. |

On the baseline Tet4 fixture the flip preserved integrated volume (0.3235671008) and minimum MMG quality (0.3148939565). Mean quality moved from 0.8088159 to 0.8088147, within the release gate's 1% allowance.

The worker selects each cap independently from the three tetrahedron corners on its unit-box plane. It gathers only the matching faces, compares them with the tagged surface, maps the source IDs through the pair table, and compares the complete triangle lists including multiplicity. The sort preserves each midpoint's opposite corner, so opposite winding and local numbering do not affect the comparison. A pair must have unique endpoints on the specified planes and satisfy the fixture's physical transform to absolute tolerance 1e-10. Wall-only contacts can have pairs even when they belong to no cap triangle.

The existing solver geometry study compares rounded transverse coordinates for linear triangles and imports SciPy and PyVista. The worker uses NumPy and integer node IDs so it can run against either installed wheel and distinguish different edge assignments sharing the same face nodes. [Solver geometry study](https://github.com/kmarchais/meshers/blob/aeb3ef09e0425d0dd9ed844aeb856ada54cd32cd/demos/solver/periodic_geometry.py#L208-L221)

## Windows results

Both supplied wheels passed all 11 original fixtures with both workers. On each wheel the old worker admitted all four mutations and the new worker rejected each one. The separate rotation control passed with 35 pairs on both wheels. Across the nine periodic fixtures the median three-axis check took 2.4 to 10.1 ms on the baseline wheel and 3.0 to 11.3 ms on the second wheel. The 70,403-element Tet10 gyroid took 9.8 ms on each.

The timing configuration was Windows 11 x64 on an Intel Core Ultra 9 285K with 24 logical CPUs and 63.4 GiB of RAM. The two supplied binaries are Windows x64 wheels. The other reachable Windows laptop used its Balanced power plan during this study. Its timing run remains pending. The JSON records alongside this report identify the installed native binary and store all five timing observations for each periodic fixture.

| Wheel source | SHA-256 of supplied wheel |
| --- | --- |
| `aeb3ef09e0425d0dd9ed844aeb856ada54cd32cd` | `bbb85214f42eaada78d42503a74bc8e0d1a94a898608af376f7299b29307096f` |
| `f9cd6efb7ba4c05ec23acd264ae9f847591f2b9b` | `d88a32c91beb4519fda10534009c9672cf79b440b79a3209558ccf0d2eae9b2f` |

## Scope and replay

The nine periodic release cases exercise all three unit translations. The graded sheet and mapped intersection run the usual nonperiodic checks. A separate cylindrical intersection uses eight cells and requests periodicity on the two angular caps. It selects those caps in parameter coordinates and rotates their physical nodes by 1/3 radian. The regression test also tries the identity transform.

The replay generates each mesh once. It passes the same float64 points, connectivity, generation time and high-water memory sample to the baseline and updated workers. For each original fixture it requires exact agreement of the result dictionaries. For each mutation it runs both complete sample validators. The five stored wall times cover the new three-axis check on the original mesh, after one discarded warm-up. Generation time is one observation per mesh. All replay cases share one process, so its high-water mark can include earlier cases. The production gate starts a fresh process per sample and records generation time and peak RSS before validation.

Use PowerShell 7 from the repository root. Set MESHERS_WHEEL to the absolute path of the wheel under test, BASELINE_WORKER to an unused absolute Python file path, and CAP_RESULT to the absolute output JSON path.

~~~powershell
git show aeb3ef09:tools/performance_worker.py > $env:BASELINE_WORKER
$env:OPENBLAS_NUM_THREADS = $env:OMP_NUM_THREADS = $env:MKL_NUM_THREADS = $env:RAYON_NUM_THREADS = '1'
uv run --no-project --python 3.12.15 --with $env:MESHERS_WHEEL --with numpy==2.5.5 --with psutil==7.2.2 --with pytest==9.2.0 -m tools.tests.replay_periodic_trace $env:BASELINE_WORKER $env:CAP_RESULT
~~~

The baseline worker is from `aeb3ef09e0425d0dd9ed844aeb856ada54cd32cd`. PR 22 at `fe3bbf3108dc2e09cc1a49cd29e10774de00a9e3` calibrates resource verdicts and their historical replay. The cap check runs inside the separately fingerprinted numerical worker. Archived scalar records retain the verdict of the worker that produced them. Rechecking a cap requires its mesh arrays. The combined driver, history, calibration and cap test suites passed 49 tests with PR 22 applied. [PR 22](https://github.com/kmarchais/meshers/pull/22)
