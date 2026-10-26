# Native density lookup and curved cap-fitting results

26 October 2026. Static `np.interp` tracing is implemented in public Meshers
Python code. It reduces generation time on the fixed linear-density Gyroid.
A separate experimental Python fitting pass admits a smaller coarse Tet10
mesh. That route sacrifices boundary accuracy and some quality quantiles.
Neither study establishes superiority over an admitted actual Microgen+MMG
density route, or matched FEA accuracy.

[Numeric findings](../../handoff/macbook/density-lookup-findings-20261026.json)
contain full precision, source hashes, measured directories and independent
replay receipts. [The dashboard](comparison-dashboard.html) shows the current
public recipes, earlier callback recipes and the experimental fitting tradeoff.

## Native interpolation removes callback overhead

The Python tracer now represents a static, finite, real, strictly increasing
interpolation table with existing native comparison, selection and arithmetic
operations. It preserves the frozen piecewise-linear table, endpoint clamping
and finite constant `left`/`right` values. Unsupported tables retain callback
fallback. It adds no dependency or native opcode. The table size is bounded by
the graph limit; eager graph evaluation still costs O(table size).

The public default uses analytic derivatives: the interval to the right at
interior knots and zero at the endpoints. `compile_field(field, gradient=False)`
retains numerical derivatives. These policies are measured separately because
changing derivatives can change the mesh.

Every table row below has one excluded warmup and five interleaved measured
pairs. All 24 attempts in each study pass admission. Compilation is charged
inside generation. Times are fresh M4 measurements with one numerical thread;
common imports, fixture loading, output archives and independent audits are
excluded for both arms. Native diagnostics report zero Python seed-evaluator
callbacks, versus about two million in each callback generation.

| Order / gradient policy | Callback median s | Native median s | Speedup | Native tets / nodes |
|---|---:|---:|---:|---:|
| Tet4, numerical | 10.893693 | 5.934124 | 1.835771x | 132,812 / 31,318 |
| Tet10, numerical | 11.300277 | 6.048265 | 1.868350x | 141,392 / 218,792 |
| Tet4, analytic default | 10.924949 | 3.947524 | 2.767545x | 132,822 / 31,330 |
| Tet10, analytic default | 11.303363 | 4.153106 | 2.721665x | 141,437 / 219,006 |

The numerical-derivative pairs have identical counts and mapped corner cells
and boundary faces. Their coordinates are not bitwise identical. Tet4 maximum
nearest-point difference is 1.2437503e-10, which fails the separate 1e-10
equivalence diagnostic; Tet10 is 2.9353244e-12. That diagnostic is retained,
not relaxed. Independent geometry observations change by about 1e-15 and
quality quantiles by about 1e-13. These facts support an evaluator gain, not
an exact-domain claim or transfer of old FE fields.

The public automatic field with numerical derivatives is the same cached
compiled object as the measured prototype selection tree. A separate 8,192
point check confirms exactly equal compiled values. The analytic default
changes the mesh counts by at most 0.032% tets and 0.098% nodes. Its geometry
and quality changes remain visible:

| Order / evaluator | Boundary error | Abs. relative volume error % | q min | q p01 | q median | q p99 |
|---|---:|---:|---:|---:|---:|---:|
| Tet4 callback | 0.0009922965 | 0.12241498 | 0.261089 | 0.528952 | 0.753412 | 0.973784 |
| Tet4 analytic native | 0.0009966085 | 0.12182071 | 0.261206 | 0.529630 | 0.752532 | 0.973559 |
| Tet10 callback | 0.0005837605 | 0.02034151 | 0.266349 | 0.533461 | 0.732800 | 0.966991 |
| Tet10 analytic native | 0.0005888778 | 0.02036104 | 0.271434 | 0.531535 | 0.733067 | 0.967437 |

Boundary error is a sampled two-sided observation normalized by cell length.
The numeric receipt retains all seven quality quantiles. Small volume
differences are not resolved accuracy wins against the empirical reference
scramble spread.

## One cap edge explains the coarse Tet10 rejection

The earlier callback coarse mesh has 25,787 tets and fails distance at
0.002020511. The fresh public-analytic uniform16 control has 25,747 tets and
45,262 nodes; it still fails, at 0.002099825. The evaluator feature alone does
not fix this coarse recipe. Both rejections remain frozen.

Localization of the new mesh finds one above-gate forward sample on a wall
edge intersecting a box cap. An independent constrained-foot calculation
still reports 0.002099661, so a projection-observer change cannot admit it.
The midpoint has one legal tangent degree of freedom. Eighteen lookup knots
lie across the edge's z span.

A saved-mesh feasibility fit moves that midpoint by 0.002003837 along its
allowed direction, below the fixed 10%-edge-length bound. It lowers the
16-Gauss-point edge objective from 0.002142696 to 0.001770301 and an independent
65-point edge maximum from 0.002170073 to 0.001801589. Whole-mesh admission then
passes at distance 0.001860206. Counts, connectivity, corners and the complete
corner-quality vector remain unchanged. This changes the curved approximation
and FE space; it is not a nested correction or a live workflow timing.

## A fresh generic fitting workflow reproduces the gain

The next experiment generates a fresh uniform16 native seed and inspects every
cap-intersection wall edge. Its marker uses 16 Gauss points, the physical field
and the unchanged distance threshold. It uses no fixed benchmark samples,
saved seed, known failing edge ID or FE field. It finds one of 694 edges.

The fitting policy allows at most eight marked edges, one pass, disjoint
incident stars and a bounded one-dimensional tangent optimization. It uses the
same scalar displacement limit, 48-iteration limit, 65-point non-regression
check and positive 125-point Jacobian checks along the move. After fitting,
the marker reports zero remaining above-threshold edges. Whole-mesh admission
uses the independent fixed forward/reverse samples and the full original gates.

The pilot and independent replay pass. A new repeated study has 12 successful
attempts: two excluded warmups and five interleaved measured pairs. An
independent replay checks both unique outputs. The fitting route's measured
time includes its Python import, detection, optimization, validity checks and
conversion. It is a hybrid experiment, not a native fitting default.

| Metric | Adaptive analytic native | Uniform16 + Python cap fit |
|---|---:|---:|
| Generation median s | 3.996870 | 1.103007 |
| Generation range s | 3.869644–4.069762 | 1.082573–1.120901 |
| Tets | 141,437 | 25,747 |
| Nodes | 219,006 | 45,262 |
| Boundary error | 0.0005888778 | 0.0018602057 |
| Abs. relative volume error % | 0.02036104 | 0.02322216 |
| q min | 0.271434 | 0.163388 |
| q p01 | 0.531535 | 0.433258 |
| q p10 | 0.614349 | 0.571614 |
| q median | 0.733067 | 0.792875 |
| q p90 | 0.897662 | 0.929516 |
| q p99 | 0.967437 | 0.969782 |
| q max | 0.999398 | 0.997678 |
| Matched FEA accuracy | Unmeasured | Unmeasured |

The candidate is 3.623612x faster with 81.796% fewer tets and 79.333% fewer
nodes. Minimum and lower-tail quality are worse, boundary error is higher,
and maximum quality is slightly worse. Median and upper quantiles improve.
Both pass the distance 0.002, absolute relative volume 1%, q 0.1, sampled
positive curved-Jacobian, facet/conformity and count gates. This bounded
linear-density grade is nonperiodic for both routes; periodic fitting is not
implemented or established by this result. The result supports a targeted
curved-fitting direction, not a universally better mesh.

## Verification and limits

Public feature commit `dfb6fe9` includes 13 interpolation checks covering values,
knots, endpoints, gradients, snapshots, unsupported cases and real generation.
The full Python/API and tools run passes 208 tests with one skipped. Updated
interpolation and Schwarz D controls pass 15 tests. Seven prior admitted
curved meshes, including the original four specimens, regenerate exactly equal
point, tet and Tri6 arrays. The native binary is unchanged at
`6394fdcd5e13c42955f96625c2e91800e5b4c2cb2817fb3349de3366dbcfab0a`.

The automatic compiler leaves the original Python callable in the caller's
variable. The initial generation observer therefore checked that callable
against itself. The worker now obtains the owned compiled field for its check.
The independent automatic-study replay compensates for the old observer with
8,963 compiled-field probes, maximum difference 4.44e-16 and analytic gradients
present. Old timing records remain unchanged.

Failed audit job 665 contains a verifier indentation error; corrected job 666
passes. Expected red feature job 669 has six failures before implementation;
green job 670 passes all 13 checks. The coarse uniform control in job 677 stays
rejected. Saved fitting job 679 is kept separate from live pilot 680, pilot
replay 681, repeated workflow 682 and repeated replay 683. Changed production
code passes lint and formatting. The recorded 973 legacy lint findings in five
frozen handoff scripts are not repaired as part of this study.

All fixed tables, specimen samples and 24 retained lattice model tensors keep
their hashes. No new FE solve runs and no reference is promoted. Best saved
native lattice disagreement remains 1.213247%, versus 0.994671% for the separate
modified actual star-refined route. Gyroid and spinodoid still lack admitted
actual comparators in the tested routes. Schwarz P remains a supported strong
generation win. The original four-case objective remains open.

The next scientific priority is a lattice error-source diagnostic that can
distinguish fixed-domain discretization from geometry/reference discrepancy.
The cap-fitting result merits broader geometry and periodic controls before
any native default or matched-stiffness claim. Another global sizing sweep is
not justified by these results.
