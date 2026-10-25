# Lattice reference assumptions and measured edge refinement

## Result, 25 October 2026

The measured-edge Tet10 refinement improves worst disagreement across all 24
retained models from 1.283822914% to 1.216217882%. It still fails the unchanged
1% target. It is not promoted to a production default, native workflow winner
or stiffness reference.

The experiment retains the current native parent's exact curved domain and
complete quadratic displacement space. It selects specific periodic edge-star
bisections by positive, block-diagonal cubic residual potential per added-node
upper bound. This differs from the historical longest-edge, curve-only and
cell-score rankings. The potential is a selection indicator, not an error bound
or a prediction of the h-refinement gain. It creates ordinary Tet10 elements;
it does not label the enriched cubic field as a native mesh.

Positive 64/125-point indicator vectors differ by 0.0122661% relative, within
the existing 1% precision gate. Edge features carry 38.5932% of the block-diagonal
potential; the remaining face-feature potential is reported, not discarded
from an asserted error bound. Shared periodic edge scores are distributed
across physical copies so their orbit sums reproduce the measured potential.

Preflight 372 admits 249,998 nodes and 163,802 tets under the unchanged caps.
Minimum corner quality is 0.113261924. Sampled distance remains
0.000343162247760018, and relative volume changes by only -1.776e-15.
The original curved maps and all six fields reproduce to 3.331e-16 over
12 probes in every child. Periodic nodes and cap facets match, and sampled
curved Jacobians stay positive. Actual sparse setup has 726,429 free rows,
58,586,609 stored entries and four AMG levels. Setup peaks at 5.372 GiB under
8 GiB and stops before CG. Preflight takes 109.984 seconds.

Job 374 completes six fresh 15-point AMG loads in 206.939 seconds, peaking at
5.772 GiB. Independent field and positive-quadrature audit takes 54.977 seconds.
The assessment takes 67.342 seconds and verifies a Galerkin energy balance
within 1.142e-9. Every stiffness correction from the original parent is negative,
with the largest 0.068247586%. All 24 models remain unchanged and scored.
The preflight, parent solve and indicator computation are not free costs;
this saved-domain investigation makes no comparative workflow timing claim.

The new selector has a frozen red control 370 and passing nested/hierarchical
controls 371, 17 tests. Solve-driver provenance controls 373 accept either
recorded preflight hash format and reject missing, conflicting or modified
parent evidence. No frozen manifest is rewritten to fit the newer format.

## Reference replay

Job 369 verifies every available original input hash and reproduces all 24
extrapolated matrices from the frozen coarse/fine tensors and declared length
ratios and orders. No required or optional listed provenance is missing.
The fine-result coefficient ranges up to 31.5900 in the closely spaced adaptive
pair. The current native worst model is the uniform-seed first-order model,
whose fine coefficient is 8. These are distinct conditioning observations.

Each original coarse/fine pair has a nonzero physical-volume change above
roundoff. Exact common domains, nested spaces and observed convergence orders
are not certified by those reports. Their stated role remains empirical
refinement sensitivity. No tensor is removed, target relaxed, new model fitted
or continuum reference promoted on this evidence.

## Scope and next action

Schwarz retains the recorded repeated matched-target win against an explicitly
modified actual Microgen route. Gyroid and spinodoid lack admitted actual
comparators in the tested routes. Lattice accuracy and reference convergence
remain open. All four fixed specimens and every geometry, quality, periodic
and resource gate remain in force.

Do not repeat this h-refinement recipe automatically. A subsequent native
change must reproduce a specific evaluation or topology defect first. A
physical-side periodic-gradient regression is being evaluated separately;
these edge-refinement results do not establish its outcome.

Receipts:

- [Reference replay](../../handoff/macbook/lattice-reference-assumptions-20261025.json)
- [Geometry and resource preflight](../../handoff/macbook/lattice26-edge-residual-preflight-20261025.json)
- [Six-load assessment](../../handoff/macbook/lattice26-edge-residual-assessment-20261025.json)
