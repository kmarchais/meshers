# Public fedoo integration on G3

`Mesh.to_fedoo()` now copies Tet4 or curved Tet10 coordinates and connectivity
into an independent fedoo mesh. The direct route passes the integration checks
on fedoo 1.0.0b2. Microgen 2.0.0b1 still rejects quadratic tetrahedra in both
`SingleMesh.from_pyvista` and `BoxMesh.from_pyvista`.

Job 9 ran the conversion tests, two solid-cell patches and the saved admitted
12-cell curved Gyroid and Schwarz P solid meshes. It used one numerical thread,
an 8 GiB worker ceiling and a 600-second timeout per stage. Peak sampled worker
RSS was 0.851 GiB. The progress audit passed. All frozen step-1 inputs retained
their hashes.

## Results

The integration uses six prescribed engineering strains and independently
integrates their stresses and mixed-load energies. The specimen comparison
uses the saved stress-controlled tensors from the Microgen workflow study.
The cube reference is the analytic isotropic material tensor, E = 2.5 and
nu = 0.25.

| Case | Nodes | Tetrahedra | Maximum relative tensor disagreement | Periodic mid-edge pairs, x/y/z |
|---|---:|---:|---:|---|
| Solid Tet4 | 189 | 768 | 1.33e-15 | 0 / 0 / 0 |
| Solid Tet10 | 1,241 | 768 | 4.44e-16 | 56 / 56 / 56 |
| Gyroid Tet10 | 19,931 | 10,832 | 4.05e-12 | 236 / 248 / 236 |
| Schwarz solid Tet10 | 19,628 | 11,976 | 1.20e-11 | 420 / 420 / 420 |

All periodic displacement jumps, including mid-edge jumps, are below 6e-17.
Periodic equivalence-class force imbalance is below 1.1e-10 relative.
Integrated stress agrees with macroscopic reactions within 1.1e-11 relative.
Mixed-load Hill–Mandel disagreement is below 4.9e-12. All 24 reduced solves
have relative residuals below 1e-12.

VTK and fedoo shape functions agree within 3.4e-16 at corresponding reference
points. Their node order agrees, including all six edge nodes. Their reference
coordinate conventions differ: a fedoo point `(xi, eta, zeta)` maps to the VTK
point `(zeta, 1-xi-eta-zeta, xi)`. Job 8 first compared the same numerical
coordinates in both conventions and failed the solid-patch audit. The test
mapping was corrected, and job 9 used a fresh run directory. Public conversion
code did not change in response to this test failure.

The experiment verifies unique edge-node ownership, matches full six-node
periodic cap faces and checks positive Jacobians at quadrature points, reference
nodes and the centroid. The smallest sampled Jacobians are 2.28e-5 for Gyroid
and 1.85e-5 for Schwarz. Sampling is not a proof of positivity everywhere.

## Reproduce

Use the existing G3 environment and local saved mesh archives from step 1.
From the repository root, select a new output directory:

```sh
qrun -l 'Tet10 public fedoo integration' --progress report \
  --interrupt none --interrupt-reason 'Immutable staged experiment' \
  --check-start --wait env \
  LD_LIBRARY_PATH=/home/kmarchais/dev/meshers/.local-work/g3-runtime/root/usr/lib/x86_64-linux-gnu \
  OPENBLAS_NUM_THREADS=1 timeout 2400 \
  .venv/bin/python experiments/adaptive/tet10_integration.py \
  --run experiments/adaptive/runs/g3-tet10-integration-next
```

The runner records source and input hashes, package versions, worker limits,
full tensors and all six displacement fields. Text records are committed under
`experiments/adaptive/runs/g3-tet10-integration-v2`. Field archives stay local.
The three focused public conversion tests passed. Ruff and formatting checks
passed for the public method and tests.

## Limits and next step

Conversion does not create boundary conditions, transfer boundary tags or
periodic transforms, or validate an arbitrary input mesh. Rotational periodicity
is outside this result. The native mesher and the step-1 admission gates are
unchanged. These tests establish integration and tensor reproduction; local
stress accuracy and continuum error remain unproven.

## Stronger genuine Microgen control, 18 October

Job 152 lifts the admitted actual Microgen VTK/MMG Schwarz mesh to straight
Tet10. Every added node is the exact midpoint of one unique original edge.
The original planar boundary and domain are unchanged. The parent uses the
explicit one-pass override from the matched Schwarz study, not Microgen's
default two-pass route. Its source and MMG-call hashes are preserved.

The lift has 230,337 nodes and 154,012 elements, within the explicitly staged
250,000-node solve limit and unchanged 8 GiB worker ceiling. Independent
geometry and periodic checks pass. Arrays returned by `Mesh.to_fedoo()` feed
all six prescribed-strain solves; no native Meshers generation is called.
This is manual elevation plus direct fedoo, not support added to Microgen's
public Tet10 converters.

Four-point integration is checked against fifteen-point fixed-field energy,
which differs by less than 9.57e-11 relative. Independent reaction disagreement
is 3.27e-11, Hill-Mandel disagreement is 2.81e-11 and periodic mid-edge jumps
are below 5.56e-17. All six reduced solves meet the strict residual check.

Worst mixed-load stiffness disagreement is 0.062583% against the primary
Schwarz reference, and 0.054006% to 0.062702% across all declared fine-window
orders and length proxies. The original linear Microgen parent spans
0.92778% to 0.99022%. Tet10 materially strengthens the competitor.

The solve takes 1,128.92 seconds; the complete worker uses 1,152.73 seconds
and peaks at 6.05 GiB. These are saved-mesh diagnostic costs, not a live
Microgen extraction/MMG/conversion/solve workflow. Do not use them to claim
a new matched-workflow speed ratio.

The stronger control also exposes a Meshers limit. Curved12 reaches 0.13559%
across the same reference models, so it does not meet a stricter 0.1% target
that this Microgen control passes. Select one current-native curved14 pilot
against that same target. The original four-fixture 1% objective remains open.
Text records are in `runs/g3-schwarz-microgen-straight-tet10-v1` and
`runs/g3-schwarz-microgen-straight-tet10-reference-v1.json`.

## Supplementary 0.1% target

The stricter target supplements the original 1% objective; it does not replace
that objective or settle any other fixture. Curved14 and curved15 both pass
geometry and independent field checks, but their model envelopes reach
0.101239% and 0.103035%. Both strictly miss 0.1%, despite passing 1%.
Do not round either into a stricter-target pass. Seven comparison-admission
and target regression cases pass, including the genuine Microgen positive
control and the curved14 negative control.

The existing curved16 reference data indicates a passing 0.093778% envelope.
Job 157 refreshes that known control with the current verified native build
and the prescribed-strain solver. It has 42,928 nodes and 27,492 elements,
H=0.000162, and an 88.01-second live workflow. Its tensor reproduces the
archived sixteen-cell input within 2.67e-11 relative. The declared envelope is
0.031469% to 0.093778%. Job 158 independently checks its fields with reaction
disagreement 1.18e-11 and periodic jumps below 5.56e-17. This also checks the
numerical reproduction of one input to the reference models, not independent
continuum convergence.

No new strict-target workflow speed ratio is claimed against job 152's
saved-mesh solve. Before a live comparison, job 159 tests whether that
Microgen Tet10 parent overresolves its interior. Change only first-pass hmax
from 0.05 to 0.1, preserving the VTK30 seed, hausd=0.0005, hgrad=1.3, explicit
one-pass route and independent admission gates. Audit the actual returned mesh
before selecting any Tet10 lift and solve.

Job 159 passes admission with 29,204 Tet4 nodes, 136,587 elements,
H=0.0009615 and q_min=0.10387. Its exact lift has 207,200 Tet10 nodes,
only 10.0% fewer than the original 230,337-node control. Job 160 runs that
coarser straight-Tet10 control with the same six-load solver and independent
audits. No accuracy or timing conclusion is assigned before reading it.

Job 160 completes the coarser control. Its declared model range is 0.053985%
to 0.062716%, so it still passes the supplementary 0.1% target. The solve
takes 1,026.88 seconds; the worker takes 1,049.97 seconds and peaks at 5.40 GiB.
Reaction disagreement is 1.47e-11, Hill-Mandel disagreement is 1.09e-11 and
periodic jumps remain below 5.56e-17. This modest coarsening reduces nodes
without losing the declared target. It remains a saved-mesh control, not a
paired live workflow measurement. Return to the unresolved lattice recovery
validation rather than spending further runs tightening Schwarz alone.

The next requested priority is a coarse FEA stress-recovery estimate for linear
mesh sizing. The new displacement archives provide integration evidence, but
they are not qualified local-stress references for that study.
