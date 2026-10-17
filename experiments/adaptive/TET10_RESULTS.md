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

The next requested priority is a coarse FEA stress-recovery estimate for linear
mesh sizing. The new displacement archives provide integration evidence, but
they are not qualified local-stress references for that study.
