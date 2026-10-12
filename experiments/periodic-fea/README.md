# Direct periodic meshing for FEA

This experiment compares Meshers' direct periodic generator with Microgen's
public VTK-to-MMG pipeline using `periodic=True`. The MMG control uses native CLI
executables. It does not use the mmgpy Python bindings.

The first batch checks a solid cube against its analytic elasticity tensor,
then tests a periodic Gyroid sheet at grids 16 and 24. Its sheet band is
`[-0.55, 0.55]` in the raw period-one Gyroid field, matching the existing
fedoo homogenization example. Microgen's offset is 1.1. Both implementations
use the same centered unit cube and analytic material set.

## Protocol

- Test direct native generation with zero, one and four optimization passes at
  grid 16. Test one pass at grid 24.
- Use the supported Microgen two-pass MMG pipeline at both grids, with periodic
  boundary preservation, `hausd=0.002`, `hsiz=1/cells` and `hgrad=1.3`.
- Keep the common geometry gates: sampled boundary distance in both directions
  at most 0.002 of the cell size, volume error plus empirical sampling uncertainty
  at most 1%, valid nonoverlapping tetrahedra, and minimum MMG quality 0.1.
- Independently require bijective opposite-face node matching and identical face
  triangulations. Combine edge and corner pairs into equivalence classes.
- Solve six periodic elasticity loads with the same fedoo assembly and solver.
  Audit displacement jumps, nodal-force equilibrium, Hill–Mandel consistency,
  reactions, symmetry, positive stiffness and the uniform-strain upper bound.
  Cross-check fedoo's strain- and stress-controlled homogenization paths.
- Preserve meshes that fail geometry or shape admission. Their stiffness may
  be computed as a diagnostic when volume validity and periodicity pass.
  Diagnostics do not qualify for accuracy or speed claims.

The homogenization code starts from `demos/solver/fedoo_homogenization.py`.
The experiment adds progress reporting, solver timing and saved displacement
fields without changing its material or solver checks. The first batch uses
single timing samples without warmup. See [measured results](RESULTS.md).

The generalized eigenvalue score compares full stiffness tensors over all
macroscopic strain combinations. A difference between two refinement levels
measures stability only. It does not certify continuum error. A fine reference
with observed convergence and geometry checks is required before admitting
2% or 1% FEA targets.

## Run the first batch

The isolated environment uses the benchmark's locked dependencies, fedoo
1.0.0b2, and the previously verified packaged PR #5 Meshers build. The prepared
manifest records dependency versions, Python source hashes and the actual native
binary hash. The prepared snapshot includes that native build's Python package
and extension. Workers import this frozen package ahead of the installed released
baseline. Each worker verifies its native binary and dependency identities
before generating a mesh.

```powershell
$candidatePackage = "C:/Users/kevin/AppData/Local/uv/cache/archive-v0/5mCyN6F85EkLo3NI/Lib/site-packages/meshers"
benchmarks/microgen/.venv/Scripts/python.exe experiments/periodic-fea/run.py --prepare experiments/periodic-fea/runs/coarse-v2 --native-package $candidatePackage
qrun -g gpu -l "Periodic Gyroid: direct Meshers vs Microgen MMG pilot" --progress report benchmarks/microgen/.venv/Scripts/python.exe experiments/periodic-fea/runs/coarse-v2/source/run.py --run experiments/periodic-fea/runs/coarse-v2
```

The job first runs the pairing and tensor-score regression checks. It then runs
one bounded worker at a time and reports a fraction, phase, ETA and live values.
Each worker has a 20-minute timeout and 16 GiB RSS limit. Each generated mesh
receives exhaustive overlap validation up to the declared 60-million-pair limit.
Reaching a validation or resource limit is a recorded failure.

`results.json` contains completed results. Each worker folder retains its log,
request, raw and canonical meshes, and periodic pair metadata. MMG folders also
retain every pass's input/output mesh and exact parameters. Archival capture
overhead is measured separately and subtracted from MMG generation samples.
Native public diagnostics and analytic sampling have their own recorded scope.
Follow-up batches also save all six displacement fields. Their stiffness, energy
and periodic force balance are reconstructed independently from those artifacts.

## Choose the next experiment from the result

First inspect whether the baseline's nodes and face triangulations stay periodic
and whether either pipeline passes the common geometry and shape gates. Then
compare the six-load-case tensors and identify the cheapest native profile.
Only after that comparison should the next finer reference or parameter trial
be queued. A lower minimum-quality threshold or a wider boundary tolerance
must never silently turn a failed configuration into a passing one.
