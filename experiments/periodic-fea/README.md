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
Each worker has a default 20-minute timeout and 16 GiB RSS limit. A prepared
fine-reference run can declare a longer timeout through `--worker-timeout`.
`--max-tetrahedra` declares a larger native element budget for a finer reference;
its default remains two million and it does not change any admission tolerance.
Each generated mesh
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

`--tasks` selects the next configurations as `workflow:cells` pairs. An optional
`case/workflow:cells` entry selects an analytic solid control. `mmg-native-seed`
measures single standard MMG adaptation of a valid one-pass Meshers seed. It uses
Microgen's required box-cap treatment and native CLI. It omits the public
pipeline's `nofem` and level-set passes. Its total cost includes native seed
generation. This algorithm control is separate from the current Microgen
adoption baseline. A failed solid-cell control stops a batch before its
expensive references run.

`periodic-schwarz-p/workflow:cells` selects a separate periodic Schwarz P pilot
with raw band [-0.6, 0.6]. It needs its own fine stiffness reference; the Gyroid
tensor cannot score this specimen. Analytic geometry probes use the selected
field. Default task selection still uses the original solid/Gyroid cases.

For `mmg-native-seed`, `--mmg-seed-cells` can set the native seed resolution
independently of the task's MMG sizing resolution. The default uses the same
resolution for both. Results record `native_seed_cells` and `mmg_sizing_cells`.
This supports testing volume coarsening while retaining a finer periodic seed
boundary. All final admission and physics checks remain active.

`mmg-native-seed-metric` uses a frozen `--mmg-size-map` instead of uniform MMG
sizing. A paired `mmg-native-seed` task in the same prepared run retains uniform
sizing. The source map's coordinates, tetrahedra, base size and file hash are
checked, then scalar sizes are matched to MEDIT vertex order. Such a map may
require prior FEA; its generation timing excludes that earlier solve and must
not be reported as a cold one-shot pipeline cost.

`gmsh-hxt-native-surface` is an experimental volume-construction control. It
extracts the full boundary from a native one-pass seed at the task resolution,
then runs Gmsh HXT with uniform edge sizing. It checks that every boundary
vertex and triangle is preserved before the common admission and FEA audits.
Preparation records the installed Gmsh adapter/library hashes and version.
Run the analytic solid control first. Generation cost includes the full native
volume seed and boundary checks; it is not a surface-only timing estimate.
`--surface-seed-cells` separates native boundary resolution from HXT interior
sizing. The analytic solid control always uses its task resolution for the seed.

Compare completed runs against the retained grid-64 tensor with:

```powershell
benchmarks/microgen/.venv/Scripts/python.exe experiments/periodic-fea/compare_reference.py experiments/periodic-fea/runs/reference-and-seed-v1/02-periodic-gyroid-meshers-pass1-64/result.json experiments/periodic-fea/runs/refinement-v1 experiments/periodic-fea/runs/reference-and-seed-v1
```

The comparison records the input artifact hashes, all six relative energy
eigenvalues and admission status. It rejects incomplete runs. A passed geometry
or physics audit alone does not establish continuum accuracy or repeatable speed.

`repeat_generation.py` reuses the frozen generation code for one warmup and
three timing samples of the native grid-40/MMG seed grid-32 pair. Every sample
must reproduce the corresponding audited canonical mesh exactly. It repeats
generation only and does not claim repeated solver measurements.
That pair was selected against the historical grid-64 reference. The completed
grid-80 comparison in [RESULTS.md](RESULTS.md#completed-grid-80-reference-check)
supersedes its accuracy qualification: native grid 40 no longer passes the
strict 1% worst-energy target against the finest measured tensor.
