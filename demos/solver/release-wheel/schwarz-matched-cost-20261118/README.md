# Schwarz P workflow cost at 1% stiffness discrepancy

Both the native mesh and the resolution-30 Microgen repair pass the sampled geometry checks and the 1% effective-stiffness threshold. The stiffness metric is the worst relative strain-energy discrepancy over mixed loads against the frozen quadratic reference.

## Fastest observed qualifying settings

Native Tet4 had the lower median wall time on 2 of 2 configurations with a finished comparison. It saved 32.3–33.5% against the fastest repaired setting. Native stiffness discrepancy was 0.992% and repaired Microgen measured 0.927–0.928%.

| CPU and system | Native wall (s) | Repaired wall (s) | Native / repaired Tet4 | Native / repaired nodes | Paired native time changes |
| --- | --- | --- | --- | --- | --- |
| [Apple M4 · macOS](apple-m4/bracket.json) | 19.7 | 29.1 | 133,730 / 154,002 | 30,070 / 32,099 | -32.9%, -32.2%, -32.3% |
| Core i7-10750H · Linux | Pending | — | — | — | — |
| Core i7-12700H · Windows | Pending | — | — | — | — |
| Core Ultra 9 285K · Windows | Pending | — | — | — | — |
| [Ryzen 7 7730U · Linux](ryzen-7-7730u/bracket.json) | 34.6 | 52.1 | 133,730 / 154,012 | 30,070 / 32,070 | -33.5%, -33.1%, -34.3% |

The first timing round compared the native finest setting (maximum edge 0.032, adaptation tolerance 0.00030) with every repaired setting that passed its screen. The fastest repaired setting beat that native mesh, which reached 0.706% stiffness discrepancy. We then tested one intermediate native setting (0.038, 0.00041) between finest and the earlier failing fine setting (0.040, 0.00045 at 1.032%). This second round used the same geometry, loads, reference, mesh budget and 1% limit.

The paired round alternates native and repaired `hmax=0.05` in AB, BA, AB order. Values above are medians of three fresh processes per setting. The three changes compare individual pairs with the repaired time as the denominator. Screening costs and the initial finest comparison remain in the raw records. [All settings and medians](summary.json)

## First-pass MMG sizing

The screen varies one parameter on the same resolution-30 VTK seed. The single MMG pass uses `nofem`, `hausd=0.0005` and `hgrad=1.3`. The interior repair keeps boundary coordinates fixed. Each column gives stiffness discrepancy in percent after the geometry check. A rejected mesh skips the solver.

| CPU and system | hmax 0.050 | 0.040 | 0.035 | 0.030 |
| --- | --- | --- | --- | --- |
| [Apple M4 · macOS](apple-m4/sizing.json) | 0.927% | 0.799% | Quality 0.09628 | 0.750% |
| [Core i7-10750H · Linux](core-i7-10750h/sizing.json) | 0.928% | 0.798% | 0.787% | 0.749% |
| Core i7-12700H · Windows | Pending | — | — | — |
| Core Ultra 9 285K · Windows | Pending | — | — | — |
| [Ryzen 7 7730U · Linux](ryzen-7-7730u/sizing.json) | 0.928% | 0.798% | 0.787% | 0.749% |


## Cost breakdown

Wall time starts before process creation and ends after the geometry audit, six strain-controlled solves, six stress-controlled solves and mechanics cross-checks. Build time includes generation, export, reading the volume back and any repair. The other stage columns show selected parts of the total. Memory is the median of the process-tree peak RSS samples from each run, polled every 0.25 seconds.

| CPU and system | Selected route | Build incl. repair (s) | Geometry audit (s) | Assembly (s) | CG, 12 loads (s) | Peak RSS (MiB) | Free DOFs |
| --- | --- | --- | --- | --- | --- | --- | --- |
| [Apple M4 · macOS](apple-m4/bracket.json) | Native 0.038 | 1.6 | 0.4 | 0.4 | 14.9 | 1279 | 80,412 |
| [Apple M4 · macOS](apple-m4/bracket.json) | Repaired hmax=0.050 | 3.8 | 0.4 | 0.5 | 21.8 | 1546 | 89,361 |
| Core i7-10750H · Linux | Pending | — | — | — | — | — | — |
| Core i7-12700H · Windows | Pending | — | — | — | — | — | — |
| Core Ultra 9 285K · Windows | Pending | — | — | — | — | — | — |
| [Ryzen 7 7730U · Linux](ryzen-7-7730u/bracket.json) | Native 0.038 | 3.4 | 0.8 | 0.6 | 25.3 | 964 | 80,412 |
| [Ryzen 7 7730U · Linux](ryzen-7-7730u/bracket.json) | Repaired hmax=0.050 | 8.4 | 0.9 | 0.8 | 37.7 | 1005 | 89,274 |

A separate replay timed the fixed-boundary repair at `hmax=0.05`. The measured times were [0.42 s on Apple M4](apple-m4/repair.json), [0.71 s on Ryzen 7 7730U](ryzen-7-7730u/repair.json).

## Held-out Schwarz phase

The twelve CG solves accounted for roughly three quarters of the earlier native finest workflow. We tested symmetric reverse Cuthill–McKee ordering to improve sparse matrix access. The experimental solver builds a new ordering on every reduced system, executes the same Jacobi-preconditioned CG, returns the displacement in original order and checks the original-system residual. The benchmark charges all of this work.

Before the experiment, we fixed the x phase at +1/8 period, the finest native mesh settings (24 background cells, maximum edge 0.032, adaptation tolerance 0.00030), CG tolerances and three paired runs in AB, BA, AB order. A periodic translation preserves the continuum homogenized tensor while moving the wall intersections with the cube faces. An independent trigonometric expression checks the translated geometry. Both arms must pass 1% and reproduce identical point and cell arrays within a configuration before we compare solver cost. The largest stiffness-entry difference must also stay within one billionth of the largest baseline entry.

| CPU and system | Baseline wall (s) | Reordered wall (s) | Wall change | Paired changes | CG before / after (s) | Ordering (s) | Stiffness discrepancy |
| --- | --- | --- | --- | --- | --- | --- | --- |
| [Apple M4 · macOS](apple-m4/phase.json) | 32.5 | 31.2 | -4.0% | -4.1%, -3.9%, -4.6% | 25.8 / 24.0 | 0.6 | 0.712% |
| [Core i7-10750H · Linux](core-i7-10750h/phase.json) | 91.0 | 85.6 | -6.0% | -6.0%, -7.0%, -6.6% | 75.2 / 68.5 | 1.4 | 0.712% |
| Core i7-12700H · Windows | Pending | — | — | — | — | — | — |
| Core Ultra 9 285K · Windows | Pending | — | — | — | — | — | — |
| [Ryzen 7 7730U · Linux](ryzen-7-7730u/phase.json) | 57.9 | 53.2 | -8.2% | -8.3%, -8.3%, -7.8% | 44.6 / 38.9 | 1.1 | 0.712% |


Across 3 completed configurations, median end-to-end time changed by -8.2% to -4.0%. The largest absolute stiffness-entry difference across the six runs on any one configuration was 3.9e-13.

Three repeats on one held-out specimen measure the effect on these configurations. CPU, operating system, architecture, sparse libraries and memory vary together. These measurements cannot isolate an operating-system effect. The paired changes show the spread of the three measurements. The experimental ordering is opt-in.

## Replay

Install the [pinned solver environment](../schwarz-series.md) and run each command with a new output directory. All three studies request one numerical thread per library. Each worker has a 750-second budget. The reproducible driver writes screen results, individual timings, geometry checks, all linear solves and sampled process-tree RSS.

```sh
python demos/solver/schwarz_trials.py sizing study-sizing
python demos/solver/schwarz_trials.py bracket study-bracket
python demos/solver/schwarz_trials.py phase study-phase
```

The raw archives include the hardware snapshot, UTC acquisition times, package versions, input hashes and driver hashes from the run. `solver_mesh_sha256` in the phase trial covers the ordered solver coordinates and connectivity. It allows direct comparison across the two solver arms without relying on HDF container metadata.
