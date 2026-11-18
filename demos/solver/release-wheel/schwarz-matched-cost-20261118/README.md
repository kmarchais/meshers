# Schwarz P workflow cost at 1% stiffness discrepancy

The resolution-30 Microgen Tet4 repair and the finest native Tet4 setting both pass the sampled geometry audit and the 1% stiffness target. The error is the worst relative strain-energy discrepancy across mixed loads against the frozen quadratic stiffness tensor. It measures agreement with that reference.

## First-pass MMG sizing

The screen varies one parameter on the same resolution-30 VTK seed. The single MMG pass uses `nofem`, `hausd=0.0005` and `hgrad=1.3`. The interior repair keeps boundary coordinates fixed. Each column gives stiffness discrepancy in percent after the geometry check. A rejected mesh skips the solver.

| CPU and system | hmax 0.050 | 0.040 | 0.035 | 0.030 |
| --- | --- | --- | --- | --- |
| [Apple M4 · macOS](apple-m4/sizing.json) | 0.927% | 0.799% | Quality 0.09628 | 0.750% |
| Core i7-10750H · Linux | Pending | — | — | — |
| Core i7-12700H · Windows | Pending | — | — | — |
| Core Ultra 9 285K · Windows | Pending | — | — | — |
| Ryzen 7 7730U · Linux | Pending | — | — | — |


## Fastest observed qualifying settings

Each timing is the median of three fresh processes on the same configuration. The native arm uses 24 background cells, maximum edge length 0.032 and adaptation tolerance 0.0003. The Microgen arm is the lowest complete median across qualifying settings in the screen. Each repeat must pass the geometry, stiffness and twelve solver checks. The order rotates by one setting each round. With four or five arms, three rounds only partly balance execution order.

| CPU and system | Route | Wall (s) | Build, repair and I/O (s) | Geometry audit (s) | Assembly (s) | CG, 12 loads (s) | Peak RSS (MiB) | Free DOFs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| [Apple M4 · macOS](apple-m4/sizing.json) | Native | 32.0 | 2.5 | 0.6 | 0.6 | 25.3 | 1659 | 115,779 |
| [Apple M4 · macOS](apple-m4/sizing.json) | Microgen hmax=0.050 | 29.2 | 3.8 | 0.4 | 0.5 | 21.8 | 1518 | 89,361 |
| Core i7-10750H · Linux | Pending | — | — | — | — | — | — | — |
| Core i7-12700H · Windows | Pending | — | — | — | — | — | — | — |
| Core Ultra 9 285K · Windows | Pending | — | — | — | — | — | — | — |
| Ryzen 7 7730U · Linux | Pending | — | — | — | — | — | — | — |


Wall time includes process startup, generation, repair where applicable, export, audit, assembly, linear solves and solver cross-checks. The stage columns show selected parts of that total. Peak RSS is the median of the per-run peaks, sampled across the worker and MMG children every 0.25 seconds. The screen is a one-time search cost and is recorded separately. [All settings and median costs](summary.json)

In a separate replay at `hmax=0.05`, time inside the fixed-boundary repair was [0.42 s on Apple M4](apple-m4/repair.json).

## Held-out Schwarz phase

The twelve CG solves accounted for roughly three quarters of the earlier native finest workflow. We tested symmetric reverse Cuthill–McKee ordering to improve sparse matrix access. The experimental solver builds a new ordering on every reduced system, executes the same Jacobi-preconditioned CG, returns the displacement in original order and checks the original-system residual. The benchmark charges all of this work.

Before the experiment, we fixed the x phase at +1/8 period, native mesh settings, CG tolerances and three paired runs in AB, BA, AB order. A periodic translation preserves the continuum homogenized tensor while moving the wall intersections with the cube faces. An independent trigonometric expression checks the translated geometry. Both arms must pass 1% and reproduce identical point and cell arrays within a configuration before we compare solver cost.

| CPU and system | Baseline wall (s) | Reordered wall (s) | Wall change | CG before / after (s) | Ordering (s) | Stiffness discrepancy |
| --- | --- | --- | --- | --- | --- | --- |
| [Apple M4 · macOS](apple-m4/phase.json) | 32.5 | 31.2 | -4.0% | 25.8 / 24.0 | 0.6 | 0.712% |
| Core i7-10750H · Linux | Pending | — | — | — | — | — |
| Core i7-12700H · Windows | Pending | — | — | — | — | — |
| Core Ultra 9 285K · Windows | Pending | — | — | — | — | — |
| Ryzen 7 7730U · Linux | Pending | — | — | — | — | — |


Three repeats on one held-out specimen measure the effect on these configurations. CPU, operating system, architecture, sparse libraries and memory vary together. These measurements cannot isolate an operating-system effect. The experimental ordering is opt-in.

## Replay

Install the [pinned solver environment](../schwarz-series.md) and run each command with a new output directory. Both drivers request one numerical thread per library. Each worker has a 750-second budget. The reproducible driver writes screen results, individual timings, geometry checks, all linear solves and sampled process-tree RSS.

```sh
python demos/solver/schwarz_trials.py sizing study-sizing
python demos/solver/schwarz_trials.py phase study-phase
```

The raw archives include the hardware snapshot, UTC acquisition times, package versions, input hashes and driver hashes from the run. `solver_mesh_sha256` in the phase trial covers the ordered solver coordinates and connectivity. It allows direct comparison across the two solver arms without relying on HDF container metadata.
