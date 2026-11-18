# Where the Schwarz P comparison misses tolerance

On all five configurations the first MMG pass moved the Schwarz P wall beyond the frozen distance limit. The historical single-pass recipe kept the wall within tolerance. Tightening the native Tet4 sizing also produced a mesh below the 1% stiffness target on every configuration.

## Microgen stages

These distances are the largest sampled distances from analytic Schwarz P wall points to the nearest mesh wall triangle, in a unit cell. The range spans configurations. Each run also records the forward distance bound, cap and component checks, and the physical coordinates of the worst sample. The distance limit is 0.002.

| Stage | Default, resolution 24 | Tuned, resolution 24 | Historical control, resolution 30 |
| --- | ---: | ---: | ---: |
| VTK seed | 0.001370 | 0.001370 | 0.000756 |
| First MMG pass | 0.005609–0.007024 | 0.005609–0.007024 | 0.000707–0.000764 |
| Export | 0.005610–0.008245 | 0.004356–0.005813 | 0.000707–0.000764 |

The two resolution-24 routes share an unmodified first pass. On Windows, the default export's worst analytic sample was at (0.206935, 0.143794, 0.423323). At that same point, the distances were 0.000082 on the seed, 0.004870 after pass one and 0.008245 after pass two. Worst locations vary across configurations. The seed passed the wall check but failed other checks, including signed orientation.

The historical control uses one MMG pass with `nofem`, `hausd=0.0005`, `hmax=0.1` and `hgrad=1.3`. The fixed-boundary interior repair raised the low element qualities and all five repaired volumes passed admission. Resolution and pass settings both changed from the default, so this control does not identify the effect of any individual setting.

## Native stiffness bracket

The Tet4 series used the same 24-cell background, material, periodic constraints and six-load reference as the [original series](../schwarz-series.md). The columns below give maximum edge length, adaptation tolerance, free displacement DOFs, worst mixed-load energy discrepancy, and local σxx at (−0.364, −0.347, −0.313) under unit xx strain.

| Case | Maximum length | Adaptation tolerance | Free DOFs | Energy discrepancy | Probe σxx |
| --- | ---: | ---: | ---: | ---: | ---: |
| Medium | 0.060 | 0.00070 | 37,689 | 1.574% | 1.73568 |
| Fine | 0.040 | 0.00045 | 74,085 | 1.032% | 1.63293 |
| Finest | 0.032 | 0.00030 | 115,779–115,782 | 0.706% | 1.63211 |

The finest mesh used 189,871–189,873 tetrahedra and 42,319–42,320 nodes, within both 250,000 limits. Every Tet4 mesh passed the independent geometry audit. The similar local stresses at the two finer levels measure sensitivity to refinement. A reference stress at the probe would be needed to assign a local error.

A separate curved Tet10 control at 12 background cells had 11,976 elements, 19,628 nodes, 53,853 free DOFs and 0.0738% stiffness discrepancy on all five configurations. The control checked curved Jacobian positivity, integrated volume and six reduced-system residuals. A full curved-wall and cap audit is still needed to admit Tet10 to the workflow timing comparison.

## Can both Tet4 routes be timed at 1%?

Four configurations completed the additional repaired Microgen screen. At resolution 30 the meshes passed geometry and missed the stiffness target, with discrepancies of 1.168% to 1.177%. At resolution 34 the first MMG pass failed on all five attempts per configuration. At resolution 36 the exported wall passed, but 30–32 tetrahedra still fell below the 0.1 quality floor after repair. The study preserved each rejection and skipped the paired timings. The Core Ultra 9 285K follow-up is queued.

For a qualifying pair, the protocol uses three fresh processes per route in alternating order. It measures the entire workflow, including generation, repair, export, geometry admission, assembly and solve. Until both routes qualify on the same configuration, timing them would answer a different question.

## Replay

Use the [pinned packages and solver setup](../schwarz-series.md). Request one numerical thread for BLAS, OpenMP, Accelerate, Rayon and VTK before starting Python. Run each case in a fresh process and use a new output path. From the repository root, for example:

```sh
python demos/solver/schwarz_diagnostics.py default default.json
python demos/solver/schwarz_diagnostics.py finest finest.json
python demos/solver/schwarz_diagnostics.py repaired-solve repaired-30.json --repair-resolution 30
```

The other stage and refinement cases are `tuned`, `single-pass-repair`, `medium`, `fine` and `curved-tet10`. For repaired stiffness, the driver also supports VTK resolutions 34 and 36. Each configuration directory contains the case records, measured hardware, solver inputs and driver hashes. Stage runs include intermediate audits. Use the separate stiffness screens for workflow costs.

| Windows | Linux | macOS |
| --- | --- | --- |
| [Core Ultra 9 285K](core-ultra-9-285k/) | [Ryzen 7 7730U](ryzen-7-7730u/) | [Apple M4](apple-m4/) |
| [Core i7-12700H](core-i7-12700h/) | [Core i7-10750H](core-i7-10750h/) | |
