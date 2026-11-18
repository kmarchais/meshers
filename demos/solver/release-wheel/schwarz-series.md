# Live Schwarz P stiffness comparison

This series asks how much time each Tet4 workflow needs to come within 1% of an archived Schwarz P effective stiffness tensor. It uses the negative cosine sum on the centered unit cube, three periodic directions, isotropic linear elasticity with E = 2.5 and Î½ = 0.25, and six unit engineering strains. The fedoo verifier repeats the six loads under stress control. Every route passes through the same [geometry audit](geometry-admission.md) before solving.

Use Python 3.12 and the [release-wheel solver environment](README.md). Install the remaining study dependencies and run the series from the repository root on an idle host. Reserve space for 27 exported volumes. Each worker runs in a new process with one numerical thread. Import time is excluded from the reported workflow time.

```sh
python -m pip install microgen==2.0.0b1 mmgpy==0.17.1 vtk==9.7.1 pyvista==0.49.0 h5py==3.16.0 psutil==7.2.2
python demos/solver/schwarz_series.py schwarz-study
```

| Route | Refinement levels | What runs |
| --- | --- | --- |
| Native | 12, 16, 24 background cells per axis | Frozen adaptive Meshers 0.2.1 recipe. |
| Microgen default MMG | Public default VTK resolution (20), then 24 and 30 | Original VTK lower skeletal solid followed by Microgen's two public MMG passes. |
| Microgen tuned MMG | VTK resolution 18, 24, 30 | Same two passes. The second receives `hausd=0.001`, `hmax=0.06`, `hgrad=1.3`. The historical single-pass recipe with interior repair remains a separate archived control. |

The worker checks the imported VTK mesher source against the original study hashes, uses the MMG executable packaged with mmgpy, and records its hash and each call. On Windows it allows Microgen's temporary mesh file to be reopened by the reader. Repetitions rotate the execution order. The default is three observations at each level and a ten-minute worker timeout. The summary lists all failures and reports median cost only for levels reaching the target in every repetition. Use the individual results to identify the measured levels on either side of the target. If no admitted level reaches it, report that outcome.

Each record includes meshing and I/O time, geometry audit time, assembly time, CG time per load, total workflow time, free displacement DOFs and sampled peak process-tree RSS. RSS is sampled every quarter second and can include shared pages in more than one process. The local quantity is Tet4 Ïƒxx at the fixed interior point (âˆ’0.364, âˆ’0.347, âˆ’0.313) under each of the six strains. Its variation across levels can be measured independently of effective stiffness.

The stiffness discrepancy is the largest absolute generalized eigenvalue of the symmetric tensor difference in the reference energy metric. The frozen [reference](schwarz-reference.json) uses fourth-order extrapolation of two quadratic discretizations. The fine-to-extrapolated difference is 0.0130% and changing the extrapolation window moves the estimate by 0.0130%. These values measure reference sensitivity. They leave continuum error and pointwise stress accuracy open. [Original reference record](https://api.github.com/repos/kmarchais/meshers/git/blobs/a2b4befa3143b904a104cef143e730eb11d6ad41)

## Windows results, 18 November 2026

All 27 workers returned measurements. None of the nine route and resolution combinations qualified for a workflow timing at the 1% stiffness target. Each native run passed geometry and completed the solver checks. Every Microgen run exceeded the 0.002 reverse wall-distance limit before the solver stage. Five of the six Microgen settings also missed the 0.1 minimum element quality. Both default resolution 20 and tuned resolution 18 exceeded the cap tolerance.

The table gives the largest sampled reference-to-mesh wall distance and the lowest element quality across three repeats. Stiffness discrepancy is relative to the extrapolated Schwarz P tensor.

| Route | Input resolution | Reverse wall distance | Minimum quality | Stiffness discrepancy |
| --- | ---: | ---: | ---: | --- |
| Microgen default | 24 | 0.008245 | 0.0224 | Skipped after geometry rejection |
| Microgen default | 30 | 0.006826 | 0.0585 | Skipped after geometry rejection |
| Microgen default | 20 (default) | 0.007045 | 0.0548 | Skipped after geometry rejection |
| Microgen tuned | 18 | 0.007486 | 0.1446 | Skipped after geometry rejection |
| Microgen tuned | 24 | 0.005813 | 0.0215 | Skipped after geometry rejection |
| Microgen tuned | 30 | 0.006736 | 0.0812 | Skipped after geometry rejection |
| Meshers | 12 | 0.000940 | 0.3802 | 2.664% |
| Meshers | 16 | 0.000794 | 0.4261 | 2.449% |
| Meshers | 24 | 0.000721 | 0.4607 | 2.304% |

Doubling native background resolution from 12 to 24 increased the final Tet4 count from 35,192 to 36,433. With the adaptive sizing parameters fixed, that sweep gave little evidence about how stiffness converges with increasing degrees of freedom. The Microgen reverse-distance failures independently confirm a wall mismatch even with the conservative forward bound removed from consideration. The next comparison should measure the VTK seed and each MMG pass separately, then vary the native size controls until the tested meshes span the stiffness target. [Series summary](schwarz-20261118/summary.json) · [All 27 worker records](schwarz-20261118/)

## Configuration replays

Across 5 configurations, 45 native workers completed the geometry checks and solves but missed the stiffness target. Geometry rejected 90 Microgen volumes before FEA. No tested level qualified for a timing at 1%. Each raw archive contains 27 worker records, including failed geometry checks. [Hardware, OS and run dates](configurations-20261118/README.md)

| CPU and raw workers | Microgen rejected by geometry | Native above target | Qualifying levels | Smallest native discrepancy |
| --- | ---: | ---: | ---: | ---: |
| [Apple M4](configurations-20261118/apple-m4-091916z-a3fb9327/schwarz.json) | 18 | 9 | 0/9 | 2.283% |
| [Core Ultra 9 285K](configurations-20261118/core-ultra-9-285k-095358z-610e8c8f/schwarz.json) | 18 | 9 | 0/9 | 2.304% |
| [Core i7-10750H](configurations-20261118/core-i7-10750h-091913z-4f4ba086/schwarz.json) | 18 | 9 | 0/9 | 2.336% |
| [Core i7-12700H](configurations-20261118/core-i7-12700h-091916z-2ebf2eee/schwarz.json) | 18 | 9 | 0/9 | 2.304% |
| [Ryzen 7 7730U](configurations-20261118/ryzen-7-7730u-091912z-fb87984d/schwarz.json) | 18 | 9 | 0/9 | 2.336% |

On the two Linux configurations, 54 initial workers exited during the Gmsh import because the OpenGL utility library was missing. Their records are in `first_attempt`. The measured retries used the distribution's libGLU and libOpenGL in the benchmark environment. Linux wheel installations also need execute permission on the MMG command-line programs.

The [stage audits and native sizing study](schwarz-diagnostics-20261118/README.md) locate the first wall failure and bracket the native stiffness target.
