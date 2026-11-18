# Release-wheel Tet4 comparison

This experiment starts with the shipped Meshers 0.2.1 CPU wheel. It replays three native cases on the unit periodic cell centered at the origin, then translates coordinates by half a unit for the existing fedoo verifier. Each case records generation time, signed volume, cell and node counts, six-load stiffness, reaction and energy checks, solver residuals and free DOFs. Solver time includes the verifier's second, stress-controlled solve. Keep the same solver and verification costs on both arms when making timing comparisons.

## Frozen evidence

[baseline.json](baseline.json) identifies its research commit and the Git blob ID of the source record for each route. Read a record with `git cat-file blob <source_blob>` in a checkout containing the research commit. The Gyroid is the raw trigonometric band ±0.55. The historical uniform/adaptive Tet4 meshes used 72,757/34,854 cells and differed by 2.92917% in worst mixed-load energy. The archived tensors are discrete reference results. Their difference measures agreement, not continuum error.

Schwarz P uses the nonpositive cosine sum. Its analytic solid fraction is 0.5. The admitted Microgen control used a resolution-30 VTK seed, one MMG pass with `nofem`, `hausd=0.0005`, `hmax=0.1` and `hgrad=1.3`, followed by interior repair. It yielded 136,511 tetrahedra and 29,192 nodes. The unmodified Schwarz control failed its quality floor. The actual resolution-40 Gyroid control retained both public MMG passes and failed the quality, volume and facet checks. Preserve these as measured failures.

The frozen cutoffs are 0.002 sampled two-sided boundary distance, 1% absolute relative volume error, minimum quality 0.1, positive element orientation, periodic nodes and triangulations, and 250,000 nodes/tets each. Gyroid's volume reference is a four-scramble Sobol estimate. The native replay checks volume, quality at generation, orientation, periodicity and mechanics. A fresh independent two-sided distance audit is still required for comparative admission. Native timings cover direct Meshers generation. Microgen construction, conversion and MMG repair must be charged in a live workflow comparison. Speeds from different machines or different stages cannot establish a win.

## Replay

Install the released wheel and solver in a clean Python 3.12 environment. Run the following commands from the repository root inside the shared compute queue, using a fresh output filename for each case. Set OPENBLAS_NUM_THREADS, OMP_NUM_THREADS and MKL_NUM_THREADS to 1 before starting Python.

```sh
python -m pip install --only-binary=:all: meshers==0.2.1 fedoo==1.0.0b2 numpy==2.5.3 scipy==1.18.1
python demos/solver/release_wheel_comparison.py gyroid_uniform gyroid-uniform.json
python demos/solver/release_wheel_comparison.py gyroid_adaptive gyroid-adaptive.json
python demos/solver/release_wheel_comparison.py schwarz_native schwarz-native.json
```

For each fresh tensor report the worst absolute generalized eigenvalue of the symmetric tensor difference relative to the reference stiffness. To compare adaptive and uniform, use uniform as reference. Preserve worker failures and raw JSON. Before changing defaults, repeat at matched physical error on an unseen specimen and charge all required preprocessing, retries and solves.

The original Microgen comparison driver, dependency manifests and settings are on research commit `9d156a2982277aed3a28fbcaeb676cc8b7236ca9`. In that research environment the following command replays the Schwarz default and repaired arms. Use the archived manifests there to provision Microgen, VTK, mmgpy, its MMG executable and the frozen fixtures. A plain release-wheel environment only supports the native commands above. [Schwarz manifest](https://api.github.com/repos/kmarchais/meshers/git/blobs/dbadf98e26353c9bb7bc5a2b799713422a0fb6f5) · [Gyroid stiffness assessment](https://api.github.com/repos/kmarchais/meshers/git/blobs/0fdf98be1e1fd5f5db25e719c5442dd1e233e5cb)

```sh
python experiments/adaptive/mesh_only_comparison.py --run schwarz-replay --case schwarz-solid --strengthened --warmups 1 --repetitions 5 --select actual-default,actual-explicit-tuned
```

The [cross-platform replay](fleet-20261118/README.md) gives three repeats per case on the configurations listed in that study, with raw timings and stiffness tensors.

## Windows release-wheel run, 17 November 2026

The original Windows run measured three cases sequentially with one numerical thread. Each passed six prescribed-strain loads and six stress-controlled verification solves. The largest recomputed relative linear residual across all 36 systems was 1.07e-12. The export test also passed on retry with a writable temporary directory.

| Tet4 mesh | Nodes | Tetrahedra | Free displacement DOFs (xx load) | Generation (s) | Homogenization and checks (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| [Uniform Gyroid](windows-0.2.1-20261117/gyroid_uniform.json) | 17,647 | 72,757 | 49,941 | 0.83 | 27.97 |
| [Adaptive Gyroid](windows-0.2.1-20261117/gyroid_adaptive.json) | 9,435 | 34,863 | 25,977 | 1.74 | 7.49 |
| [Adaptive Schwarz P](windows-0.2.1-20261117/schwarz_native.json) | 8,700 | 35,989 | 23,145 | 0.97 | 5.13 |

Adaptive Gyroid used 52.1% fewer tetrahedra. The worst relative strain-energy change was +2.93861% across mixed loads against the fresh uniform tensor, compared with 2.92917% in the archived study. This is the largest absolute generalized eigenvalue of the symmetric tensor difference with uniform stiffness as the metric. It measures disagreement between two discrete solutions. Establishing which is closer to the physical solution requires an independently converged reference.

Times are single observations in one process, in table order. The solver column includes the verification solves and audits inside homogenization. Independent two-sided distance admission and repeated equal-error workflow timings remain open. No live Microgen/MMG arm ran in this job.

Run the [independent geometry audit](geometry-admission.md) on regenerated exports before using solver results in a comparison.

The [live Schwarz P refinement series](schwarz-series.md) measures all three workflows against the same stiffness reference.
