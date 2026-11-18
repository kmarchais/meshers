# Geometry admission for the periodic Tet4 study

Use the Python 3.12 release-wheel environment in the [replay instructions](README.md). Install the observer with `python -m pip install pyvista==0.49.0 vtk==9.7.1 h5py==3.16.0`. From the repository root, run each command with a fresh output directory.

```sh
python demos/solver/audit_release_wheel.py gyroid_uniform audit-gyroid-uniform
python demos/solver/audit_release_wheel.py gyroid_adaptive audit-gyroid-adaptive
python demos/solver/audit_release_wheel.py schwarz_native audit-schwarz-native
```

Each directory contains the exported VTKHDF volume and its audit. A failed gate returns exit code 1 and retains the measurements. To check a Tet4 volume from another producer on the same centered unit cell, run `python demos/solver/periodic_geometry.py schwarz mesh.vtu audit.json`. Use `gyroid` for the band with levels ±0.55.

The observer recovers exterior faces from tetrahedral incidence. It checks face orientation, opposite-cap triangle matches, signed volume, minimum signed MMG quality and the node and element budgets. It reports face-connected components in the cut unit cell and their signed volumes. A seeded 2,048-point Sobol set checks solid and void occupancy away from the 0.002 tolerance layer. This tests component coverage at the sampled locations. It can miss a feature smaller than the sample spacing.

For reverse distance the observer intersects 2,048 lines per axis with each analytic wall and measures distance to the closest triangle on that wall. It samples each cap separately using 2,048 planar points before material filtering. For forward distance it checks vertices, edge midpoints and triangle centers against the matching analytic wall or cap. Its shortest axis-parallel displacement to a valid root is an upper bound and can reject a mesh near the cutoff. Each direction records its worst location. The measurements do not bound the unsampled surface.

The study uses 0.002 sampled distance, 1% relative volume error, 0.1 minimum quality and 250,000 nodes and tetrahedra each. Schwarz P has exact volume 0.5. The Gyroid volume estimate is 0.3564434, with a range of 0.0001907 across four scrambled Sobol estimates of 262,144 points each. [Frozen volume estimates](https://api.github.com/repos/kmarchais/meshers/git/blobs/e846e43479c155bc776dd2ff65ef61ceb5e171bb)

## Windows release-wheel observations

All three exported volumes met the frozen sampled gates in the 17 November run. All periodic triangle traces matched and all interior and exterior witnesses agreed. Components here count shared tetrahedral faces inside the cut box, before identifying periodic images.

| Mesh | Largest forward upper bound | Largest reverse distance | Components |
| --- | ---: | ---: | ---: |
| [Uniform Gyroid](geometry-20261117/gyroid_uniform.json) | 0.001930 | 0.001247 | 3 |
| [Adaptive Gyroid](geometry-20261117/gyroid_adaptive.json) | 0.001807 | 0.001225 | 3 |
| [Adaptive Schwarz P](geometry-20261117/schwarz_native.json) | 0.001216 | 0.000794 | 1 |

## Configuration replays

All 12 of 12 regenerated exports passed the sampled independent geometry checks. The records include the signed volume of each face-connected component, the worst locations on each wall and cap, and the periodic face comparison. The audit uses the same numerical limits on every configuration. [Configuration and runtime details](configurations-20261118/README.md)

| CPU and raw audits | Admitted | Uniform Gyroid tets | Adaptive Gyroid tets | Schwarz P tets |
| --- | ---: | ---: | ---: | ---: |
| [Apple M4](configurations-20261118/apple-m4-091916z-a3fb9327/geometry.json) | 3/3 | 72,757 | 34,854 | 35,994 |
| [Core i7-10750H](configurations-20261118/core-i7-10750h-091913z-4f4ba086/geometry.json) | 3/3 | 72,757 | 34,863 | 35,989 |
| [Core i7-12700H](configurations-20261118/core-i7-12700h-091916z-2ebf2eee/geometry.json) | 3/3 | 72,757 | 34,863 | 35,989 |
| [Ryzen 7 7730U](configurations-20261118/ryzen-7-7730u-091912z-fb87984d/geometry.json) | 3/3 | 72,757 | 34,863 | 35,989 |
