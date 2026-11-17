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

The study uses 0.002 sampled distance, 1% relative volume error, 0.1 minimum quality and 250,000 nodes and tetrahedra each. Schwarz P has exact volume 0.5. The Gyroid volume estimate is 0.3564434, with a range of 0.0001907 across four scrambled Sobol estimates of 262,144 points each. [Frozen volume estimates](https://github.com/kmarchais/meshers/blob/9d156a2982277aed3a28fbcaeb676cc8b7236ca9/experiments/adaptive/runs/g3-microgen-defaults-v4/gyroid-specimen.json)
