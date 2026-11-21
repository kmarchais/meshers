# Periodic bench baseline

`baseline.json` holds the records of the regression suite of `periodic_bench.py`: 33 configurations of periodic implicit solids, meshed on the lattice path and with adaptation. Every case runs at 16 cells with a tolerance of 0.02. The cases most sensitive to where the lattice falls also run at 12 and 24 cells, and two smooth ones at 24 cells with a tolerance of 0.01.

Compare the current build with it from the repository root:

```sh
uv run --with scipy python demos/solver/periodic_bench.py --suite --mode both --baseline demos/solver/periodic-bench/baseline.json
```

The run fails when a configuration does worse than its baseline record: a mesh now rejected, a limit newly missed, a minimum quality lower by more than 0.03, a wall distance larger by more than 15% of the tolerance, a volume error larger by more than 0.2%, 15% more tetrahedra, or a new recourse to defeaturing or to the lattice fallback. A configuration the baseline rejected cannot regress.

The suite takes about twenty minutes on one core. Lattice meshes repeat exactly on one machine. Adaptive meshes of symmetric solids can differ between machines, because rounding breaks their ties, so regenerate the baseline on the machine that compares against it: run the same command with `--output` and a new file name instead of `--baseline`.
