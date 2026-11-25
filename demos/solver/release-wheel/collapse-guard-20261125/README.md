# Tet collapse guard order

The directed collapse from vertex u to v already rejects a replacement edge longer than its target. Checking the neighbors of u first avoids building the shell and testing the links for those rejected collapses. Accepted collapses still run every test before editing the mesh.

On the earlier base `aeb3ef09`, a sampled profile of the adaptive six-mode Fourier fixture attributed 44.3–45.0% of directed-collapse time to the volume-link block. More than half of the attempts reaching the old length check failed there. The following release runs use the merged base `5ed0ae95` and guard change `2ac044a1`.

## Measurements

Wall times and peak memory are medians of 13 measured runs per build. The percentage and its 95% interval use the ratios within pairs.

| CPU | Fixture | Base (s) | Candidate (s) | Paired wall change [95% interval] | Base / candidate peak memory (MiB) |
| --- | --- | ---: | ---: | ---: | ---: |
| Ryzen 7 7730U | Adaptive Fourier | 2.527 | 2.387 | -5.6% [-6.5%, -4.4%] | 104.71 / 104.59 |
| Ryzen 7 7730U | Uniform Fourier | 1.780 | 1.783 | +0.5% [-0.7%, +1.5%] | 90.22 / 90.15 |
| Ryzen 7 7730U | Adaptive lattice | 1.346 | 1.260 | -6.5% [-7.6%, -5.4%] | 76.79 / 78.71 |
| Core Ultra 9 285K | Adaptive Fourier | 1.764 | 1.683 | -6.2% [-7.4%, -4.4%] | 81.31 / 81.35 |
| Core Ultra 9 285K | Uniform Fourier | 0.910 | 0.913 | +1.0% [-1.0%, +2.3%] | 63.06 / 63.35 |
| Core Ultra 9 285K | Adaptive lattice | 0.925 | 0.850 | -7.0% [-10.0%, -5.6%] | 58.83 / 58.79 |

Both adaptive Fourier intervals clear the 1.5% criterion. The uniform Fourier point estimates are small slowdowns, and both intervals include zero. These timings cover single-threaded generation on these fixtures.

Each host passed all 132 ordered-mesh and release checks. The [Linux record](linux.json) and [Windows record](windows.json) give mesh digests, numerical checks, memory and build details. The [Linux CSV](linux.csv) and [Windows CSV](windows.csv) give individual samples.

## Replay

Use the same interpreter, NumPy, Rust toolchain and maturin version for both builds on a host. Each host's JSON records those versions, the native extension SHA-256 values and the hardware. The mesher and validator use float64. Build without profiling. The commands below reproduce the Linux cycle from the checkout containing the comparison driver, with its shared compute slot reserved. The Windows run used the default wheel compatibility and logical CPU 10.

```sh
set -eu
trial=$(mktemp -d)
uv venv --python 3.12.16 "$trial/env"
py="$trial/env/bin/python"
uv pip install --python "$py" numpy==2.5.3
for side in base candidate; do
    case "$side" in
        base) ref=5ed0ae95223eddede731bfcee8901d7d7d75ae72 ;;
        candidate) ref=2ac044a1988ede1f0f455969c8b2cf235220bf6b ;;
    esac
    mkdir "$trial/$side"
    git archive "$ref" Cargo.toml Cargo.lock LICENSE README.md crates | tar -x -C "$trial/$side"
    uvx --from maturin==1.16.0 maturin build --locked --release --compatibility linux \
        --manifest-path "$trial/$side/crates/meshers-python/Cargo.toml" \
        --interpreter "$py" --out "$trial/$side-wheel"
    uv run --no-project --python "$py" python -m zipfile -e \
        "$trial/$side-wheel/"*.whl "$trial/$side-packages"
done
taskset -c 2 uv run --no-project --python "$py" python -m tools.collapse_compare "$trial"
```

Each fixture starts with a discarded warmup pair, then seven base/base pairs. The base/candidate cycle starts with another warmup pair and measures 13 pairs, reversing order on odd blocks. Each sample starts a fresh one-thread process. Its wall timer covers generation, including field preparation. Peak resident memory includes imports, and the release validator runs on the returned arrays after timing.

The CSV files preserve execution order. Block -1 is warmup, side A is always the base, and side B uses the build named in the comparison column. For a case's wall change, take the median of the 13 within-pair B/A ratios and subtract one. For the 95% interval, resample pairs within each case 4,000 times with seed zero. Widen the percentile bounds around the point estimate by at least 1% on either side, or by 1.96 times the sample standard deviation of the seven A/A log ratios divided by the square root of 13, whichever is larger.

The predeclared timing criterion requires the entire adaptive Fourier interval to be at least 1.5% faster on each host, with at most a 3% median slowdown on either control. The release memory limit is 15% plus 8 MiB. Each run must also reproduce its host's base digests for point coordinates, tetrahedra, surface triangles, boundary tags and all three ordered periodic-pair arrays, and exactly match the release validator's mesh metrics.

## Schwarz P solver check

Use the candidate checkout's `native-candidate` release gate for both wheels. It writes a 24-cell adaptive Tet4 mesh to VTKHDF, reads it back, and solves six strain-controlled and six stress-controlled periodic load cases. The energy criterion is a maximum discrepancy of 1% against the reference stiffness tensor.

```sh
uv pip install --python "$py" fedoo==1.0.0b2 vtk==9.7.1 pyvista==0.49.0 scipy==1.18.1
for side in base candidate; do
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=1 \
    PYTHONPATH="$trial/$side-packages" taskset -c 2 \
        uv run --no-project --python "$py" python demos/solver/schwarz_workflow.py \
        native-candidate "$trial/$side-solver"
done
```

On each host, both wheels produced 42,320 points and 189,873 tetrahedra. Their exported HDF5 datasets and attributes matched exactly, as did the reopened solver mesh, geometry audits, stiffness tensor and every numerical solver diagnostic. The HDF5 files themselves had different byte hashes. Both wheels reached 0.708% energy discrepancy on each host. [Linux solver results](linux-solver.json) and [Windows solver results](windows-solver.json).

The solver records include hashes of the geometry audit and FEA output. To reproduce them, drop `seconds` and keys ending in `_seconds` recursively from the FEA object. Serialize each object with Python's `json.dumps(obj, sort_keys=True, separators=(",", ":"))` and hash its UTF-8 bytes with SHA-256. The VTKHDF records hash each dataset and attribute with the same dtype, shape and byte-order scheme as the mesh sample.
