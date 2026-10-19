"""Refine the fixed lattice's independent volume estimate without rewriting gates."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import sys
import time

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def deadline(_signal, _frame):
    raise TimeoutError("Volume replicate exceeded its 300-second stage cap")


def half_cube_count(points):
    import numpy as np

    return int(np.count_nonzero(points[:, 0] < 0))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (workflow.MEMORY, workflow.MEMORY))
    signal.signal(signal.SIGALRM, deadline)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    np, _, _, microgen, _, scorecard = workflow.imports()
    import psutil
    from scipy.stats import qmc

    control = qmc.Sobol(3, scramble=False).random_base2(1) - .5
    assert half_cube_count(control) == 1
    source = workflow.HERE / "runs/g3-microgen-defaults-v4/smooth-lattice-specimen.json"
    old = json.loads(source.read_text())
    case, shape, _ = workflow.specimen("smooth-lattice")
    probes = np.random.default_rng(20261019).uniform(-.5, .5, (4096, 3))
    np.testing.assert_allclose(case.function(*probes.T), shape.evaluate(*probes.T), atol=1e-12, rtol=1e-12)
    candidates = [workflow.HERE / "runs" / name / "result.json" for name in (
        "g3-lattice-curved24-symmetric-v1", "g3-lattice-adapt24-tol0028-v1",
        "g3-lattice-live-energy-floor025-v1/final", "g3-lattice-microgen-straight-tet10-hausd001-amg-v2")]
    paths = [Path(__file__), Path(workflow.__file__), Path(scorecard.__file__), source, *candidates,
             Path(microgen.shape.implicit_ops.__file__)]
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "case": "smooth-lattice",
        "memory_bytes": workflow.MEMORY, "threads": 1, "replicate_timeout_s": 300,
        "seeds": [0, 1, 2, 3], "batch_points": 2**18, "nested_powers": list(range(18, 25)),
        "maximum_points_per_replicate": 2**24, "total_points": 4 * 2**24,
        "maximum_final_spread": 3e-5, "maximum_final_mean_prefix_shift": 1e-5,
        "definition": old["definition"], "corners_in_order": np.asarray(scorecard.CORNERS).tolist(),
        "versions": {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "microgen")},
        "scope": "Independent QMC volume audit only; frozen outputs and all scientific admission gates unchanged"})
    process = psutil.Process()
    peak = process.memory_info().rss
    rows = []
    with Job(parts=[f"seed-{seed}" for seed in range(4)]) as job:
        for seed in range(4):
            with job.part(f"seed-{seed}") as part:
                started = time.perf_counter()
                signal.alarm(300)
                try:
                    engine = qmc.Sobol(3, scramble=True, seed=seed)
                    material, half_cube, completed, checkpoints = 0, 0, 0, []
                    for batch in range(64):
                        points = engine.random(2**18) - .5
                        values = case.function(*points.T)
                        assert np.isfinite(values).all()
                        material += int(np.count_nonzero((values >= case.levels[0]) & (values <= case.levels[1])))
                        half_cube += half_cube_count(points)
                        completed += len(points)
                        if completed & (completed - 1) == 0:
                            assert half_cube == completed // 2
                            estimate = material / completed
                            if completed == 2**18:
                                assert estimate == old["volume_sobol_estimates"][seed], "Frozen prefix changed"
                            checkpoints.append({"power": completed.bit_length() - 1,
                                "points": completed, "material_points": material, "volume": estimate})
                            workflow.save(args.run / f"seed-{seed}.json", {"seed": seed, "checkpoints": checkpoints})
                        rss = process.memory_info().rss
                        peak = max(peak, rss)
                        assert rss <= workflow.MEMORY
                        elapsed = time.perf_counter() - started
                        part.update((batch + 1) / 64, f"Seed {seed}: {completed:,} / {2**24:,} points",
                                    volume=material / completed, points_per_s=completed / max(elapsed, 1e-9),
                                    rss_gib=rss / 1024**3, peak_gib=peak / 1024**3)
                    row = {"seed": seed, "volume": material / completed,
                        "checkpoints": checkpoints, "seconds": time.perf_counter() - started,
                        "frozen_prefix_bitwise_identical": True, "half_cube_control_exact": True}
                    workflow.save(args.run / f"seed-{seed}.json", row)
                    rows.append(row)
                    job.result(f"seed-{seed}", volume=row["volume"], seconds=row["seconds"],
                               peak_gib=peak / 1024**3)
                finally:
                    signal.alarm(0)
    levels = []
    for power in range(18, 25):
        estimates = np.array([next(checkpoint["volume"] for checkpoint in row["checkpoints"]
                                  if checkpoint["power"] == power) for row in rows])
        levels.append({"power": power, "points_per_replicate": 2**power, "estimates": estimates.tolist(),
            "mean": float(estimates.mean()), "spread": float(np.ptp(estimates)),
            "replicate_standard_deviation": float(estimates.std(ddof=1)),
            "estimated_standard_error_of_mean": float(estimates.std(ddof=1) / 2)})
    final = levels[-1]
    final_shift = abs(final["mean"] - levels[-2]["mean"])
    minimum, maximum = min(final["estimates"]), max(final["estimates"])
    comparisons = []
    for path in candidates:
        row = json.loads(path.read_text())
        volume = row["geometry"]["volume"]
        comparisons.append({"source": str(path), "mesh_volume": volume,
            "frozen_relative_volume_error": row["geometry"]["relative_volume_error"],
            "relative_volume_error_at_new_mean": volume / final["mean"] - 1,
            "relative_error_across_replicate_range": [volume / maximum - 1, volume / minimum - 1],
            "passes_one_percent_volume_gate_across_replicate_range":
                max(abs(volume / maximum - 1), abs(volume / minimum - 1)) <= .01})
    workflow.save(args.run / "result.json", {"levels": levels, "frozen_reference_volume": old["volume"],
        "new_mean_minus_frozen_reference": final["mean"] - old["volume"],
        "final_mean_prefix_shift": final_shift,
        "predeclared_volume_precision_screen_passed": final["spread"] <= 3e-5 and final_shift <= 1e-5,
        "comparisons": comparisons, "peak_rss_bytes": peak,
        "scope": "Empirical nested scrambled-QMC stability, not a rigorous volume interval, stiffness reference, local-stress certificate or permission to modify frozen results"})
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
