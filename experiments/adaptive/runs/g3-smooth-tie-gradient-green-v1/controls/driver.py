"""Replay admitted geometry controls after changing pending-wall settling."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                  RAYON_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    controls = [
        ("g3-schwarz-quality-floor-control-v1", ["--no-adapt"], True),
        ("g3-gyroid-curved9-v1", ["--no-adapt"], True),
        ("g3-lattice-adapt24-verified-v1", ["--symmetric-union", "--clip-void", ".02"], False),
        ("g3-spinodoid-curved32-v1", ["--compiled-spinodoid"], False),
    ]
    import meshers._meshers
    import numpy as np

    sources = [Path(__file__), Path(workflow.__file__), workflow.HERE / "admission_pilot.py",
               workflow.ROOT / "crates/meshers-core/src/tet_adapt.rs",
               Path(meshers._meshers.__file__)]
    inputs = {str(path): workflow.digest(path) for path in sources}
    for name, _, _ in controls:
        source = workflow.HERE / "runs" / name
        row = json.loads((source / "result.json").read_text())
        for path in (source / "result.json", source / "mesh.npz", source / "driver.py",
                     workflow.HERE / "runs/g3-microgen-defaults-v4" / (row["case"] + "-samples.npz"),
                     workflow.HERE / "runs/g3-microgen-defaults-v4" / (row["case"] + "-specimen.json")):
            inputs[str(path)] = workflow.digest(path)
    workflow.save(args.run / "manifest.json", {"inputs": inputs, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 600, "scope": "Geometry regression; no new FE or speed claim"})
    results = []
    with Job(parts=[name for name, _, _ in controls]) as job:
        for name, flags, require_identical in controls:
            source = workflow.HERE / "runs" / name
            previous = json.loads((source / "result.json").read_text())
            assert previous["qualified"] and previous["geometry"]["geometry_admitted"]
            output = args.run / name
            output.mkdir()
            for suffix in ("-samples.npz", "-specimen.json"):
                filename = previous["case"] + suffix
                shutil.copyfile(workflow.HERE / "runs/g3-microgen-defaults-v4" / filename, output / filename)
            command = [sys.executable, str(workflow.HERE / "admission_pilot.py"),
                "--run", str(output), "--case", previous["case"], "--worker",
                "--cells", str(previous["cells"]), "--order", str(previous["order"]),
                "--min-length", str(previous["min_length"]),
                "--native-tolerance", str(previous["native_tolerance"]),
                "--optimization-passes", str(previous.get("optimization_passes", 1)),
                "--adapt-tolerance", str(previous.get("adapt_tolerance", .002)), *flags]
            with job.part(name) as part:
                record = workflow.bounded(command, output / "worker.log", part, timeout=600)
                workflow.save(output / "worker.json", record)
                assert record["returncode"] == 0, f"Read {output / 'worker.log'}"
                row = json.loads((output / "result.json").read_text())
                assert row["qualified"] and row["periodic"]["periodic_triangles_match"]
                with np.load(source / "mesh.npz") as old, np.load(output / "mesh.npz") as new:
                    identical = all(np.array_equal(old[key], new[key]) for key in ("points", "tetrahedra", "surface"))
                assert identical or not require_identical, "Nonadaptive control changed"
                results.append({"source": str(source), "result": str(output / "result.json"),
                    "qualified": row["qualified"], "bitwise_identical": identical,
                    "previous_nodes": previous["nodes"], "nodes": row["nodes"],
                    "previous_distance": previous["geometry"]["sampled_two_sided_max"],
                    "distance": row["geometry"]["sampled_two_sided_max"]})
                workflow.save(args.run / "result.json", {"controls": results})
                job.result(name, qualified=int(row["qualified"]), identical=int(identical), nodes=row["nodes"],
                           seconds=record["wall_s"])
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


if __name__ == "__main__":
    main()
