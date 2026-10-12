"""Repeat generation only, requiring exact identity with previously audited meshes."""

import argparse
import importlib.metadata
import json
import os
import statistics
import sys
import threading
import time
from pathlib import Path

for variable in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[variable] = "1"

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import report
from run import digest, native_seed_control, write_json


def repeat(output, references):
    started = time.monotonic()
    live = {"completed": 0, "message": "Verify frozen generation sources", "elapsed": 0}
    stopped = threading.Event()

    def heartbeat():
        while not stopped.is_set():
            report(live["completed"] / 8, live["message"], eta_s=max(0, 170 - (time.monotonic() - started)),
                   seconds=round(time.monotonic() - started), samples=live["completed"])
            stopped.wait(2)

    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    try:
        manifest = json.loads((output / "manifest.json").read_text())
        for relative, expected in manifest["source_sha256"].items():
            if digest(output / "source" / relative) != expected:
                raise ValueError(f"Frozen source changed: {relative}")
        for path, expected in manifest["dependency_sources_sha256"].items():
            if digest(path) != expected:
                raise ValueError(f"Dependency changed: {path}")
        for name, version in manifest["packages"].items():
            if importlib.metadata.version(name) != version:
                raise ValueError(f"Dependency version changed: {name}")

        import meshers._meshers
        import numpy as np
        from microgen_benchmark.backends import generate
        from microgen_benchmark.cases import Case
        from microgen_benchmark.geometry import canonical_mesh

        if digest(meshers._meshers.__file__) != manifest["native_sha256"]:
            raise ValueError("Imported native binary differs from the frozen candidate")
        if len(references) != 2:
            raise ValueError("Expected two admitted periodic-Gyroid artifacts")
        controls = []
        for path in references:
            result = json.loads((path / "result.json").read_text())
            if result["case"] != "periodic-gyroid":
                raise ValueError(f"Expected the shared periodic Gyroid case: {path}")
            if not all(result.get(k) for k in ["homogenization_audited", "geometry_eligible", "quality_eligible", "periodic_eligible"]):
                raise ValueError(f"Reference is not admitted: {path}")
            if digest(path / "mesh.npz") != result["mesh_sha256"]:
                raise ValueError(f"Audited mesh changed: {path}")
            controls.append((path, result, np.load(path / "mesh.npz")))
        case = Case("periodic-gyroid", "periodic-gyroid", "gyroid", 1.1,
                    "Periodic unit cell", "Shared analytic field")
        document = {"state": "running", "manifest_sha256": digest(output / "manifest.json"), "rows": [],
                    "scope": "One warmup and three generation samples; exact canonical mesh identity to previously audited FEA artifacts; no repeated solver timing"}
        for repetition in range(4):
            for index in ([0, 1] if repetition % 2 == 0 else [1, 0]):
                path, reference, saved = controls[index]
                workflow, cells = reference["workflow"], reference["cells"]
                directory = output / f"repeat-{repetition}-{workflow}"
                directory.mkdir()
                live["message"] = f"{'Warmup' if repetition == 0 else 'Sample ' + str(repetition)}: {workflow} {cells}"
                if workflow == "mmg-native-seed":
                    points, tets, metrics = native_seed_control(case, cells, manifest["config"], directory, lambda **v: None)
                else:
                    points, tets, metrics = generate(case, workflow, cells, manifest["config"], lambda **v: None)
                points, tets, _, _ = canonical_mesh(points, tets)
                if not np.array_equal(points, saved["points"]) or not np.array_equal(tets, saved["tetrahedra"]):
                    raise ValueError(f"Generated mesh differs from audited artifact: {workflow}, repeat {repetition}")
                document["rows"].append({"workflow": workflow, "cells": cells, "warmup": repetition == 0,
                                         "mesh_s": metrics["mesh_s"], "tetrahedra": len(tets),
                                         "audited_artifact": str(path), "audited_mesh_sha256": reference["mesh_sha256"], "exact_mesh_identity": True})
                live["completed"] += 1
                write_json(output / "generation-repeats.json", document)
        document["medians_s"] = {workflow: statistics.median(row["mesh_s"] for row in document["rows"]
                                     if row["workflow"] == workflow and not row["warmup"])
                                     for workflow in {row["workflow"] for row in document["rows"]}}
        document.update(state="complete", elapsed_s=time.monotonic() - started)
        write_json(output / "generation-repeats.json", document)
        print(json.dumps(document["medians_s"]), flush=True)
        report(1, "Repeated generation matches both audited meshes", eta_s=0, samples=8)
    finally:
        stopped.set()
        thread.join()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--references", type=Path, nargs=2, required=True)
    args = parser.parse_args()
    repeat(args.run.resolve(), [path.resolve() for path in args.references])
