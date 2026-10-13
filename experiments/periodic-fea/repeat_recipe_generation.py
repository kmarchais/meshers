"""Compare live generation recipes while retaining each admitted build and config.

One warmup and three samples per recipe run in alternating order. Each sample
must reproduce the exact coordinates and connectivity of its admitted artifact.
Imports, process startup and the independent identity check are outside timing.
"""

import argparse
import importlib.metadata
import json
import os
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

from run import digest, write_json


def prepare(output, references):
    from reference_sensitivity import read

    output.mkdir(parents=True, exist_ok=False)
    source = output / "source"
    source.mkdir()
    for name in [Path(__file__).name, "run.py", "reference_sensitivity.py", "periodic.py"]:
        shutil.copy2(Path(__file__).parent / name, source / name)
    recipes, inputs = [], {}
    for index, artifact in enumerate(references):
        row, _, _ = read(artifact / "result.json", geometry_gate=.002)
        if row["workflow"] not in {"meshers-native-body", "mmg-native-seed"}:
            raise ValueError("Expected admitted native-body and standard-MMG recipes")
        original = artifact.parent
        manifest = json.loads((original / "manifest.json").read_text())
        target = output / "recipes" / str(index)
        shutil.copytree(original / "source", target, ignore=shutil.ignore_patterns("__pycache__"))
        for name, sha in manifest["source_sha256"].items():
            if digest(target / name) != sha:
                raise ValueError(f"Changed archived recipe source: {name}")
        paths = [artifact / "result.json", artifact / "mesh.npz", original / "manifest.json",
                 original / "results.json", original / "artifact-verification.json"]
        inputs.update({str(path): digest(path) for path in paths})
        recipes.append({"artifact": str(artifact), "source": str(target), "manifest": manifest,
                        "row": row})
    if {r["row"]["workflow"] for r in recipes} != {"meshers-native-body", "mmg-native-seed"}:
        raise ValueError("Need one native-body and one MMG recipe")
    if recipes[0]["row"]["case"] != recipes[1]["row"]["case"]:
        raise ValueError("Recipes must mesh the same specimen")
    case = recipes[0]["row"]["case"]
    if recipes[0]["manifest"]["specimens"][case] != recipes[1]["manifest"]["specimens"][case]:
        raise ValueError("Recipe field definitions differ")
    hashes = {str(p.relative_to(output)): digest(p) for folder in [source, output / "recipes"]
              for p in folder.rglob("*") if p.is_file()}
    write_json(output / "manifest.json", {"state": "prepared", "recipes": recipes,
                                         "source_sha256": hashes, "input_sha256": inputs,
                                         "repetitions": 4, "memory_limit_mib": 16384,
                                         "worker_timeout_s": 1200, "scope": __doc__})
    print(f"Prepared two live recipes: {len(hashes)} source and {len(inputs)} input hashes")


def sample(request_path):
    request = json.loads(request_path.read_text())
    recipe = request["recipe"]
    source = Path(recipe["source"])
    sys.path.insert(0, str(source))
    import importlib.util
    import numpy as np
    import meshers._meshers
    from microgen_benchmark.backends import generate
    from microgen_benchmark.cases import Case
    from microgen_benchmark.geometry import canonical_mesh

    manifest, row = recipe["manifest"], recipe["row"]
    if digest(meshers._meshers.__file__) != manifest["native_sha256"]:
        raise ValueError("Wrong native recipe build")
    for name, version in manifest["packages"].items():
        if importlib.metadata.version(name) != version:
            raise ValueError(f"Changed recipe dependency: {name}")
    for name, sha in manifest["dependency_sources_sha256"].items():
        if digest(name) != sha:
            raise ValueError(f"Changed recipe dependency source: {name}")
    specimen = manifest["specimens"][row["case"]]
    if (source / "fixtures.py").is_file():
        from fixtures import make_case
        case = make_case(row["case"], specimen)
    else:
        # Older archived workers constructed the benchmark Case directly.
        if specimen.get("microgen_callable") or specimen["function"] == "spinodoid_spectral":
            raise ValueError("Older recipe has no helper for this Microgen fixture")
        case = Case(row["case"], row["case"], specimen["function"], specimen["offset"],
                    "Periodic unit cell", "Shared analytic field")
    config = dict(manifest["config"])
    config["geometry_tolerance"] = config.get("native_geometry_tolerance", config["geometry_tolerance"])
    if row["workflow"] == "meshers-native-body":
        config["interior_refine_nodes"] = config["body_refine_nodes"]
        points, tets, metrics = generate(case, "meshers-pass1", row["cells"], config)
    else:
        spec = importlib.util.spec_from_file_location("recipe_run", source / "run.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        points, tets, metrics = module.native_seed_control(case, row["cells"], config,
                                                         request_path.parent, lambda **v: None)
    points, tets, _, _ = canonical_mesh(points, tets)
    with np.load(Path(recipe["artifact"]) / "mesh.npz") as saved:
        if not np.array_equal(points, saved["points"]) or not np.array_equal(tets, saved["tetrahedra"]):
            raise ValueError("Live recipe differs from its admitted mesh")
    write_json(request_path.parent / "sample.json", {"workflow": row["workflow"],
               "native_sha256": manifest["native_sha256"], "generation_s": metrics["generation_s"],
               "mesh_s": metrics["mesh_s"], "tetrahedra": len(tets), "nodes": len(points),
               "exact_mesh_identity": True, "audited_mesh_sha256": row["mesh_sha256"]})


def repeat(output):
    import psutil

    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import report

    report(0, "Verify frozen live generation recipes", eta_s=210)
    manifest = json.loads((output / "manifest.json").read_text())
    for name, sha in manifest["source_sha256"].items():
        if digest(output / name) != sha:
            raise ValueError(f"Changed frozen source: {name}")
    for name, sha in manifest["input_sha256"].items():
        if digest(name) != sha:
            raise ValueError(f"Changed admitted input: {name}")
    document = {"state": "running", "manifest_sha256": digest(output / "manifest.json"),
                "scope": __doc__, "rows": []}
    total = manifest["repetitions"] * len(manifest["recipes"])
    started = time.monotonic()
    for repetition in range(manifest["repetitions"]):
        for index in ([0, 1] if repetition % 2 == 0 else [1, 0]):
            recipe = manifest["recipes"][index]
            directory = output / f"sample-{repetition}-{index}"
            directory.mkdir()
            request = directory / "request.json"
            write_json(request, {"recipe": recipe})
            env = os.environ.copy()
            env["PYTHONPATH"] = recipe["source"]
            for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
                env[name] = "1"
            peak, worker_started = 0, time.monotonic()
            with (directory / "worker.log").open("w", encoding="utf-8") as log:
                process = subprocess.Popen([recipe["manifest"]["interpreter"], "-u", str(Path(__file__)),
                                            "--sample", str(request)], env=env, stdout=log, stderr=subprocess.STDOUT)
                monitored = psutil.Process(process.pid)
                try:
                    while process.poll() is None:
                        try:
                            peak = max(peak, sum(p.memory_info().rss for p in [monitored, *monitored.children(recursive=True)]))
                        except psutil.NoSuchProcess:
                            pass
                        done = len(document["rows"])
                        elapsed = time.monotonic() - worker_started
                        report(done / total, f"{'Warmup' if repetition == 0 else 'Sample ' + str(repetition)}: {recipe['row']['workflow']}",
                               eta_s=(time.monotonic() - started) / max(done, 1) * (total - done),
                               rss_mib=round(peak / 2**20), seconds=round(elapsed))
                        if elapsed > manifest["worker_timeout_s"] or peak > manifest["memory_limit_mib"] * 2**20:
                            raise RuntimeError("Generation worker exceeded unchanged timeout/memory limits")
                        time.sleep(1)
                finally:
                    if process.poll() is None:
                        try:
                            children = monitored.children(recursive=True)
                        except psutil.NoSuchProcess:
                            children = []
                        for child in children:
                            try:
                                child.kill()
                            except psutil.NoSuchProcess:
                                pass
                        process.kill()
                        process.wait()
            if process.returncode:
                raise RuntimeError(f"Recipe failed: inspect {directory / 'worker.log'}")
            row = json.loads((directory / "sample.json").read_text())
            row.update(warmup=repetition == 0, repetition=repetition, peak_rss_mib=peak / 2**20)
            document["rows"].append(row)
            write_json(output / "generation-repeats.json", document)
    for name, sha in manifest["source_sha256"].items():
        if digest(output / name) != sha:
            raise ValueError(f"Frozen source changed during repeats: {name}")
    for name, sha in manifest["input_sha256"].items():
        if digest(name) != sha:
            raise ValueError(f"Admitted input changed during repeats: {name}")
    document["medians_s"] = {workflow: statistics.median(row["generation_s"] for row in document["rows"]
                             if row["workflow"] == workflow and not row["warmup"])
                             for workflow in {row["workflow"] for row in document["rows"]}}
    document.update(state="complete", inputs_unchanged=True, elapsed_s=time.monotonic() - started)
    write_json(output / "generation-repeats.json", document)
    report(1, "Both live recipes reproduce their admitted meshes in all eight runs", eta_s=0, samples=total)
    print(json.dumps(document["medians_s"]), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare", type=Path)
    mode.add_argument("--run", type=Path)
    mode.add_argument("--sample", type=Path)
    parser.add_argument("--references", type=Path, nargs=2)
    args = parser.parse_args()
    if args.prepare:
        if not args.references:
            parser.error("--prepare requires two --references")
        prepare(args.prepare.resolve(), [p.resolve() for p in args.references])
    elif args.run:
        repeat(args.run.resolve())
    else:
        sample(args.sample.resolve())
