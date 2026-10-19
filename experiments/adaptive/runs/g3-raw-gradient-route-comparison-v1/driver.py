"""Compare raw-value lattice routes before and after the exact-tie derivative fix."""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job
import microgen_workflow as workflow

run = Path(__file__).parent
assert not (run / "result.json").exists()
baseline = ROOT / ".local-work/g3-builds/baseline-1ac72d9/crates/meshers-python/python"
original_expression = baseline / "meshers/expression.py"
current_expression = ROOT / "crates/meshers-python/python/meshers/expression.py"
import meshers._meshers
sources = [Path(__file__), original_expression, current_expression,
    baseline / "meshers/_meshers.abi3.so", Path(meshers._meshers.__file__),
    ROOT / "experiments/adaptive/admission_pilot.py", Path(workflow.__file__)]
hashes = {str(path): workflow.digest(path) for path in sources}
native_hash = workflow.digest(Path(meshers._meshers.__file__))
assert workflow.digest(baseline / "meshers/_meshers.abi3.so") == native_hash
workflow.save(run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
    "worker_timeout_s": 300, "threads": 1, "native_sha256": native_hash,
    "baseline_expression_sha256": workflow.digest(original_expression),
    "candidate_expression_sha256": workflow.digest(current_expression),
    "provenance_note": "Inner observer-source paths refer to the current checkout; these explicit package checks distinguish the actual loaded derivative policy",
    "scope": "Three geometry-only controls; no algebraic symmetric-minimum rewrite, FE or timing victory"})
arms = [("before-clamped", True, True), ("after-clamped", False, True), ("after-fully-raw", False, False)]
rows = []
with Job(parts=[name for name, _, _ in arms]) as job:
    for name, before, clamped in arms:
        expression = original_expression if before else current_expression
        bootstrap = ("import sys,hashlib;from pathlib import Path;"
            f"sys.path.insert(0,{str(ROOT / 'experiments/adaptive')!r});"
            "import meshers.expression as expression;"
            f"assert hashlib.sha256(Path(expression.__file__).read_bytes()).hexdigest()=={workflow.digest(expression)!r};"
            "import admission_pilot;admission_pilot.main()")
        progress = str(run / (name + "-progress.json"))
        command = ["env", "OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1", "RAYON_NUM_THREADS=1",
            "PUEUE_PROGRESS_FILE=" + progress, "GPU_PROGRESS_FILE=" + progress]
        if before:
            command.append("PYTHONPATH=" + str(baseline))
        command.extend([sys.executable, "-c", bootstrap, "--run", str(run / name),
            "--case", "smooth-lattice", "--cells", "24", "--min-length", ".002",
            "--order", "2", "--native-tolerance", "1", "--adapt-tolerance", ".004"])
        if clamped:
            command.extend(["--clip-void", ".02"])
        with job.part(name) as part:
            record = workflow.bounded(command, run / (name + ".log"), part, timeout=300)
            workflow.save(run / (name + "-worker.json"), record)
            assert record["returncode"] == 0, f"Read {name}.log"
            manifest = json.loads((run / name / "manifest.json").read_text())
            assert manifest["native_sha256"] == native_hash
            row = json.loads((run / name / "result.json").read_text())
            rows.append({"arm": name, "actual_expression_sha256": workflow.digest(expression),
                "qualified": row["qualified"], "status": row["status"], "nodes": row.get("nodes"),
                "elements": row.get("elements"), "geometry": row.get("geometry"), "error": row.get("error")})
            workflow.save(run / "result.json", {"arms": rows,
                "scope": "Actual loaded derivative-policy comparison; scientific admission unchanged; no full-tensor or repeated timing claim"})
            job.result(name, qualified=int(row["qualified"]), nodes=row.get("nodes", 0), seconds=record["wall_s"])
assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())
