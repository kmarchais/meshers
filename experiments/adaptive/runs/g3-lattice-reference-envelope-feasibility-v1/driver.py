"""Read all retained lattice models without changing any admission reference."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

import microgen_workflow as workflow
from cubic_refinement import MODEL_REPORTS
from job_parts import worker_part
from reference_envelope import assess

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def worker(run):
    models, provenance = [], []
    for filename in MODEL_REPORTS:
        path = workflow.HERE / "runs" / filename
        report = json.loads(path.read_text())
        assert report["case"] == "smooth-lattice" and len(report["models"]) == 8
        assert report["stiffness_target"] == .01
        for index, model in enumerate(report["models"]):
            models.append(model["stiffness"])
            provenance.append({"source": str(path), "index": index, "length_proxy": model["length_proxy"],
                "assumed_order": model["assumed_order"]})
    assert len(models) == 24
    result = assess(models, target=.01)
    result.update(case="smooth-lattice", provenance=provenance, references_unchanged=True,
        admission_threshold_unchanged=True, mesh_candidate_created=False, stiffness_solve_performed=False)
    workflow.save(run / "result.json", result)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args.run)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), *[workflow.HERE / "runs" / name for name in MODEL_REPORTS]]
    paths.extend(workflow.HERE / name for name in ("reference_envelope.py", "test_reference_envelope.py", "cubic_refinement.py",
        "microgen_workflow.py", "job_parts.py", "test_job_parts.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 60, "models": 24, "target": .01,
        "maximum_ordered_pairs": 552, "contradiction_margin": 1e-10, "no_admission_change": True})
    with Job(parts=["controls", "envelope"]) as job:
        for stage in ("controls", "envelope"):
            with worker_part(job, stage) as part:
                arguments = (["-m", "pytest", "-q", "experiments/adaptive/test_reference_envelope.py", "experiments/adaptive/test_job_parts.py"]
                    if stage == "controls" else [__file__, "--run", str(args.run), "--worker"])
                record = workflow.bounded([sys.executable, *arguments], args.run / (stage + ".log"), part, timeout=60)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                job.result(stage, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
