"""Fresh common-admission replay without rewriting frozen geometry records."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

import numpy as np

import microgen_workflow as workflow
from job_parts import worker_part

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

SOURCES = {
    "invalid-gyroid": workflow.HERE / "runs/g3-gyroid-first-pass-fem43-v1",
    "admitted-schwarz": workflow.HERE / "runs/g3-schwarz-tet10-live-amg-repeated-v1/repeat-1/meshers-tet10",
}


def worker(run, name):
    source = SOURCES[name]
    row = json.loads((source / "result.json").read_text())
    archive = np.load(source / "mesh.npz")
    measured = workflow.geometry(source, row["case"], archive["points"], archive["tetrahedra"], archive["surface"])
    for key, value in row["geometry"].items():
        if isinstance(value, float):
            np.testing.assert_allclose(measured[key], value, atol=1e-12, rtol=1e-10)
        else:
            assert measured[key] == value, key
    expected = name == "admitted-schwarz"
    assert measured["oriented_facet_consistency_passed"] is expected
    assert measured["geometry_admitted"] is expected
    workflow.save(run / (name + ".json"), {"source": str(source), "legacy_metrics_reproduced": True,
        "geometry": measured, "frozen_source_unchanged": True})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", choices=tuple(SOURCES))
    args = parser.parse_args()
    if args.worker:
        worker(args.run, args.worker)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__)]
    paths.extend(workflow.HERE / name for name in ("microgen_workflow.py", "facet_incidence_audit.py",
        "test_facet_incidence_audit.py", "test_facet_admission.py", "job_parts.py", "test_job_parts.py",
        "verify_schwarz_tet10_repeated.py", "test_verify_schwarz_tet10_repeated.py"))
    for source in SOURCES.values():
        row = json.loads((source / "result.json").read_text())
        paths.extend(source / name for name in ("result.json", "mesh.npz", row["case"] + "-samples.npz", row["case"] + "-specimen.json"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 120, "measurement_source_ref": "2f0363f",
        "numerical_admission_thresholds_unchanged": True, "frozen_sources_unmodified": True})
    with Job(parts=["controls", *SOURCES, "frozen-verification"]) as job:
        for stage in ("controls", *SOURCES, "frozen-verification"):
            with worker_part(job, stage) as part:
                if stage == "controls":
                    arguments = ["-m", "pytest", "-q", "experiments/adaptive/test_facet_incidence_audit.py",
                        "experiments/adaptive/test_facet_admission.py", "experiments/adaptive/test_job_parts.py",
                        "experiments/adaptive/test_verify_schwarz_tet10_repeated.py"]
                elif stage == "frozen-verification":
                    arguments = [str(workflow.HERE / "verify_schwarz_tet10_repeated.py"), "--run",
                        str(workflow.HERE / "runs/g3-schwarz-tet10-live-amg-repeated-v1"), "--output",
                        str(args.run / "frozen-verification.json"), "--source-ref", "2f0363f"]
                else:
                    arguments = [__file__, "--run", str(args.run), "--worker", stage]
                record = workflow.bounded([sys.executable, *arguments], args.run / (stage + ".log"), part, timeout=120)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                job.result(stage, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
