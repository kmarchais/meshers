"""Audit an archived seed or guide before the bounded body-refinement prototype."""
import json
import os
import sys
from pathlib import Path

for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[name] = "1"
sys.path.insert(0, str(Path.home() / "bin"))
import pueue_progress
from run import digest, run, write_json
from verify import verify_run
from body_refine import verify_preserved_boundary


def main():
    output = Path(__file__).resolve().parents[1]
    manifest = json.loads((output / "manifest.json").read_text())
    report = pueue_progress.report
    report(0, "Verify frozen body-refinement inputs", eta_s=60)
    for name, expected in manifest["source_sha256"].items():
        assert digest(output / "source" / name) == expected, name
    for name, expected in manifest["preflight_input_sha256"].items():
        assert digest(name) == expected, name
    input_run = Path(manifest.get("preflight_audit_run") or manifest["guide_audit_run"])
    audit_path = input_run / "artifact-verification.json"
    audit_hash = digest(audit_path)
    audit = verify_run(input_run, write_report=False)
    assert audit["verified_rows"] == 1 and audit["rows"][0]["fields"] == 6
    assert digest(audit_path) == audit_hash
    write_json(output / "guide-readonly-revalidation.json", {
        "audit": audit, "unchanged_audit_sha256": audit_hash,
        "manifest_sha256": digest(output / "manifest.json"),
    })
    pueue_progress.report = lambda fraction, message, **kw: report(.04 + .94 * fraction, message, **kw)
    run(output)
    row = json.loads((output / "results.json").read_text())["rows"][0]
    if row["status"] != "ok" or not all(row[k] for k in ["geometry_eligible", "quality_eligible", "periodic_eligible"]):
        raise RuntimeError(f"Guided prototype failed admission: {row.get('error', row['status'])}")
    import numpy as np
    with np.load(output / "source/body-seed/mesh.npz") as seed, np.load(output / row["artifact"] / "mesh.npz") as body:
        boundary = verify_preserved_boundary(seed["points"], seed["surface"], body["points"], body["surface"])
    for name, expected in manifest["preflight_input_sha256"].items():
        assert digest(name) == expected, name
    write_json(output / "guided-completion.json", {
        "state": "complete", "manifest_sha256": digest(output / "manifest.json"),
        "retained_boundary": boundary, "original_inputs_unchanged": True,
        "scope": "Python adaptation of archived inputs; no generation speed claim",
    })
    report(1, "Body prototype admitted with retained boundary and input audit unchanged", eta_s=0)


if __name__ == "__main__":
    main()
