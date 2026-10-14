"""Revalidate the original archive read-only, then run one bounded FEA replay."""
import json
import os
from pathlib import Path
import sys

for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[name] = "1"
sys.path.insert(0, str(Path.home() / "bin"))
import pueue_progress
from run import digest, run, write_json
from verify import verify_run


def main():
    output = Path(__file__).resolve().parents[1]
    manifest = json.loads((output / "manifest.json").read_text())
    report = pueue_progress.report
    report(0, "Check frozen geometry and FEA inputs", eta_s=1800)

    def hashes():
        for name, expected in manifest["source_sha256"].items():
            assert digest(output / "source" / name) == expected, name
        for name, expected in manifest["preflight_input_sha256"].items():
            assert digest(name) == expected, name

    hashes()
    provenance = json.loads((output / "source/geometry-seed/input-provenance.json").read_text())
    original = Path(provenance["original_seed_run"])
    audit_path = original / "artifact-verification.json"
    audit_hash = digest(audit_path)
    audit = verify_run(original, write_report=False)
    assert audit["verified_rows"] == 1 and audit["rows"][0]["fields"] == 6
    assert digest(audit_path) == audit_hash
    write_json(output / "readonly-revalidation.json", {"audit": audit, "unchanged_audit_sha256": audit_hash})
    pueue_progress.report = lambda fraction, message, **kw: report(.04 + .94*fraction, message, **kw)
    run(output)
    document = json.loads((output / "results.json").read_text())
    row = document["rows"][0]
    if row["status"] != "ok" or not all(row[k] for k in ["geometry_eligible", "quality_eligible", "periodic_eligible"]):
        raise RuntimeError(f"Geometry FEA replay failed: {row.get('error', row['status'])}")
    assert document["artifact_verification"]["verified_rows"] == 1
    assert row["geometry_replay"]["exact_qualified_mesh_arrays"]
    hashes()
    write_json(output / "completion.json", {
        "state": "complete", "manifest_sha256": digest(output / "manifest.json"),
        "original_inputs_unchanged": True, "geometry_replay": row["geometry_replay"],
        "scope": "Six-load FEA qualification of saved coarsened geometry; accuracy comparison pending",
    })
    report(1, "Coarsened geometry FEA and saved-field audit complete", eta_s=0)


if __name__ == "__main__":
    main()
