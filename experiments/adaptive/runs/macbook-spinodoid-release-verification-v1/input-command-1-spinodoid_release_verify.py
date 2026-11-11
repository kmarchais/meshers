"""Replay selected spinodoid artifacts and verify candidate/baseline provenance."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import microgen_workflow as w


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    a = p.parse_args()
    a.run.mkdir(parents=True, exist_ok=False)
    study = w.HERE/"runs/macbook-spinodoid-release-repeated-v1/study"
    report = w.ROOT/"handoff/macbook/spinodoid-release-comparison-20261111.json"
    output = a.run/"array-verification.json"
    subprocess.run([sys.executable,str(w.HERE/"tpms_tet4_count_verify.py"),
        "--selection",str(w.ROOT/"handoff/macbook/spinodoid-release-final-selection-20261111.json"),
        "--studies",str(study),"--report",str(report),"--output",str(output),
        "--arms","before,native,actual"],check=True)
    arrays = json.loads(output.read_text())
    assert all(n == 1 for n in arrays["repeated_array_set_counts"].values())
    candidate = w.ROOT/".local-work/macbook-spinodoid-release-native-v2"
    receipt = json.loads((candidate/"receipt.json").read_text())
    for name,digest in receipt["inputs"].items():
        assert w.digest(w.ROOT/name) == digest == w.digest(candidate/"sources"/name),name
    binary = next((candidate/"python/meshers").glob("_meshers*.so"))
    assert w.digest(binary) == receipt["native_binary_sha256"]
    baseline = json.loads((w.ROOT/".local-work/macbook-spinodoid-release-baseline-python/receipt.json").read_text())
    assert all(w.digest(Path(path)) == digest for path,digest in baseline["files"].items())
    before_hash = next(digest for path,digest in baseline["files"].items() if path.endswith(".so"))
    rows = json.loads((study/"summary.json").read_text())["rows"]
    for row in rows:
        expected = receipt["native_binary_sha256"] if row["arm"] == "native" else before_hash
        assert row["native_binary_sha256"] == expected
        assert row["generated_field_identity_checked"] == 4096 if row["arm"] != "actual" else True
    provenance = w.ROOT/"handoff/macbook/implicit-tet4-verification-20261111.json"
    old = json.loads(provenance.read_text())
    assert old["passed"] and all(w.digest(Path(path)) == digest for path,digest in old["inputs"].items())
    regression_path = w.HERE/"runs/macbook-spinodoid-release-regression-v1/study/summary.json"
    regression = json.loads(regression_path.read_text())
    assert regression["complete"] and {r["case"] for r in regression["rows"]} == {"gyroid","schwarz-solid","smooth-lattice"}
    assert all(r["qualified"] and r["native_binary_sha256"] == receipt["native_binary_sha256"] for r in regression["rows"])
    tests = w.HERE/"runs/macbook-spinodoid-release-python-tests-v1"
    assert json.loads((tests/"worker.json").read_text())["returncode"] == 0
    w.save(a.run/"verification.json",dict(passed=True, arrays={k:v for k,v in arrays.items()
        if k not in ("artifacts","distinct_meshes")},
        candidate_binary_sha256=receipt["native_binary_sha256"], baseline_binary_sha256=before_hash,
        build_receipt_sha256=w.digest(candidate/"receipt.json"), field_provenance_sha256=w.digest(provenance),
        regression_sha256=w.digest(regression_path), python_tests_log_sha256=w.digest(tests/"worker.log"),
        report_sha256=w.digest(report), source_sha256=w.digest(Path(__file__)),
        scope="Native source snapshot, old/new binaries, actual Tet4 arrays, frozen physical inputs, replayed metrics and original-case admission. Mac arm64 only; no FEA or cross-platform publication certification."))
    print("Verified spinodoid before/after artifacts, binary provenance and original-study controls")


if __name__ == "__main__":
    main()
