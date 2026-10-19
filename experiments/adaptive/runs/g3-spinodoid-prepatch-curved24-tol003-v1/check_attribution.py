"""Verify native-binary attribution without treating the observer source as baseline source."""

import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import microgen_workflow as workflow

run = Path(__file__).parent
candidate = ROOT / "experiments/adaptive/runs/g3-spinodoid-partial-settle-curved24-tol003-solve-v1"
baseline_package = ROOT / ".local-work/g3-builds/baseline-090a71c/crates/meshers-python/python/meshers"
paths = [Path(__file__), run / "manifest.json", run / "result.json",
         candidate / "manifest.json", candidate / "result.json",
         baseline_package / "_meshers.abi3.so"]
for filename in ("__init__.py", "expression.py"):
    old = baseline_package / filename
    current = ROOT / "crates/meshers-python/python/meshers" / filename
    assert workflow.digest(old) == workflow.digest(current)
    paths.extend((old, current))
inputs = {str(path): workflow.digest(path) for path in paths}
baseline_manifest = json.loads((run / "manifest.json").read_text())
new_manifest = json.loads((candidate / "manifest.json").read_text())
baseline = json.loads((run / "result.json").read_text())
new = json.loads((candidate / "result.json").read_text())
expected_binary = "d41ccb980fbf3339a86e8f9a195e675abaa3910be67de8589dca138fd401b700"
assert baseline_manifest["native_sha256"] == baseline_manifest["native_sha256_after"] == expected_binary
assert workflow.digest(baseline_package / "_meshers.abi3.so") == expected_binary
assert new_manifest["native_sha256"] == "908e326211792142f83f854ad2696cd879d927bd81deacf38eb6a0cb63bc4dc1"
assert not baseline["qualified"] and new["qualified"]
for key in ("case", "cells", "order", "min_length", "native_tolerance", "optimization_passes", "adapt_tolerance"):
    assert baseline[key] == new[key]
source_commit = subprocess.check_output(["git", "rev-parse", "090a71c"], cwd=ROOT, text=True).strip()
source = subprocess.check_output(["git", "show", "090a71c:crates/meshers-core/src/tet_adapt.rs"], cwd=ROOT)
source_hash = hashlib.sha256(source).hexdigest()
assert source_hash == "0c698bd89e31ea692b57971bf90082b5ee61c1fff98d2bccd8b096f9487d1f82"
output = run / "attribution.json"
assert not output.exists()
workflow.save(output, {"inputs": inputs, "baseline_source_commit": source_commit,
    "baseline_tet_adapt_sha256": source_hash,
    "observer_tet_adapt_sha256": baseline_manifest["tet_adapt_sha256"],
    "observer_source_note": "The runner records current workspace source, not the archived binary's source",
    "baseline_distance": baseline["geometry"]["sampled_two_sided_max"],
    "candidate_distance": new["geometry"]["sampled_two_sided_max"],
    "scope": "Matched recipe geometry attribution; no FE or repeated timing claim for rejected baseline"})
assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())
print(json.dumps({"baseline_admitted": baseline["qualified"], "candidate_admitted": new["qualified"]}))
