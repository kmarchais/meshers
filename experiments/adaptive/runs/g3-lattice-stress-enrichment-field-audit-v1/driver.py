"""Check saved enrichment-control periodic jumps before selecting a live pilot."""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import microgen_workflow as workflow
from stress_recovery import periodic_displacement_jump_error
import stress_recovery

np, *_ = workflow.imports()
source = ROOT / "experiments/adaptive/runs/g3-lattice-stress-enrichment-v1"
output = Path(__file__).with_name("result.json")
assert not output.exists()
row = json.loads((source / "smooth-lattice.json").read_text())
workers = json.loads((source / "workers.json").read_text())
assert workers["smooth-lattice"]["returncode"] == 0
assert row["energy_identity_relative_error"] < 1e-7
saved = np.load(source / "smooth-lattice-recovery.npz")
checks = {}
inputs = [source / name for name in ("workers.json", "smooth-lattice.json", "smooth-lattice-recovery.npz")]
for route in ("linear", "straight"):
    path = source / ("smooth-lattice-control-" + route + ".npz")
    inputs.append(path)
    mesh = np.load(path)
    jump = periodic_displacement_jump_error(mesh["points"], mesh["surface"], saved[route + "_displacement"])
    assert jump < 1e-9
    recorded = row["field_consistency"][route]
    assert max(recorded["reaction_relative_error"], recorded["hill_mandel_relative_error"]) < 1e-7
    assert len(row["solver_records"][route]) == 6
    assert max(record["relative_residual"] for record in row["solver_records"][route]) < 2e-9
    assert np.linalg.eigvalsh(np.asarray(row["stiffness"][route])).min() > 0
    checks[route] = {"periodic_jump_max_recomputed": jump,
        "recorded_worker_field_consistency": recorded}
rankings = row["top_20_percent_elements_enrichment_energy_capture"]
selected = rankings["recovery"] > max(rankings["volume"], rankings["stress-magnitude"]) and row["density_spearman"] > 0.5
assert bool(selected) == row["live_sizing_pilot_worth_testing"]
workflow.save(output, {"inputs": {str(path): workflow.digest(path) for path in inputs},
    "driver_sha256": workflow.digest(Path(__file__)),
    "helper_sha256": workflow.digest(Path(stress_recovery.__file__)), "checks": checks,
    "ranking_selector_passed": bool(selected), "field_checks_passed": True, "live_pilot_eligible": bool(selected),
    "scope": "Recomputed saved-field jumps and completed-worker stress/energy checks; diagnostic controls are not continuum references or benchmark rankings"})
