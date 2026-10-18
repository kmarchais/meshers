"""Screen a cut-biased edge median and select a volume-based activation scale."""

import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import microgen_workflow as workflow

root = ROOT / "experiments/adaptive/runs"
control = root / "g3-lattice-sizing-minbound-v1/g3-lattice-sizing-minbound-v1-uniform/result.json"
edge_report = root / "g3-lattice-sizing-edge-scale-v1/result.json"
output = Path(__file__).with_name("result.json")
assert not output.exists()
row = json.loads(control.read_text())
assert row["qualified"] and row["geometry"]["geometry_admitted"]
median_edge = json.loads(edge_report.read_text())["selected_nominal_length"]
volume, elements = row["geometry"]["volume"], row["elements"]
equivalent_edge = (12 * volume / (math.sqrt(2) * elements))**(1 / 3)
regular_volume_at_split_threshold = math.sqrt(2) / 12 * (1.33 * median_edge)**3
workflow.save(output, {"inputs": {str(path): workflow.digest(path) for path in (control, edge_report)},
    "driver_sha256": workflow.digest(Path(__file__)), "unweighted_edge_median": median_edge,
    "regular_tet_count_estimate_at_unweighted_scale": volume / regular_volume_at_split_threshold,
    "tetrahedra_cap": workflow.CAP, "selected_nominal_length": equivalent_edge,
    "selection_rule": "Use the regular-tetrahedron edge corresponding to the admitted control's mean element volume; do not use a cut-biased count median",
    "scope": "Input-feasibility screening and one preselected scale only; regular-cell count is an estimate, not a guarantee about the adaptive algorithm or stiffness"})
