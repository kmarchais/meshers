"""Separate element-mass and density ranking on immutable recovery indicators."""

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import microgen_workflow as workflow

source = ROOT / "experiments/adaptive/runs/g3-lattice-recovery-localization-v1"
output = Path(__file__).with_name("result.json")
assert not output.exists()
data = np.load(source / "indicators.npz")
mass, volumes = data["mass"], data["volumes"]
previous = json.loads((source / "result.json").read_text())
assert mass.shape == volumes.shape and len(mass) == previous["elements"]
assert np.isfinite(mass).all() and (mass >= 0).all() and (volumes > 0).all()
np.testing.assert_allclose(mass.sum(), previous["recovery_trace_indicator"], atol=1e-14, rtol=1e-14)
rows = []
for ranking, values in (("element_mass", mass), ("density", mass / volumes), ("volume", volumes)):
    ordered = np.argsort(values)
    for fraction in (0.05, 0.1, 0.2, 0.5):
        selected = ordered[-int(np.ceil(fraction * len(mass))):]
        rows.append({"ranking": ranking, "element_fraction": fraction,
            "volume_fraction": float(volumes[selected].sum() / volumes.sum()),
            "recovery_mass_fraction": float(mass[selected].sum() / mass.sum())})
workflow.save(output, {"inputs": {str(source / name): workflow.digest(source / name)
    for name in ("result.json", "indicators.npz")}, "driver_sha256": workflow.digest(Path(__file__)),
    "rows": rows, "scope": "Saved recovery-indicator concentration only, not observed enrichment or continuum error"})
