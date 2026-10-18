"""Check whether retained reference models can share any one-percent tensor."""

from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import microgen_workflow as workflow


def pair_witness(first, second):
    values, vectors = eigh(first, second)
    index = -1 if values[-1] >= 1 / values[0] else 0
    strain = vectors[:, index]
    strain = strain / np.linalg.norm(strain)
    energies = [float(strain @ tensor @ strain) for tensor in (first, second)]
    high = int(energies[1] > energies[0])
    ratio = energies[high] / energies[1 - high]
    lower_bound = (ratio - 1) / (ratio + 1)
    return lower_bound, high, strain, energies


for first_scale, second_scale in ((1.0, 1.0), (1.01, 1.0), (1.0, 1.04), (0.96, 1.0)):
    measured, *_ = pair_witness(first_scale * np.eye(6), second_scale * np.eye(6))
    expected = abs(first_scale - second_scale) / (first_scale + second_scale)
    assert abs(measured - expected) < 1e-14

output = Path(__file__).with_name("result.json")
assert not output.exists()
sources = [ROOT / "experiments/adaptive/runs" / name for name in (
    "g3-lattice-reference-sensitivity-v2.json",
    "g3-lattice-microgen-reference-comparison-v1.json",
    "g3-lattice-finer-reference-comparison-v1.json")]
models = []
for source in sources:
    report = json.loads(source.read_text())
    assert report["case"] == "smooth-lattice"
    for index, model in enumerate(report["models"]):
        tensor = np.asarray(model["stiffness"])
        assert tensor.shape == (6, 6) and np.isfinite(tensor).all()
        assert np.linalg.eigvalsh(tensor).min() > 0
        models.append((tensor, {"source": str(source), "model_index": index,
            "length_proxy": model["length_proxy"], "assumed_order": model["assumed_order"]}))
target = 0.01
worst_bound, worst = -1.0, None
for first, second in combinations(models, 2):
    lower_bound, high, strain, energies = pair_witness(first[0], second[0])
    if lower_bound > worst_bound:
        pair = (first, second)
        intervals = [[(1 - target) * energy, (1 + target) * energy] for energy in energies]
        worst_bound = lower_bound
        worst = {"higher_model": pair[high][1], "lower_model": pair[1 - high][1],
            "engineering_strain_in_saved_tensor_basis": strain.tolist(),
            "higher_quadratic_form": energies[high], "lower_quadratic_form": energies[1 - high],
            "higher_model_candidate_interval": intervals[high], "lower_model_candidate_interval": intervals[1 - high],
            "interval_gap_at_target": intervals[high][0] - intervals[1 - high][1]}
assert worst is not None
incompatible = worst_bound > target + 1e-12
assert (worst["interval_gap_at_target"] > 0) == incompatible
workflow.save(output, {"inputs": {str(source): workflow.digest(source) for source in sources},
    "driver_sha256": workflow.digest(Path(__file__)), "models": len(models), "target": target,
    "necessary_common_relative_error_lower_bound": worst_bound,
    "pairwise_obstruction_to_target": bool(incompatible), "worst_pair_witness": worst,
    "derivation": "A target tensor must satisfy (1-eps)*R <= A <= (1+eps)*R for every reference R. A pair with quadratic-form ratio r requires eps >= (r-1)/(r+1).",
    "scope": "Necessary compatibility condition for retained sensitivity models only; not an optimal global radius, a mesher impossibility theorem, or a continuum error bound"})
