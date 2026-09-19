"""Render the recorded displacement and energy convergence measurements."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter

ROOT = Path(__file__).resolve().parents[2]


def main():
    """Write a standalone scientific figure from the solver's raw JSON results."""
    report = json.loads((ROOT / "docs/solver-convergence-results.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for case in report["cases"]:
        rows = case["refinement"]
        h = np.array([r["maximum_edge_length"] for r in rows])
        for axis, metric in zip(axes, ["l2_displacement_error", "energy_error"], strict=True):
            axis.loglog(h, [r[metric] for r in rows], "o-", label=case["geometry"].title())
    for axis, title, order, anchor in zip(
        axes, ["Displacement L2 error", "Elastic energy error"], [2, 1], [0.003, 0.3], strict=True
    ):
        h = np.array([0.12, 0.23])
        axis.loglog(h, anchor * (h / 0.12) ** order, "--", color="0.55", label=f"Slope {order}")
        axis.set_title(title)
        axis.set_xlabel("Maximum physical edge length")
        axis.set_ylabel("Absolute error")
        axis.grid(visible=True, which="both", alpha=0.2)
        axis.set_xticks([0.12, 0.16, 0.20, 0.24], labels=["0.12", "0.16", "0.20", "0.24"])
        axis.xaxis.set_minor_formatter(NullFormatter())
        axis.legend(fontsize=9)
    fig.savefig(ROOT / "docs/site/assets/solver-convergence.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
