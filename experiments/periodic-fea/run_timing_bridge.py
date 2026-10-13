"""Frozen live matched-node timing and public Microgen application checks."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path.home() / "bin"))
import pueue_progress
from run import digest, write_json
from repeat_recipe_generation import repeat

def main():
    output = Path(__file__).resolve().parents[1]
    manifest = json.loads((output / "manifest.json").read_text())
    report = pueue_progress.report
    pueue_progress.report = lambda fraction, message, **kw: report(.8 * fraction, message, **kw)
    repeat(output)
    import verify_microgen_bridge as bridge
    bridge.report = lambda fraction, message, **kw: report(.8 + .15 * fraction, message, **kw)
    checked = bridge.verify(list(map(Path, manifest["bridge_artifacts"])))
    write_json(output / "microgen-bridge.json", checked)
    import numpy as np
    from microgen.shape import surface_functions
    from microgen_benchmark.cases import Case
    probes = np.random.default_rng(20261013).uniform(-.5, .5, (2048, 3))
    field = Case("schwarz-p", "schwarz-p", "schwarz_p", 1.2, "", "")
    microgen_values = surface_functions.schwarz_p(*(2 * np.pi * probes).T) / .6
    difference = float(np.abs(field.field(*probes.T) - microgen_values).max())
    assert difference < 1e-14
    write_json(output / "microgen-field-equivalence.json", {
        "probes": len(probes), "maximum_absolute_difference": difference,
        "raw_band": [-.6, .6], "microgen_source_sha256": digest(manifest["microgen_surface_functions"]),
        "scope": "Existing Schwarz P analytic fixture versus the installed Microgen callable at fixed width"})
    for n, h in manifest["input_sha256"].items():
        assert digest(n) == h, n
    for n, h in manifest["source_sha256"].items():
        assert digest(output / n) == h, n
    write_json(output / "completion.json", {"state": "complete", "manifest_sha256": digest(output / "manifest.json"),
                                           "input_hashes_unchanged": True, "exact_public_microgen_roundtrips": len(checked["rows"])})
    report(1, "Matched-node live timing and public Microgen checks passed", eta_s=0)

if __name__ == "__main__":
    main()
