"""Freeze TPMS thicknesses and a seeded spectrum from Microgen's spinodoid source."""

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
from microgen.shape import surface_functions
from scipy.stats import qmc

from fixtures import TPMS
from run import digest, write_json


def prepare(output, microgen_source):
    output.mkdir(parents=True, exist_ok=False)
    samples = qmc.Sobol(3, scramble=True, seed=20261013).random_base2(16) - .5
    specimens = {"solid-cell": {"function": "solid", "offset": 1}}
    formula_path = Path(surface_functions.__file__)
    for name in TPMS:
        raw = getattr(surface_functions, name)(*(2 * np.pi * samples.T))
        half_width = float(np.quantile(np.abs(raw), .35))
        specimens[f"periodic-{name.replace('_', '-')}"] = {
            "function": name, "offset": 2 * half_width, "microgen_callable": True,
            "sampled_density_target": .35,
        }
    spectral_path = microgen_source / "microgen/shape/_frep_grf.py"
    spec = importlib.util.spec_from_file_location("frozen_microgen_grf", spectral_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    field = module._FrepGRF.from_fft_grf(k0=4 * np.pi, bandwidth=3,
                                       anisotropy=(1, 1, 1), cell_size=1,
                                       resolution=16, seed=20261010, mode_amplitude_threshold=.05)
    spectrum = output / "spinodoid-spectrum.npz"
    np.savez_compressed(spectrum, modes=field.modes, amplitudes=field.amplitudes,
                        phases=field.phases, threshold=field.threshold)
    specimens["periodic-spinodoid"] = {"function": "spinodoid_spectral", "offset": 1,
                                       "microgen_callable": True, "solid": "raw field >= 0"}
    document = {
        "specimens": specimens, "tpms_formula_source": str(formula_path),
        "tpms_formula_sha256": digest(formula_path), "spinodoid_source": str(spectral_path),
        "spinodoid_source_sha256": digest(spectral_path), "spinodoid_spectrum": spectrum.name,
        "spinodoid_spectrum_sha256": digest(spectrum), "spinodoid_mode_count": len(field.modes),
        "spinodoid_parameters": {"k0": float(4 * np.pi), "bandwidth": 3, "anisotropy": [1, 1, 1],
                                  "spectral_resolution": 16, "seed": 20261010, "amplitude_threshold": .05},
        "scope": "Fixture catalogue only; no meshing, FEA qualification or superiority claim",
    }
    write_json(output / "specimens.json", document)
    print(f"Frozen {len(TPMS)} TPMS fields and one {len(field.modes)}-mode spinodoid", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--microgen-source", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.output.resolve(), args.microgen_source.resolve())
