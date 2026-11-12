"""Frozen representatives of all pinned Microgen TPMS fields and solid parts."""

import argparse
import json
from pathlib import Path

import microgen_workflow as w
import numpy as np
from scipy.stats import qmc

from microgen_catalog import catalog
import scorecard


def field(spec):
    from microgen.shape import surface_functions

    function = getattr(surface_functions, spec["surface_function"])
    half = spec["offset"] / 2

    def raw(x, y, z):
        return function(2 * np.pi * x, 2 * np.pi * y, 2 * np.pi * z)

    if spec["part"] == "sheet":
        return scorecard.Case(spec["name"], raw, (-half, half))

    def solid(x, y, z):
        return half + (raw(x, y, z) if spec["part"] == "lower skeletal" else -raw(x, y, z))

    return scorecard.Case(spec["name"], solid, (-np.inf, 0.0))


def case(name, root):
    return field(json.loads((root / name / (name + "-definition.json")).read_text()))


def actual_shape(name, root, resolution=None):
    import microgen
    from microgen.shape import surface_functions

    spec = json.loads((root / name / (name + "-definition.json")).read_text())
    options = {} if resolution is None else dict(resolution=resolution)
    return microgen.Tpms(
        getattr(surface_functions, spec["surface_function"]), offset=spec["offset"], **options
    ), spec["part"]


def prepare(root):
    import microgen.shape.surface_functions as functions

    inventory = catalog()
    root.mkdir(parents=True, exist_ok=False)
    w.save(
        root / "protocol.json",
        dict(
            catalog=inventory,
            target_volume_fraction=0.30,
            calibration_seed=20261110,
            calibration_points=2**18,
            independent_volume_scrambles=4,
            points_per_scramble=2**20,
            wall_projection_starts=2**15,
            cap_candidates_per_face=2**12,
            rule="One Sobol quantile fixes each raw-field offset once; both arms use it at every resolution. No per-mesh density recalibration. Original study unchanged.",
            mesh_policy="First screen: native uniform16 Tet10 with explicit optional fit; actual public default Tet4, then read failures before modified routes. Counts and order remain explicit; no same-order speed claim from this screen.",
            caps=dict(nodes=w.CAP, tetrahedra=w.CAP),
            distance=0.002,
            relative_volume=0.01,
            minimum_quality=0.1,
            source_hashes={
                str(f): w.digest(f)
                for f in (Path(__file__), Path(functions.__file__), w.HERE / "scorecard.py")
            },
        ),
    )
    definitions = []
    for entry in inventory["rows"]:
        if entry["family"] != "periodic fields":
            continue
        function_name = entry["api"].rsplit(".", 1)[1]
        function = getattr(functions, function_name)
        part = entry["part"]
        name = "tpms-" + function_name.replace("_", "-") + "-" + part.split()[0]
        probes = qmc.Sobol(3, scramble=True, seed=20261110).random_base2(18) - 0.5
        raw = function(*(2 * np.pi * probes).T)
        offset = float(
            2 * np.quantile(abs(raw), 0.30)
            if part == "sheet"
            else -2 * np.quantile(raw, 0.30)
            if part == "lower skeletal"
            else 2 * np.quantile(raw, 0.70)
        )
        spec = dict(name=name, surface_function=function_name, part=part, offset=offset)
        directory = root / name
        directory.mkdir()
        w.save(directory / (name + "-definition.json"), spec)
        target = field(spec)
        shape, _ = actual_shape(name, root)
        probes = np.random.default_rng(20261110).uniform(-0.5, 0.5, (4096, 3))
        raw = shape.raw_field(*probes.T)
        expected = (
            raw if part == "sheet" else offset / 2 + (raw if part == "lower skeletal" else -raw)
        )
        np.testing.assert_allclose(target.function(*probes.T), expected, atol=1e-12, rtol=1e-12)
        periodic_errors = []
        for axis in range(3):
            lo, hi = probes.copy(), probes.copy()
            lo[:, axis], hi[:, axis] = -0.5, 0.5
            periodic_errors.append(
                float(abs(target.function(*lo.T) - target.function(*hi.T)).max())
            )
        assert max(periodic_errors) < 1e-10
        estimates = []
        for seed in range(4):
            sobol = qmc.Sobol(3, scramble=True, seed=seed)
            count = 0
            for _ in range(4):
                values = target.function(*(sobol.random(2**18) - 0.5).T)
                count += int(((values >= target.levels[0]) & (values <= target.levels[1])).sum())
            estimates.append(count / 2**20)
        walls = scorecard.wall_samples(target, count=2**15)
        assert len(walls) >= 4096, name
        uv = qmc.Sobol(2, scramble=True, seed=119).random_base2(12) - 0.5
        caps = []
        for axis in range(3):
            for side in (-0.5, 0.5):
                points = np.empty((len(uv), 3))
                points[:, axis] = side
                points[:, [i for i in range(3) if i != axis]] = uv
                values = target.function(*points.T)
                caps.append(points[(values >= target.levels[0]) & (values <= target.levels[1])])
        np.savez_compressed(
            directory / (name + "-samples.npz"), walls=walls, caps=np.concatenate(caps)
        )
        w.save(
            directory / (name + "-specimen.json"),
            dict(
                volume=float(np.mean(estimates)),
                volume_sobol_estimates=estimates,
                volume_scramble_spread=float(np.ptp(estimates)),
                volume_samples_per_scramble=2**20,
                wall_samples=len(walls),
                definition=spec,
                periodic_axes=[True] * 3,
                periodic_field_errors=periodic_errors,
                field_identity_points=4096,
            ),
        )
        definitions.append(spec)
        print(
            f"Prepared {len(definitions)}/48: {name}, volume {np.mean(estimates):.6f}", flush=True
        )
    w.save(
        root / "result.json",
        dict(
            passed=True,
            cases=definitions,
            inputs={str(f): w.digest(f) for f in root.rglob("*") if f.is_file()},
            original_study_unchanged=True,
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--run", type=Path, required=True)
    a = p.parse_args()
    w.imports()
    prepare(a.run)
