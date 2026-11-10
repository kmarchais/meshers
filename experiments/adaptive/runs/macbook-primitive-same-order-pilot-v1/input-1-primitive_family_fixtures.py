"""Five frozen basic solids with analytic volumes and actual Microgen fields."""

import argparse
import json
from pathlib import Path

import microgen_workflow as w
import numpy as np
from scipy.stats import qmc

import scorecard


DEFINITIONS = [
    dict(api="Box", parameters=dict(dim=[0.6, 0.5, 0.4]), volume=0.12),
    dict(api="Capsule", parameters=dict(height=0.4, radius=0.15),
         volume=np.pi * 0.15**2 * 0.4 + 4 * np.pi * 0.15**3 / 3),
    dict(api="Cylinder", parameters=dict(height=0.6, radius=0.23),
         volume=np.pi * 0.23**2 * 0.6),
    dict(api="Ellipsoid", parameters=dict(radii=[0.35, 0.25, 0.18]),
         volume=4 * np.pi * 0.35 * 0.25 * 0.18 / 3),
    dict(api="Sphere", parameters=dict(radius=0.35), volume=4 * np.pi * 0.35**3 / 3),
]


def definition(name, root):
    return json.loads((root / name / (name + "-definition.json")).read_text())


def actual_shape(name, root, resolution=None):
    import microgen

    spec = definition(name, root)
    return getattr(microgen, spec["api"])(**spec["parameters"]), "solid"


def case(name, root):
    spec = definition(name, root)
    p = spec["parameters"]
    if spec["api"] == "Sphere":
        def field(x, y, z):
            return np.sqrt(x*x + y*y + z*z) - p["radius"]
    elif spec["api"] == "Ellipsoid":
        a, b, c = p["radii"]
        def field(x, y, z):
            return np.sqrt((x/a)**2 + (y/b)**2 + (z/c)**2) - 1
    elif spec["api"] == "Capsule":
        def field(x, y, z):
            dx = np.maximum(np.abs(x) - p["height"]/2, 0)
            return np.sqrt(dx*dx + y*y + z*z) - p["radius"]
    elif spec["api"] == "Cylinder":
        def field(x, y, z):
            radial = np.sqrt(y*y + z*z) - p["radius"]
            axial = np.abs(x) - p["height"]/2
            return (np.sqrt(np.maximum(radial, 0)**2 + np.maximum(axial, 0)**2)
                    + np.minimum(np.maximum(radial, axial), 0))
    else:
        assert spec["api"] == "Box"
        a, b, c = np.array(p["dim"]) / 2
        def field(x, y, z):
            u, v, t = np.abs(x)-a, np.abs(y)-b, np.abs(z)-c
            return (np.sqrt(np.maximum(u, 0)**2 + np.maximum(v, 0)**2 + np.maximum(t, 0)**2)
                    + np.minimum(np.maximum(u, np.maximum(v, t)), 0))
    return scorecard.Case(name, field, (-np.inf, 0), spec["volume"])


def retained_field(name, root):
    """Exact same negative sublevel set, with explicit intersecting branches."""
    spec = definition(name, root)
    p = spec["parameters"]
    if spec["api"] == "Box":
        a, b, c = np.array(p["dim"]) / 2
        def field(x, y, z):
            return np.maximum(np.abs(x)-a, np.maximum(np.abs(y)-b, np.abs(z)-c))
    else:
        assert spec["api"] == "Cylinder"
        def field(x, y, z):
            return np.maximum(np.sqrt(y*y+z*z)-p["radius"], np.abs(x)-p["height"]/2)
    return field


def prepare(root):
    root.mkdir(parents=True, exist_ok=False)
    entries = []
    sources = {str(Path(__file__)): w.digest(Path(__file__))}
    for item in DEFINITIONS:
        name = "primitive-" + item["api"].lower()
        spec = dict(item, name=name, part="solid", surface_function=item["api"], periodic=False)
        folder = root / name
        folder.mkdir()
        w.save(folder / (name + "-definition.json"), spec)
        target = case(name, root)
        actual, _ = actual_shape(name, root)
        import inspect
        source = Path(inspect.getfile(type(actual)))
        sources[str(source)] = w.digest(source)
        points = qmc.Sobol(3, scramble=True, seed=20261110).random_base2(12) - 0.5
        np.testing.assert_allclose(target.function(*points.T), actual._func(*points.T),
                                   rtol=1e-12, atol=1e-12)
        walls = scorecard.wall_samples(target, count=2**15)
        # Sample exact sharp edges too; random wall probes can miss them.
        if item["api"] == "Box":
            half = np.array(item["parameters"]["dim"]) / 2
            extra = []
            for axis in range(3):
                other = [i for i in range(3) if i != axis]
                for u in (-1, 1):
                    for v in (-1, 1):
                        edge = np.zeros((257, 3))
                        edge[:, axis] = np.linspace(-half[axis], half[axis], 257)
                        edge[:, other] = [u*half[other[0]], v*half[other[1]]]
                        extra.append(edge)
            walls = np.vstack((walls, *extra))
        if item["api"] == "Cylinder":
            angle = np.linspace(0, 2*np.pi, 4096, endpoint=False)
            circle = np.column_stack((np.zeros_like(angle), item["parameters"]["radius"] * np.cos(angle),
                                      item["parameters"]["radius"] * np.sin(angle)))
            circles = []
            for sign in (-1, 1):
                ring = circle.copy()
                ring[:, 0] = sign * item["parameters"]["height"] / 2
                circles.append(ring)
            walls = np.vstack((walls, *circles))
        assert len(walls) >= 32768 and np.max(abs(target.function(*walls.T))) < 1e-9
        np.savez_compressed(folder / (name + "-samples.npz"), walls=walls, caps=np.empty((0, 3)))
        w.save(folder / (name + "-specimen.json"), dict(
            definition=spec, volume=item["volume"], volume_reference="analytic", volume_scramble_spread=0,
            wall_samples=len(walls), field_identity_points=4096, periodic_axes=[False]*3))
        entries.append(spec)
    w.save(root / "protocol.json", dict(
        sources=sources, unit_domain=[[-0.5]*3, [0.5]*3], cases=entries,
        actual_bounds="Public shape-specific bounds, enclosing the same complete solid",
        scope="Five implicit basic solids; no CAD route or arbitrary pose claim. Original study unchanged."))
    w.save(root / "result.json", dict(passed=True, cases=entries,
        inputs={str(f): w.digest(f) for f in root.rglob("*") if f.is_file()}))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--run", type=Path, required=True)
    a = p.parse_args()
    w.imports()
    prepare(a.run)
