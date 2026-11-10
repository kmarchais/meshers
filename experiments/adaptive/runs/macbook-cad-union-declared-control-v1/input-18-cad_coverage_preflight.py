"""Inspect one exact CAD lattice fixture without generating a mesh.

Run through macbook_stage.py in the shared queue. CAD availability is reported
separately from completion, and this preflight never claims mesh admission.
"""

import argparse
import importlib.metadata
import importlib.util
import inspect
from pathlib import Path

import microgen_workflow as workflow
import numpy as np


RADIUS = 0.15
LENGTH = 1.0
BOUNDS = (-0.5, 0.5) * 3
VERTICES = [[-0.5, 0, 0], [0.5, 0, 0], [0, -0.5, 0], [0, 0.5, 0]]
PAIRS = [[0, 1], [2, 3]]


def field(x, y, z):
    """Negative sublevel set is exactly the union of two flat cylinders."""
    along_x = np.maximum(np.sqrt(y * y + z * z) - RADIUS, np.abs(x) - LENGTH / 2)
    along_y = np.maximum(np.sqrt(x * x + z * z) - RADIUS, np.abs(y) - LENGTH / 2)
    return np.minimum(along_x, along_y)


def periodic_field(x, y, z):
    """Same material in BOUNDS, with negative traces inside material caps."""
    return np.minimum(np.sqrt(y * y + z * z) - RADIUS,
                      np.sqrt(x * x + z * z) - RADIUS)


def make_shape():
    import microgen

    return microgen.CustomLattice(
        base_vertices=np.array(VERTICES, dtype=float),
        strut_vertex_pairs=np.array(PAIRS, dtype=int),
        strut_radius=RADIUS,
        strut_heights=LENGTH,
        cell_size=LENGTH,
        strut_joints=False,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    args.run.mkdir(parents=True, exist_ok=False)
    import gmsh
    import meshers
    import microgen
    import microgen.cad
    from meshers.intersection import _periodic_affine, _retained_constraints

    shape = make_shape()
    # Compare independent membership inequalities to the installed public
    # Cylinder fields with exactly the same centers, lengths and rotations.
    rng = np.random.default_rng(20261110)
    points = rng.uniform(-0.6, 0.6, size=(4096, 3))
    fields = []
    for center, orientation, height in zip(
        shape.strut_centers, shape.strut_rotations, shape.strut_heights, strict=True
    ):
        cylinder = microgen.Cylinder(
            center=tuple(center), orientation=orientation, height=height, radius=RADIUS
        )
        fields.append(cylinder._func(*points.T))
    public_union = np.minimum.reduce(fields)
    expected = field(*points.T)
    assert np.array_equal(public_union < 0, expected < 0)
    clipped_periodic = np.maximum(
        periodic_field(*points.T), np.max(np.abs(points), axis=1) - 0.5
    )
    assert np.array_equal(clipped_periodic < 0, expected < 0)

    # The two lateral walls meet on x = +/-y, x² + z² = r².
    # These samples remain separate from every existing specimen and sample set.
    angle = np.linspace(0, 2 * np.pi, 2048, endpoint=False)
    ridges = np.vstack(
        [np.column_stack((RADIUS * np.cos(angle), sign * RADIUS * np.cos(angle),
                          RADIUS * np.sin(angle))) for sign in (-1, 1)]
    )
    assert np.max(np.abs(field(*ridges.T))) < 1e-14
    cap_points = rng.uniform(-0.5, 0.5, size=(1024, 3))
    for axis in range(3):
        low = cap_points.copy()
        high = cap_points.copy()
        low[:, axis] = -0.5
        high[:, axis] = 0.5
        np.testing.assert_array_equal(periodic_field(*low.T), periodic_field(*high.T))
    samples = args.run / "field-and-ridge-probes.npz"
    np.savez_compressed(samples, identity_points=points, ridge_points=ridges,
                        cap_probe_points=cap_points)

    compiled = meshers.compile_field(periodic_field)
    retained = _retained_constraints(compiled)
    analytic_volume = 2 * np.pi * RADIUS**2 * LENGTH - 16 * RADIUS**3 / 3
    result = {
        "preflight_completed": True,
        "mesh_admission": "not tested",
        "fixture": {
            "api": "microgen.CustomLattice", "base_vertices": VERTICES,
            "strut_vertex_pairs": PAIRS, "strut_radius": RADIUS,
            "strut_heights": LENGTH, "cell_size": LENGTH, "strut_joints": False,
            "bounds": list(BOUNDS), "periodic_axes": [True, True, True],
            "nonempty_material_cap_axes": [True, True, False],
            "analytic_volume": float(analytic_volume),
            "volume_formula": "2*pi*r^2*L - 16*r^3/3",
            "native_field": "min(sqrt(y*y+z*z)-r, sqrt(x*x+z*z)-r), restricted to BOUNDS",
            "periodic_semantics": "Cylinders end on the cell faces and remain inside the other axes. Periodic splitting leaves them unchanged. The final box clip changes no material.",
        },
        "checks": {"public_cylinder_membership_points": len(points),
                   "ridge_points": len(ridges), "paired_cap_probes_per_axis": len(cap_points)},
        "native_retained_constraint_count": None if retained is None else len(retained),
        "native_periodic_affine_eligible": bool(
            _periodic_affine(compiled, (True, True, True), np.array(BOUNDS))),
        "versions": {"microgen": importlib.metadata.version("microgen"),
                     "gmsh": gmsh.__version__},
        "cad": {"available": importlib.util.find_spec("OCP") is not None},
        "scope": "One CustomLattice representative. No geometry or performance coverage claim for the other twelve lattice classes, CAD-only solids, or mapped TPMS.",
    }
    try:
        microgen.cad.require_cad()
    except ImportError as error:
        result["cad"].update(status="blocked", error=str(error))
    else:
        # If CAD has been installed before this preflight runs, check the real
        # Microgen solid. Do not replace it with a Gmsh construction.
        from OCP.BRepCheck import BRepCheck_Analyzer
        from OCP.BRepClass3d import BRepClass3d_SolidClassifier
        from OCP.TopAbs import TopAbs_IN, TopAbs_ON
        from OCP.gp import gp_Pnt

        cad = shape.generate_cad()
        volume = cad.volume()
        valid = BRepCheck_Analyzer(cad.wrapped).IsValid()
        np.testing.assert_allclose(volume, analytic_volume, rtol=1e-8, atol=1e-10)
        assert valid
        classifiers = [BRepClass3d_SolidClassifier(solid)
                       for solid in microgen.cad.enumerate_solids(cad)]
        checked = points[:512]
        cad_inside = []
        for point in checked:
            states = []
            for classifier in classifiers:
                classifier.Perform(gp_Pnt(*point.tolist()), 1e-9)
                states.append(classifier.State())
            assert TopAbs_ON not in states
            cad_inside.append(TopAbs_IN in states)
        np.testing.assert_array_equal(cad_inside, field(*checked.T) < 0)
        step = args.run / "crossed-flat-cylinders.step"
        cad.export_step(str(step))
        result["cad"].update(status="verified", valid=bool(valid), volume=volume,
                             membership_points=len(checked), solid_count=len(classifiers),
                             relative_volume_error=abs(volume / analytic_volume - 1),
                             step_file=str(step.resolve()), step_sha256=workflow.digest(step))

    source_files = {Path(__file__).resolve(), Path(inspect.getfile(type(shape))),
                    Path(inspect.getfile(type(shape).__mro__[1])),
                    Path(inspect.getfile(microgen.Cylinder)),
                    Path(inspect.getfile(microgen.cad)),
                    Path(inspect.getfile(microgen.periodic_split_and_translate)),
                    Path(inspect.getfile(meshers.compile_field)),
                    Path(inspect.getfile(_retained_constraints)),
                    Path(inspect.getfile(meshers.generate))}
    result["inputs"] = {str(path): workflow.digest(path) for path in sorted(source_files)}
    result["samples"] = {str(samples.resolve()): workflow.digest(samples)}
    workflow.save(args.run / "result.json", result)
    print("Exact crossed-cylinder fixture checked; CAD " + result["cad"]["status"] +
          "; mesh admission not tested.", flush=True)


if __name__ == "__main__":
    main()
