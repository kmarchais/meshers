"""Inventory the pinned actual Microgen API; untested families stay explicit."""

import ast
import importlib.metadata
from pathlib import Path

import microgen_workflow as w


def catalog():
    distribution = importlib.metadata.distribution("microgen")
    package = Path(distribution.locate_file("microgen"))
    shapes = package / "shape"
    surface = shapes / "surface_functions.py"
    functions = [
        n.name
        for n in ast.parse(surface.read_text()).body
        if isinstance(n, ast.FunctionDef) and not n.name.startswith("_")
    ]
    assert len(functions) == 16
    rows = []
    for name in functions:
        for part in ("sheet", "lower skeletal", "upper skeletal"):
            rows.append(
                dict(
                    id=f"tpms/{name}/{part}",
                    family="periodic fields",
                    api=f"microgen.surface_functions.{name}",
                    part=part,
                    native_status="not tested in expanded fixed protocol",
                    actual_route="Tpms.generate_volume_mesh then public periodic MMG workflow",
                    actual_status="not tested in expanded fixed protocol",
                    periodic=True,
                )
            )
    for name in (
        "Box",
        "Capsule",
        "Cylinder",
        "Ellipsoid",
        "Sphere",
        "ExtrudedPolygon",
        "Polyhedron",
    ):
        cad_only = name in ("ExtrudedPolygon", "Polyhedron")
        rows.append(
            dict(
                id=f"primitive/{name}",
                family="basic solids",
                api=f"microgen.{name}",
                native_status="not tested in expanded fixed protocol",
                actual_route="CAD/Gmsh"
                if cad_only
                else "public implicit volume mesh then MMG; retain separate CAD route",
                actual_status="not tested in expanded fixed protocol",
                periodic=False,
            )
        )
    lattice_init = ast.parse((shapes / "strut_lattice/__init__.py").read_text())
    lattices = sorted(
        {
            a.name
            for n in lattice_init.body
            if isinstance(n, ast.ImportFrom)
            for a in n.names
            if a.name not in ("AbstractLattice",)
        }
    )
    for name in lattices:
        rows.append(
            dict(
                id=f"strut/{name}",
                family="strut lattices",
                api=f"microgen.{name}",
                native_status="exact CAD/implicit equivalence fixture required",
                actual_route="mesh_for_fem: CAD/STEP/Gmsh; public lattice has no implicit field",
                actual_status="not tested",
                periodic=True,
            )
        )
    for name, periodic in (
        ("CylindricalTpms", "rotational or axial where defined"),
        ("SphericalTpms", False),
        ("Infill", False),
        ("Sweep", False),
        ("GradedInfill", False),
    ):
        rows.append(
            dict(
                id=f"mapped/{name}",
                family="mapped and clipped TPMS",
                api=("microgen.shape.tpms." if name in ("Sweep", "GradedInfill") else "microgen.")
                + name,
                native_status="matching physical-domain fixture required",
                actual_status="not tested",
                actual_route="class-specific public volume generation then appropriate MMG route",
                periodic=periodic,
            )
        )
    for name in (
        "Spinodoid isotropic",
        "Spinodoid anisotropic",
        "linear thickness grading",
        "linear density grading",
        "periodic grading",
        "implicit ordered smooth-union lattice",
        "implicit booleans and custom fields",
    ):
        rows.append(
            dict(
                id=f"field/{name}",
                family="spectral, graded and composed fields",
                api="microgen.Shape / Spinodoid / Tpms / implicit_ops",
                native_status="existing fixed specimen tested; parameter coverage incomplete"
                if name != "implicit booleans and custom fields"
                else "representative operator fixtures required",
                actual_status="see unchanged original comparison; expanded parameters untested",
                actual_route="public implicit volume generation then MMG",
                periodic="per frozen specimen",
            )
        )
    paths = [
        package / "__init__.py",
        shapes / "__init__.py",
        surface,
        shapes / "shape.py",
        shapes / "tpms.py",
        shapes / "strut_lattice/__init__.py",
        shapes / "strut_lattice/abstract_lattice.py",
    ]
    return dict(
        version=distribution.version,
        date="2026-11-10",
        rows=rows,
        surface_function_count=len(functions),
        strut_class_count=len(lattices),
        inputs={str(p): w.digest(p) for p in paths},
        official_api="https://3mah.github.io/microgen-docs/microgen.shape.html",
        scope="Finite representatives of each supported family and TPMS part, not every possible parameter or custom field. Surface-only zero-thickness outputs need a separate surface/shell comparison, not tetrahedral volume admission.",
        protocol=dict(
            fixed_original_study_unchanged=True,
            maximum_nodes=250000,
            maximum_tetrahedra=250000,
            required_boundary_distance=0.002,
            required_absolute_relative_volume=0.01,
            required_corner_quality=0.1,
            periodic_rule="All material-cap nodes, complete existing constraint classes and Tri3/Tri6 facets must match",
            curve_rule="Positive sampled curved Jacobians and coherent quadratic associations",
            comparisons="Same host, numerical threads1, admitted actual routes, same order and reported node/DOF costs",
        ),
    )


if __name__ == "__main__":
    w.save(w.ROOT / "handoff/macbook/microgen-coverage-20261110.json", catalog())
