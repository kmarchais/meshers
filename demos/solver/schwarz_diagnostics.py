# Copyright (c) 2026 Kevin Marchais

"""Locate Schwarz P wall errors and measure sensitivity to native sizing on the pinned wheel."""

import argparse
import datetime as dt
import json
import sys
import tempfile
from importlib.metadata import version
from pathlib import Path

import fedoo as fd
import numpy as np
import pyvista as pv
from fedoo_homogenization import AuditedCG, ReorderedCG
from historical_schwarz_repair import repaired_grid
from microgen.external import MmgError
from periodic_geometry import BASELINE, audit
from release_wheel_comparison import generate, hardware_snapshot
from schwarz_workflow import PACKAGES, REFERENCE, SchwarzSettings, digest, live_microgen, run
from scipy.linalg import eigvalsh
from vtkmodules.vtkCommonCore import vtkSMPTools

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.performance_gate import cell_geometry

SIZES = {
    "medium": (0.06, 0.0007),
    "fine": (0.04, 0.00045),
    "near-target": (0.038, 0.00041),
    "finest": (0.032, 0.0003),
}
CASES = (
    "default",
    "tuned",
    "single-pass-repair",
    "repaired-solve",
    *SIZES,
    "curved-tet10",
    "shifted-baseline",
    "shifted-rcm",
)


def stage_study(case, folder, row, baseline):
    """Audit the seed, intermediate volumes and the returned mesh on identical reference samples."""
    control = case == "single-pass-repair"
    walls = {}
    stages = {}
    row["stages"] = stages

    def observe(name, mixed):
        grid = mixed.extract_cells(np.flatnonzero(mixed.celltypes == pv.CellType.TETRA))
        stages[name] = audit(grid, "schwarz", baseline)
        skin = grid.extract_surface(algorithm="dataset_surface")
        tri = skin.regular_faces
        coords = skin.points[tri]
        caps = np.zeros(len(tri), dtype=bool)
        for side in (-0.5, 0.5):
            caps |= (abs(coords - side) < 1e-10).all(axis=1).any(axis=1)
        walls[name] = pv.PolyData.from_regular_faces(skin.points, tri[~caps])

    path = live_microgen(
        SchwarzSettings(
            "single-pass-repair" if control else f"microgen-{case}",
            30 if control else 24,
        ),
        folder,
        row,
        baseline["limits"]["max_nodes_or_elements"],
        observe=observe,
    )
    grid = pv.read(path)
    observe("export", grid)
    if control:
        observe("repaired", repaired_grid(grid, row))
    last = next(reversed(stages))
    loc = stages[last]["boundary"]["wall 0"]["reference_to_mesh"]["location"]
    row["tracked_reference_xyz"] = loc
    row["distance_at_tracked_point"] = {}
    for stage, wall in walls.items():
        _, nearest = wall.find_closest_cell(loc, return_closest_point=True)
        row["distance_at_tracked_point"][stage] = float(np.linalg.norm(nearest - loc))
    row["status"] = "audited"


def quadratic_control(row, baseline):
    """Solve a coarse curved mesh as a stiffness control.

    Raises:
        ValueError: The volume, count, symmetry or load checks fail.
    """
    cfg = dict(baseline["routes"]["schwarz_native"]["config"], cells=12, adapt=False, order=2)
    row["config"] = cfg
    mesh = generate("schwarz_native", cfg, baseline["limits"])
    pts, cells = mesh.points, mesh.tetrahedra
    row.update(nodes=len(pts), elements=len(cells))
    if max(len(pts), len(cells)) > baseline["limits"]["max_nodes_or_elements"]:
        msg = "Quadratic mesh exceeds the node or element cap"
        raise ValueError(msg)
    volumes = [cell_geometry(pts[batch]).sum() for batch in np.array_split(cells, 32) if len(batch)]
    row["integrated_volume"] = float(sum(volumes))
    if abs(row["integrated_volume"] / 0.5 - 1) > baseline["limits"]["relative_volume_error"]:
        msg = "Curved volume outside the Schwarz limit"
        raise ValueError(msg)
    fd.ModelingSpace("3D")
    fem = fd.Mesh(pts + 0.5, cells.copy(), "tet10")
    material = fd.constitutivelaw.ElasticIsotrop(2.5, 0.25)
    assembly = fd.Assembly.create(fd.weakform.StressEquilibrium(material), fem)
    solver = AuditedCG()
    tensor = np.asarray(
        fd.homogen.get_homogenized_stiffness(
            assembly,
            meshperio=True,
            solver=solver,
            rigid_body_constraint="pin",
        ),
    )
    if len(solver.records) != 6 or not np.allclose(tensor, tensor.T, rtol=1e-8, atol=1e-9):
        msg = "Quadratic six-load stiffness verification failed"
        raise ValueError(msg)
    reference = np.asarray(json.loads(REFERENCE.read_text())["stiffness"])
    reference = (reference + reference.T) / 2
    row.update(
        stiffness=tensor.tolist(),
        loads=solver.records,
        relative_energy_discrepancy=float(
            abs(eigvalsh((tensor + tensor.T) / 2 - reference, reference)).max(),
        ),
        status="quadratic-control",
    )


def measure(case, parent, row, *, repair_resolution, repair_hmax):
    """Run one diagnostic under the pinned geometry and solver settings.

    Raises:
        RuntimeError: The installed packages or VTK backend differ from the study.
    """
    vtkSMPTools.SetBackend("Sequential")
    if row["versions"] != PACKAGES or vtkSMPTools.GetBackend() != "Sequential":
        msg = "Install the pinned packages with sequential VTK support"
        raise RuntimeError(msg)
    baseline = json.loads(BASELINE.read_text())
    with tempfile.TemporaryDirectory(dir=parent) as temp:
        if case in SIZES or case.startswith("shifted-"):
            size, tolerance = SIZES.get(case, SIZES["finest"])
            run(
                SchwarzSettings(
                    "native",
                    24,
                    native_sizing={"max_length": size, "adapt_tolerance": tolerance},
                    shifted=case.startswith("shifted-"),
                ),
                Path(temp),
                row,
                solver_factory=ReorderedCG if case == "shifted-rcm" else AuditedCG,
            )
        elif case == "repaired-solve":
            run(
                SchwarzSettings("microgen-repaired", repair_resolution, repair_hmax=repair_hmax),
                Path(temp),
                row,
            )
        elif case == "curved-tet10":
            quadratic_control(row, baseline)
        else:
            stage_study(case, Path(temp), row, baseline)


def main():
    """Write one independent result, preserving failed and over-budget cases.

    Raises:
        MmgError: The recorded MMG command failed.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=CASES)
    parser.add_argument("output", type=Path)
    parser.add_argument("--repair-resolution", type=int, choices=(30, 34, 36))
    parser.add_argument(
        "--repair-hmax",
        type=float,
        help="Maximum edge length in the first MMG pass",
    )
    args = parser.parse_args()
    if args.repair_hmax is not None and (
        args.case != "repaired-solve" or not 0 < args.repair_hmax <= 1
    ):
        parser.error("--repair-hmax requires repaired-solve and a length in (0, 1]")
    if args.repair_resolution is not None and args.case != "repaired-solve":
        parser.error("--repair-resolution applies to repaired-solve")
    row = {
        "case": args.case,
        "recorded_at": dt.datetime.now(dt.UTC).isoformat(),
        "hardware": hardware_snapshot(),
        "versions": {n: version(n) for n in PACKAGES},
        "baseline_sha256": digest(BASELINE),
        "reference_sha256": digest(REFERENCE),
    }
    with args.output.open("x", encoding="utf8") as out:
        try:
            measure(
                args.case,
                args.output.parent,
                row,
                repair_resolution=args.repair_resolution or 30,
                repair_hmax=args.repair_hmax or 0.1,
            )
        except MmgError:
            row.update(status="failed", error="MmgError: MMG failed with the recorded arguments")
            raise
        except Exception as exc:
            row.update(status="failed", error=f"{type(exc).__name__}: {exc}")
            raise
        finally:
            json.dump(row, out, indent=2, allow_nan=False)
            out.write("\n")


if __name__ == "__main__":
    main()
