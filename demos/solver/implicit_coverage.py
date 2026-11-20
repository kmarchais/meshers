"""Replay historical implicit solids with the installed Meshers build and frozen geometry audit."""

import argparse
import datetime as dt
import hashlib
import json
import logging
import sys
import time
from importlib import import_module
from importlib.metadata import version
from pathlib import Path

import meshers
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from benchmark_hardware import hardware_snapshot

ARCHIVE = "602dc900cd028d205b36e29dcdb5b4a9cf4f581a"
ORACLE_DIGEST = "084219b2bb0b03f9c35a408aea3f233b38fbc48e01d771068c69ab15ae23f734"
# Solid fraction and empirical Sobol uncertainty from the archive's
# native-difficult-fields-36 and gyroid-clean-controls calibration records.
VOLUMES = {
    "gyroid": (0.16107845306396484, 0.00027370452880859375),
    "schwarz-p": (0.17131328582763672, 0.00016117095947265625),
    "neovius": (0.22656726837158203, 0.00019550323486328125),
    "graded-gyroid": (0.16110610961914062, 0.000141143798828125),
    "thin-gyroid": (0.08042430877685547, 0.0004425048828125),
    "implicit-lattice": (0.24582576751708984, 0.00031375885009765625),
    "spinodoid": (0.502223014831543, 0.00043487548828125),
}
LIMITS = {"geometry_tolerance": 0.01, "minimum_quality": 0.1, "volume_tolerance": 0.1}
TET_BUDGET = 250_000
PAIR_BUDGET = 2_000_000
ENERGY_TARGET = 0.02


_LOG = logging.getLogger(__name__)


def load_reference(folder):
    """Load the archived case definitions and checks after verifying their source files.

    Returns:
        The benchmark package with its case, geometry and elasticity modules loaded.

    Raises:
        ValueError: The source files differ from the published benchmark.
    """
    package = folder.resolve() / "microgen_benchmark"
    digest = hashlib.sha256()
    names = (
        "__init__",
        "cases",
        "elasticity",
        "geometry",
        "linear_solve",
        "preview",
        "qualification",
        "reference",
    )
    for name in names:
        digest.update((package / f"{name}.py").read_bytes().replace(b"\r\n", b"\n") + b"\0")
    if digest.hexdigest() != ORACLE_DIGEST:
        msg = f"Benchmark source mismatch. Use the Microgen benchmark at {ARCHIVE}"
        raise ValueError(msg)
    sys.path.insert(0, str(package.parent))
    for name in names[1:]:
        module = import_module(f"microgen_benchmark.{name}")
        if Path(module.__file__).resolve().parent != package:
            msg = "A different benchmark package is already loaded"
            raise ValueError(msg)
    return import_module("microgen_benchmark")


def mechanics(solve, mesh, surface, scale):
    """Return the archived quadratic-load energy score and the 2% decision."""
    try:
        metrics, _ = solve(mesh.points, mesh.tetrahedra, surface, scale)
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
        return {"status": "solver-failed", "reason": str(exc)}
    keys = (
        "worst_energy_error",
        "patch_max_error",
        "free_residual",
        "free_dofs",
        "assembly_s",
        "solve_s",
    )
    if not all(np.isfinite(metrics[key]) for key in keys):
        return {"status": "solver-failed", "reason": "Nonfinite solver measurement"}
    return {
        "status": "qualified"
        if metrics["worst_energy_error"] <= ENERGY_TARGET
        else "outside-target",
        **{key: metrics[key] for key in keys},
    }


def measure(case, request, oracle, *, fea):
    """Generate one mesh, then check its returned connectivity against the independent reference.

    Returns:
        Measured counts, admission diagnostics and outcome, including rejected attempts.
    """
    start = time.perf_counter()
    try:
        mesh = meshers.generate(
            case.field,
            bounds=(-0.5, 0.5) * 3,
            band=case.band,
            threads=1,
            minimum_quality=LIMITS["minimum_quality"],
            geometry_tolerance=LIMITS["geometry_tolerance"],
            max_tetrahedra=TET_BUDGET,
            **request,
        )
    except meshers.MeshingError as exc:
        return {
            "status": "native-rejected",
            "reason": str(exc),
            "generation_seconds": time.perf_counter() - start,
        }
    row = {
        "generation_seconds": time.perf_counter() - start,
        "tetrahedra": len(mesh.tetrahedra),
        "nodes": len(mesh.points),
    }
    if not (0 < len(mesh.tetrahedra) <= TET_BUDGET and np.isfinite(mesh.points).all()):
        return {**row, "status": "geometry-rejected", "reason": "Invalid cell count or coordinates"}
    probes = oracle.reference.boundary_samples(case)
    try:
        metrics, surface = oracle.geometry.validate(
            mesh.points,
            mesh.tetrahedra,
            case,
            PAIR_BUDGET,
            probes,
        )
    except ValueError as exc:
        return {**row, "status": "geometry-rejected", "reason": str(exc)}
    except RuntimeError as exc:
        return {**row, "status": "validation-incomplete", "reason": str(exc)}
    if any(isinstance(v, float) and not np.isfinite(v) for v in metrics.values()):
        return {**row, "status": "validation-incomplete", "reason": "Nonfinite audit measurement"}
    volume, uncertainty = VOLUMES[case.id]
    metrics["relative_volume_error"] = abs(metrics["solid_fraction"] / volume - 1)
    admission = oracle.qualification.admission(metrics, LIMITS, uncertainty / volume)
    accepted = admission["geometry_eligible"] and admission["quality_eligible"]
    if fea and accepted:
        row["fea"] = mechanics(oracle.elasticity.solve, mesh, surface, case.scale)
    return {
        **row,
        **metrics,
        "status": "admitted" if accepted else "geometry-rejected",
        "admission_reasons": admission["admission_reasons"],
        "quality_eligible": admission["quality_eligible"],
        "relative_volume_error_with_uncertainty": admission[
            "relative_volume_error_with_uncertainty"
        ],
    }


def main():
    """Record one attempt.

    Returns:
        A failing exit code for incomplete validation or a missed required admission.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference",
        type=Path,
        required=True,
        help="Archived benchmarks/microgen directory",
    )
    parser.add_argument("--case", choices=VOLUMES, required=True)
    parser.add_argument("--cells", type=int, default=8)
    parser.add_argument(
        "--passes",
        type=int,
        default=4,
        help="Number of native quality optimization passes",
    )
    parser.add_argument("--adapt", action="store_true")
    parser.add_argument(
        "--adapt-tolerance",
        type=float,
        help="Internal adaptation target in physical units",
    )
    parser.add_argument(
        "--fea",
        action="store_true",
        help="Solve the archived 2% energy target after geometry admission",
    )
    parser.add_argument(
        "--require-admitted",
        action="store_true",
        help="Exit with failure on rejection",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--snap", type=float, help="Maximum lattice snap as a fraction of spacing")
    args = parser.parse_args()
    if args.adapt_tolerance is not None and (
        not args.adapt or not np.isfinite(args.adapt_tolerance) or args.adapt_tolerance <= 0
    ):
        parser.error("--adapt-tolerance requires --adapt and a finite positive target")
    oracle = load_reference(args.reference)
    case = next(c for c in oracle.cases.CASES if c.id == args.case)
    if args.case == "spinodoid" and version("microgen") != "2.0.0b1":
        parser.error("The archived spinodoid requires microgen==2.0.0b1")
    request = {"cells": args.cells, "optimize_passes": args.passes, "adapt": args.adapt}
    if args.snap is not None:
        if not 0 <= args.snap <= 0.2:
            parser.error("--snap must be between 0 and 0.2")
        request["snap"] = args.snap
    if args.adapt_tolerance is not None:
        request["adapt_tolerance"] = args.adapt_tolerance
    row = {
        "case": args.case,
        "request": request,
        "band": case.band,
        "bounds": [-0.5, 0.5] * 3,
        "threads": 1,
        "tet_budget": TET_BUDGET,
        "overlap_pair_budget": PAIR_BUDGET,
        "limits": LIMITS,
        "energy_target": ENERGY_TARGET if args.fea else None,
        "reference_commit": ARCHIVE,
        "utc": dt.datetime.now(dt.UTC).isoformat(),
        "hardware": hardware_snapshot(),
        "packages": {n: version(n) for n in ("meshers", "numpy", "scipy", "pyvista")},
        "native_sha256": hashlib.sha256(
            Path(import_module("meshers._meshers").__file__).read_bytes(),
        ).hexdigest(),
        "status": "incomplete",
    }
    if args.case == "spinodoid":
        row["packages"]["microgen"] = version("microgen")
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(row, stream, indent=2)
    try:
        row.update(measure(case, request, oracle, fea=args.fea))
    except Exception as exc:
        row.update(status="worker-failed", error_type=type(exc).__name__)
        raise
    finally:
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        with temporary.open("x", encoding="utf-8") as stream:
            json.dump(row, stream, indent=2, allow_nan=False)
            stream.write("\n")
        temporary.replace(args.output)
    _LOG.info(
        "%s, %s cells, %s passes, adapt=%s: %s",
        case.id,
        args.cells,
        args.passes,
        args.adapt,
        row["status"],
    )
    fea_status = row.get("fea", {}).get("status")
    qualified = row["status"] == "admitted" and (not args.fea or fea_status == "qualified")
    return int(
        row["status"] == "validation-incomplete"
        or fea_status == "solver-failed"
        or (args.require_admitted and not qualified),
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    sys.exit(main())
