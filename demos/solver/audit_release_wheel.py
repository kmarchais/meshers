# Copyright (c) 2026 Kevin Marchais

"""Generate, export and audit a frozen Meshers Tet4 case without running the solver."""

import argparse
import datetime as dt
import hashlib
import json
from importlib.metadata import version
from pathlib import Path

import pyvista as pv
from periodic_geometry import BASELINE, audit
from release_wheel_comparison import generate, hardware_snapshot


def main() -> int:
    """Return the admission exit code after saving the exported volume and its audit.

    Raises:
        RuntimeError: The installed wheel version differs from the frozen baseline.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("gyroid_uniform", "gyroid_adaptive", "schwarz_native"))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    baseline = json.loads(BASELINE.read_text(encoding="utf8"))
    if version("meshers") != baseline["release"]:
        msg = "Expected the pinned release wheel"
        raise RuntimeError(msg)
    args.output.mkdir()
    row = {
        "recorded_at": dt.datetime.now(dt.UTC).isoformat(),
        "case": args.case,
        "hardware": hardware_snapshot(),
        "versions": {pkg: version(pkg) for pkg in ("meshers", "vtk", "pyvista", "numpy", "scipy")},
    }
    volume = args.output / "volume.vtkhdf"
    try:
        mesh = generate(args.case, baseline["routes"][args.case]["config"], baseline["limits"])
        mesh.write_vtkhdf(volume)
        row["mesh_sha256"] = hashlib.sha256(volume.read_bytes()).hexdigest()
        grid = pv.read(volume)
        row.update(audit(grid, args.case.split("_")[0], baseline))
    except Exception as exc:
        row.update(sampled_admission=False, error=str(exc))
        raise
    finally:
        (args.output / "audit.json").write_text(
            json.dumps(row, indent=2, allow_nan=False) + "\n",
            encoding="utf8",
        )
    return 0 if row["sampled_admission"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
