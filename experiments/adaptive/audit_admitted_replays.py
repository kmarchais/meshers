"""Transfer exact corner-incidence checks to every archived Schwarz timing arm."""

import argparse
import json
import os
from pathlib import Path
import resource
import sys

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

import numpy as np

import microgen_workflow as workflow
from facet_incidence_audit import audit
from job_parts import worker_part

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

RUNS = workflow.HERE / "runs"
LINEAGE = RUNS / "g3-facet-incidence-lineage-controls-v1"
PHASES = ["warmup", *[f"repeat-{number}" for number in range(1, 6)]]
GROUPS = {
    "linear-microgen": (RUNS / "g3-schwarz-matched-v1/repeat-1-microgen/mesh.npz",
        [RUNS / f"g3-schwarz-matched-v1/{phase}-microgen/mesh.npz" for phase in PHASES], "schwarz-microgen-linear.json"),
    "linear-study-native": (RUNS / "g3-schwarz-matched-v1/repeat-1-meshers/mesh.npz",
        [RUNS / f"g3-schwarz-matched-v1/{phase}-meshers/mesh.npz" for phase in PHASES], None),
    "strong-microgen": (RUNS / "g3-schwarz-microgen-hmax01-seed-v1/mesh.npz",
        [RUNS / f"g3-schwarz-tet10-live-amg-repeated-v1/{phase}/microgen-parent/mesh.npz" for phase in PHASES], "schwarz-microgen-tet10-parent.json"),
    "strong-native": (RUNS / "g3-schwarz-tet10-live-amg-repeated-v1/repeat-1/meshers-tet10/mesh.npz",
        [RUNS / f"g3-schwarz-tet10-live-amg-repeated-v1/{phase}/meshers-tet10/mesh.npz" for phase in PHASES], "schwarz-native.json"),
}


def arrays(path):
    with np.load(path) as archive:
        return tuple(archive[name] for name in ("points", "tetrahedra", "surface"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    resource.setrlimit(resource.RLIMIT_AS, (workflow.MEMORY, workflow.MEMORY))
    inputs, records = {}, {}
    with Job(parts=list(GROUPS)) as job:
        for name, (canonical, paths, audit_name) in GROUPS.items():
            with worker_part(job, name) as part:
                reference = arrays(canonical)
                if audit_name:
                    audit_path = LINEAGE / audit_name
                    checked = json.loads(audit_path.read_text())
                    assert checked["source_sha256"] == workflow.digest(canonical)
                    inputs[str(audit_path)] = workflow.digest(audit_path)
                else:
                    checked = audit(*reference)
                assert checked["oriented_facet_consistency_passed"] and checked["nonpositive_corner_determinants"] == 0
                rows = []
                for index, path in enumerate(paths):
                    current = arrays(path)
                    identical = all(np.array_equal(actual, expected) for actual, expected in zip(current, reference, strict=True))
                    measured = checked if identical else audit(*current)
                    assert measured["oriented_facet_consistency_passed"] and measured["nonpositive_corner_determinants"] == 0
                    inputs[str(path)] = workflow.digest(path)
                    record = {"source": str(path), "exact_arrays_match_audited_control": identical,
                        "oriented_corner_incidence_passed": True}
                    if name == "strong-microgen":
                        lifted = path.parent.parent / "microgen-tet10/mesh.npz"
                        child_points, child_cells, child_surface = arrays(lifted)
                        np.testing.assert_array_equal(child_points[:len(current[0])], current[0])
                        np.testing.assert_array_equal(child_cells[:, :4], current[1])
                        np.testing.assert_array_equal(child_surface[:, :3], current[2])
                        inputs[str(lifted)] = workflow.digest(lifted)
                        record.update(manual_lift_preserves_parent_corners=True, lifted_mesh=str(lifted))
                    rows.append(record)
                    part.update((index + 1) / len(paths), "Checking immutable replay arrays", meshes_checked=index + 1,
                        control_conflicts=checked["same_orientation_shared_faces"])
                records[name] = {"audited_control": str(canonical), "control": checked, "replays": rows}
                job.result(name, replays=len(rows), passed=1)
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())
    workflow.save(args.output, {"verifier_sha256": workflow.digest(Path(__file__)), "inputs": inputs,
        "groups": records, "passed": True, "archived_arm_meshes": 24, "manual_lifts_checked": 6,
        "scope": "Exact replay or direct audit of necessary oriented corner-face consistency; does not certify curved or nonlocal intersection absence; no frozen qualification flags modified"})


if __name__ == "__main__":
    main()
