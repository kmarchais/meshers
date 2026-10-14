"""Replay admitted geometry through the existing full periodic FEA worker."""
import json
from pathlib import Path
import shutil

import numpy as np
from run import digest, write_json


def qualification(directory):
    row = json.loads((directory / "geometry-result.json").read_text())
    completion = json.loads((directory / "completion.json").read_text())
    manifest = json.loads((directory / "manifest.json").read_text())
    if (completion["state"] != "complete" or completion["worker_exit_code"] != 0
            or completion["stopped"] is not None or not row["admitted_geometry"]
            or row["stage"] != "geometry-only" or row["fea_performed"]):
        raise ValueError("Geometry archive is not a completed admitted geometry-only trial")
    if (not completion.get("rss_scope", "").startswith("Sum of owned worker and recursive")
            or completion["peak_rss_mib"] > manifest["config"]["memory_limit_mib"]
            or completion["worker_seconds"] > manifest["config"]["worker_timeout_s"]):
        raise ValueError("Geometry archive lacks a passing worker-tree resource check")
    if (digest(directory / "manifest.json") != completion["manifest_sha256"]
            or digest(directory / "mesh.npz") != row["mesh_sha256"]):
        raise ValueError("Geometry archive hash changed")
    if any(not row[k] for k in ["geometry_eligible", "quality_eligible", "periodic_eligible",
                                "oriented_caps_unchanged", "retained_coordinates_exact"]):
        raise ValueError("Geometry archive failed qualification")
    if (manifest["config"]["geometry_tolerance"] != .002
            or manifest["config"]["volume_tolerance"] != .01
            or manifest["config"]["minimum_quality"] != .1):
        raise ValueError("Geometry archive did not use the common gates")
    return row, manifest


def freeze_geometry(directory, target):
    directory = Path(directory).resolve()
    row, manifest = qualification(directory)
    for group, root in [("source_sha256", directory / "source"), ("input_sha256", Path())]:
        for name, expected in manifest[group].items():
            if digest(root / name) != expected:
                raise ValueError(f"Changed geometry input: {name}")
    seed = json.loads(Path(manifest["seed_result"]).read_text())
    if seed["status"] != "ok" or seed["case"] != row["case"]:
        raise ValueError("Original geometry seed is not an admitted FEA record")
    target.mkdir()
    for name in ["mesh.npz", "geometry-result.json", "completion.json", "manifest.json"]:
        shutil.copy2(directory / name, target / name)
    shutil.copy2(manifest["seed_mesh"], target / "original-mesh.npz")
    write_json(target / "input-provenance.json", {
        "specimen": manifest["specimen"], "case": seed["case"], "cells": seed["cells"],
        "original_seed_run": str(Path(manifest["seed_result"]).parent.parent),
        "geometry_archive": str(directory), "original_mesh_sha256": seed["mesh_sha256"],
        "scope": "New FEA on geometry-only archive; no prior stiffness or nested-energy claim",
    })


def replay_geometry(directory, task, specimen):
    row, _ = qualification(directory)
    provenance = json.loads((directory / "input-provenance.json").read_text())
    if (task["case"] != row["case"] or task["case"] != provenance["case"]
            or task["cells"] != provenance["cells"] or specimen != provenance["specimen"]):
        raise ValueError("Geometry replay differs from the qualified specimen")
    with np.load(directory / "mesh.npz") as mesh:
        return mesh["points"].copy(), mesh["tetrahedra"].copy(), {
            "generation_s": 0., "mesh_s": 0., "setup_s": 0., "prototype_only": True,
            "generation_timing_scope": "Saved coarsened geometry replay; no meshing speed claim",
            "input_mesh_sha256": row["mesh_sha256"],
        }


def check_replayed_geometry(directory, points, tets, surface):
    from run_surface_coarsening import caps
    qualification(directory)
    provenance = json.loads((directory / "input-provenance.json").read_text())
    if digest(directory / "original-mesh.npz") != provenance["original_mesh_sha256"]:
        raise ValueError("Original seed mesh changed")
    with np.load(directory / "mesh.npz") as saved, np.load(directory / "original-mesh.npz") as original:
        for key, actual in [("points", points), ("tetrahedra", tets), ("surface", surface)]:
            np.testing.assert_array_equal(actual, saved[key])
        retained = saved["retained_vertex_ids"]
        np.testing.assert_array_equal(points, original["points"][retained])
        np.testing.assert_array_equal(caps(points, surface, retained),
                                      caps(original["points"], original["surface"], np.arange(len(original["points"]))))
    return {"exact_qualified_mesh_arrays": True, "retained_coordinates_exact": True,
            "oriented_caps_unchanged": True, "retained_nodes": len(points)}
