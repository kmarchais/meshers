"""Geometry-only admission cannot stand in for an audited FEA solution."""
import json

import numpy as np
import pytest

from geometry_replay import qualification, replay_geometry, check_replayed_geometry
from microgen_benchmark.geometry import boundary
from run import digest, write_json
from test_body_refine import cube


@pytest.fixture
def archive(tmp_path):
    points, tets = cube()
    surface = boundary(points, tets)
    np.savez(tmp_path / "mesh.npz", points=points, tetrahedra=tets, surface=surface,
             retained_vertex_ids=np.arange(len(points)))
    np.savez(tmp_path / "original-mesh.npz", points=points, tetrahedra=tets, surface=surface)
    write_json(tmp_path / "manifest.json", {"config": {
        "geometry_tolerance": .002, "volume_tolerance": .01, "minimum_quality": .1,
        "memory_limit_mib": 16384, "worker_timeout_s": 900}})
    write_json(tmp_path / "geometry-result.json", {
        "case": "solid-cell", "stage": "geometry-only", "fea_performed": False,
        "admitted_geometry": True, "mesh_sha256": digest(tmp_path / "mesh.npz"),
        **{k: True for k in ["geometry_eligible", "quality_eligible", "periodic_eligible",
                             "oriented_caps_unchanged", "retained_coordinates_exact"]}})
    write_json(tmp_path / "completion.json", {
        "state": "complete", "worker_exit_code": 0, "stopped": None,
        "manifest_sha256": digest(tmp_path / "manifest.json"), "peak_rss_mib": 100,
        "worker_seconds": 10, "rss_scope": "Sum of owned worker and recursive children"})
    write_json(tmp_path / "input-provenance.json", {
        "case": "solid-cell", "cells": 4, "specimen": {"function": "solid", "offset": 1},
        "original_mesh_sha256": digest(tmp_path / "original-mesh.npz")})
    return tmp_path


def change(path, **values):
    record = json.loads(path.read_text())
    write_json(path, record | values)


@pytest.mark.parametrize("state", ["failed", "rejected"])
def test_failed_geometry_is_rejected_even_with_passing_row(archive, state):
    change(archive / "completion.json", state=state)
    with pytest.raises(ValueError, match="not a completed"):
        qualification(archive)


def test_launcher_only_memory_check_is_rejected(archive):
    change(archive / "completion.json", rss_scope="Launcher RSS only")
    with pytest.raises(ValueError, match="worker-tree"):
        qualification(archive)


def test_changed_mesh_is_rejected(archive):
    with (archive / "mesh.npz").open("ab") as stream:
        stream.write(b"changed")
    with pytest.raises(ValueError, match="hash changed"):
        qualification(archive)


def test_replay_checks_specimen_and_preserves_cap_orientation(archive):
    specimen = {"function": "solid", "offset": 1}
    task = {"case": "solid-cell", "cells": 4}
    with pytest.raises(ValueError, match="specimen"):
        replay_geometry(archive, task, specimen | {"offset": 2})
    points, tets, metrics = replay_geometry(archive, task, specimen)
    assert metrics["generation_s"] == 0
    surface = boundary(points, tets)
    assert check_replayed_geometry(archive, points, tets, surface)["oriented_caps_unchanged"]
    with pytest.raises(AssertionError):
        check_replayed_geometry(archive, points, tets, surface[:, ::-1])
