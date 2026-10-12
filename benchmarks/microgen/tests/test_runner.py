import importlib.metadata
import json
import subprocess
import sys
import threading
import time

import pytest

from microgen_benchmark.report import render
from microgen_benchmark.runner import (
    failure_catalog,
    provenance,
    stop_tree,
    verify_sources,
    worker,
    write_json,
)


def test_html_payload_cannot_close_script():
    html = render({"error": "</script><script>alert(1)</script>"})
    assert "\\u003c/script>" in html
    assert "<script>alert(1)" not in html
    assert "__BENCHMARK_DATA__" not in html


def test_json_rejects_nonfinite_values(tmp_path):
    with pytest.raises(ValueError):
        write_json(tmp_path / "result.json", {"invalid": float("nan")})


@pytest.mark.skipif(sys.platform != "win32", reason="Windows read handles deny rename")
def test_json_write_waits_for_windows_reader(tmp_path):
    path = tmp_path / "results.json"
    write_json(path, {"old": True})
    reader = path.open(encoding="utf-8")
    assert json.load(reader) == {"old": True}
    release = threading.Timer(0.2, reader.close)
    release.start()
    try:
        write_json(path, {"new": True})
    finally:
        release.join()
        reader.close()
    assert json.loads(path.read_text()) == {"new": True}
    assert not path.with_suffix(".json.tmp").exists()


@pytest.mark.skipif(sys.platform != "win32", reason="Windows read handles deny rename")
def test_json_persistent_windows_lock_fails_without_corrupting_old_result(tmp_path):
    path = tmp_path / "results.json"
    write_json(path, {"old": True})
    with path.open(encoding="utf-8"):
        start = time.perf_counter()
        with pytest.raises(PermissionError):
            write_json(path, {"new": True})
        assert 4.9 <= time.perf_counter() - start < 7
        assert json.loads(path.read_text()) == {"old": True}


def test_changed_source_cannot_mix_results():
    with pytest.raises(RuntimeError, match="source changed"):
        verify_sources({"cases.py": "not-the-actual-hash"})


def test_editable_native_binary_is_recorded_without_distribution_binary(monkeypatch):
    original = importlib.metadata.distribution

    class EditableDistribution:
        version = "0.1.0"
        files = []

    monkeypatch.setattr(
        importlib.metadata,
        "distribution",
        lambda name: EditableDistribution() if name == "meshers" else original(name),
    )
    package = provenance({})["packages"]["meshers"]
    assert package["native_sha256"] == {}
    assert package["loaded_native"]["path"].endswith((".pyd", ".so"))
    assert len(package["loaded_native"]["sha256"]) == 64


def test_worker_rejects_different_native_binary_before_generation(tmp_path):
    from microgen_benchmark.cases import CASES

    request = tmp_path / "request.json"
    write_json(
        request,
        {
            "case": CASES[0].record(),
            "workflow": "meshers-fast",
            "cells": 4,
            "config": {},
            "meshers_native_sha256": "different-build",
        },
    )
    worker(request)
    result = json.loads((tmp_path / "result.json").read_text())
    assert result["status"] == "failed"
    assert "native binary differs" in result["error"]
    assert "mesh_s" not in result


def test_failed_worker_keeps_owner_phase_and_attempt_time(tmp_path, monkeypatch):
    from microgen_benchmark import backends
    from microgen_benchmark.cases import CASES

    def rejected(case, workflow, cells, config, progress):
        progress(stage="Meshers native generation", failure_owner="Meshers", setup_s=0.01)
        raise ValueError("retained test rejection")

    monkeypatch.setattr(backends, "generate", rejected)
    request = tmp_path / "request.json"
    write_json(
        request, {"case": CASES[0].record(), "workflow": "meshers-fast", "cells": 4, "config": {}}
    )
    worker(request)
    result = json.loads((tmp_path / "result.json").read_text())
    assert result["status"] == "failed"
    assert result["stage"] == "Meshers native generation"
    assert result["failure_owner"] == "Meshers"
    assert result["failed_stage_elapsed_s"] >= 0
    assert "mesh_s" not in result and "worst_energy_error" not in result


def test_failure_catalog_keeps_exact_reproduction_policy():
    config = {
        "threads": 2,
        "timeout_s": 73,
        "geometry_tolerance": 0.007,
        "volume_tolerance": 0.08,
        "minimum_quality": 0.12,
        "max_tetrahedra": 12345,
        "pair_cap": 67890,
        "mmg_boundary": "fixed",
    }
    row = {
        "case": "spinodoid",
        "workflow": "meshers-fast",
        "cells": 12,
        "repeat": 0,
        "status": "failed",
        "error": "retained",
        "preview": ["display only"],
    }
    catalog = failure_catalog(
        {"provenance": {"config": config}, "manifest_sha256": "frozen", "rows": [row], "cases": []}
    )
    failure = catalog["failures"][0]
    assert "preview" not in failure
    assert "--geometry-tolerance 0.007" in failure["reproduce"]
    assert "--mmg-boundary fixed" in failure["reproduce"]
    assert "--minimum-quality 0.12" in failure["reproduce"]


def test_owned_tree_cleanup(tmp_path):
    import psutil

    pid_file = tmp_path / "child.json"
    code = "import subprocess,sys,time,json; from pathlib import Path; p=subprocess.Popen([sys.executable,'-c','import time;time.sleep(60)']); Path(sys.argv[1]).write_text(json.dumps(p.pid));time.sleep(60)"
    process = subprocess.Popen([sys.executable, "-c", code, str(pid_file)])
    try:
        for _ in range(100):
            if pid_file.exists():
                break
            time.sleep(0.05)
        assert pid_file.exists()
        child = psutil.Process(json.loads(pid_file.read_text()))
        stop_tree(process)
        assert process.poll() is not None
        assert not child.is_running() or child.status() == psutil.STATUS_ZOMBIE
    finally:
        stop_tree(process)
