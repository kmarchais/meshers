import json
import subprocess
import sys
import time

import pytest

from microgen_benchmark.report import render
from microgen_benchmark.runner import stop_tree, verify_sources, write_json


def test_html_payload_cannot_close_script():
    html = render({"error": "</script><script>alert(1)</script>"})
    assert "\\u003c/script>" in html
    assert "<script>alert(1)" not in html
    assert "__BENCHMARK_DATA__" not in html


def test_json_rejects_nonfinite_values(tmp_path):
    with pytest.raises(ValueError):
        write_json(tmp_path / "result.json", {"invalid": float("nan")})


def test_changed_source_cannot_mix_results():
    with pytest.raises(RuntimeError, match="source changed"):
        verify_sources({"cases.py": "not-the-actual-hash"})


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
