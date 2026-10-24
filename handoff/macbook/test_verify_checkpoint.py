"""Source history lookup must find exact bytes and reject mismatched blobs."""

import hashlib
from types import SimpleNamespace

import verify_checkpoint as verifier


def test_measured_python_is_found_after_head_changes(monkeypatch):
    verifier.committed_python_ref.cache_clear()
    expected = hashlib.sha256(b"measured source").hexdigest()
    def show(command, *, cwd, capture_output):
        assert command[:2] == ["git", "show"]
        return SimpleNamespace(returncode=0, stdout=b"measured source" if command[2].startswith("measured:") else b"different")
    monkeypatch.setattr(verifier.subprocess, "run", show)
    monkeypatch.setattr(verifier.subprocess, "check_output", lambda command, **_: b"measured\n")
    assert verifier.committed_python_ref("experiments/adaptive/example.py", expected, "legacy") == "measured"
    assert verifier.committed_python_ref("experiments/adaptive/example.py", "0" * 64, "legacy") is None
