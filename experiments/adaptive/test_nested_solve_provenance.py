"""Preflight format compatibility must retain exact parent evidence."""

import json

import pytest

import microgen_workflow as workflow
from nested_solve import verify_parent_inputs


def test_parent_provenance_supports_both_formats_and_rejects_missing_changed_or_conflicting(tmp_path):
    parent,prepared = tmp_path/'parent',tmp_path/'prepared'
    parent.mkdir();prepared.mkdir()
    inputs = {}
    for name in ('result.json','mesh.npz','fields.npz'):
        path = parent/name
        path.write_bytes(b'fixed parent input')
        inputs[str(path)] = workflow.digest(path)
    manifest,analysis = prepared/'manifest.json',prepared/'analysis-inputs.json'
    manifest.write_text(json.dumps(dict(inputs=inputs)))
    verify_parent_inputs(prepared,parent)
    manifest.write_text(json.dumps(dict(inputs={})))
    with pytest.raises(AssertionError,match='Missing parent provenance'):
        verify_parent_inputs(prepared,parent)
    analysis.write_text(json.dumps(inputs))
    verify_parent_inputs(prepared,parent)
    (parent/'fields.npz').write_bytes(b'different field')
    with pytest.raises(AssertionError,match='Changed parent input'):
        verify_parent_inputs(prepared,parent)
    (parent/'fields.npz').write_bytes(b'fixed parent input')
    manifest.write_text(json.dumps(dict(inputs={str(parent/'fields.npz'):'conflicting'})))
    with pytest.raises(AssertionError,match='Conflicting parent provenance'):
        verify_parent_inputs(prepared,parent)
