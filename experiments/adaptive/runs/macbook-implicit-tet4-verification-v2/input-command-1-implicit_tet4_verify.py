"""Replay the Tet4 artifacts and field identities for the fixed implicit study."""
import argparse
import ast
import inspect
import json
from pathlib import Path
import subprocess
import sys

import microgen_workflow as w


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    a = p.parse_args()
    a.run.mkdir(parents=True, exist_ok=False)
    outputs = {}
    for label, pilot, correction, repeated in (
        ("default", "count-pilot", "count-correction", "default-repeated"),
        ("modified", "modified-pilot", "modified-correction", "modified-repeated")):
        output = a.run/(label+"-verification.json")
        subprocess.run([sys.executable, str(w.HERE/"tpms_tet4_count_verify.py"),
            "--selection", str(w.ROOT/f"handoff/macbook/implicit-tet4-{label}-final-20261111.json"),
            "--studies", *[str(w.HERE/f"runs/macbook-implicit-tet4-{suffix}-v1/study")
                for suffix in (pilot, correction, repeated)],
            "--report", str(w.ROOT/f"handoff/macbook/implicit-tet4-{label}-comparison-20261111.json"),
            "--output", str(output)], check=True)
        outputs[label] = json.loads(output.read_text())

    np, _, meshers, microgen, *_ = w.imports()
    import comparison_fixtures as fixtures
    import primitive_family_fixtures as primitives
    selection = json.loads((w.ROOT/"handoff/macbook/implicit-tet4-default-final-20261111.json").read_text())
    controls, resolved_sources = [], []
    package = Path(microgen.__file__).parent
    inputs = {str(path):w.digest(path) for path in (
        package/"shape/shape.py", package/"shape/implicit_ops.py", package/"shape/spinodoid.py",
        package/"shape/surface_functions.py", Path(fixtures.__file__), Path(primitives.__file__))}
    transfer = w.ROOT.parent/"meshers-handoff-20261024/data-manifest.json"
    transferred = json.loads(transfer.read_text())
    relative = "experiments/adaptive/runs/g3-spinodoid-compiled32-v1/fourier-coefficients.json"
    spectral_receipt = next(f for repo in transferred["repositories"] if repo["name"] == "meshers"
        for f in repo["files"] if f["path"] == relative)
    assert w.digest(w.ROOT/relative) == spectral_receipt["sha256"]
    inputs[str(w.ROOT/relative)] = spectral_receipt["sha256"]
    receipt_path = w.ROOT/"handoff/macbook/comparison-fixture-receipt-20261026.json"
    receipt = json.loads(receipt_path.read_text())
    anisotropic = w.HERE/"runs/macbook-comparison-fixtures-v1/spinodoid-anisotropic-coefficients.json"
    assert w.digest(anisotropic) == receipt["anisotropic_coefficients_sha256"]
    inputs[str(anisotropic)] = receipt["anisotropic_coefficients_sha256"]
    inputs[str(receipt_path)] = w.digest(receipt_path)
    points = np.random.default_rng(20261111).uniform(-.5, .5, (4096, 3))
    for entry in selection["cases"]:
        name, root = entry["case"], Path(entry["fixtures"])
        frozen = json.loads((root/"result.json").read_text())
        if name.startswith("primitive-"):
            expected = {Path(path):h for path,h in frozen["inputs"].items()}
            target = primitives.case(name, root)
            shape, _ = primitives.actual_shape(name, root)
            values = shape.evaluate(*points.T)
            protocol = json.loads((root/"protocol.json").read_text())
            expected.update({Path(path):h for path,h in protocol["sources"].items()})
            original = Path(primitives.__file__)
            if w.digest(original) != expected[original]:
                archived = w.HERE/"runs/macbook-primitive-family-prepare-v1/input-command-1-primitive_family_fixtures.py"
                assert w.digest(archived) == expected.pop(original)
                old_tree = ast.parse(archived.read_text())
                tree = ast.parse(original.read_text())
                tree.body = [n for n in tree.body if not isinstance(n, ast.FunctionDef) or n.name != "retained_field"]
                assert ast.dump(tree) == ast.dump(old_tree), "Frozen primitive definitions or preparation changed"
                expected[archived] = w.digest(archived)
                resolved_sources.append(dict(original=str(original), archived=str(archived),
                    historical_sha256=w.digest(archived), current_sha256=w.digest(original),
                    proof="Full module AST identical after removing the later retained_field addition; formatting ignored"))
            if name in ("primitive-box", "primitive-cylinder"):
                retained = primitives.retained_field(name, root)(*points.T)
                np.testing.assert_array_equal(retained <= 0, target.function(*points.T) <= 0)
        else:
            expected = ({root/path:h for path,h in frozen["hashes"].items()} if name.endswith("density-linear") else
                {root/name/path:h for path,h in frozen["cases"][name]["fixture_hashes"].items()})
            target = fixtures.case(name, root)
            shape, _ = fixtures.actual_shape(name, root)
            if name.startswith("spinodoid"):
                values = shape.evaluate(*(points+.5).T)
            elif name.startswith("graded-"):
                values = shape.raw_field(*points.T)/(shape._offset_func(*points.T)/2)
            else:
                assert name == "smooth-lattice"
                values = shape.evaluate(*points.T)
        for path, digest in expected.items():
            assert path.exists() and w.digest(path) == digest, path
            inputs[str(path)] = digest
        native = target.function(*points.T)
        np.testing.assert_allclose(values, native, atol=1e-12, rtol=1e-12)
        source = Path(inspect.getfile(type(shape)))
        inputs[str(source)] = w.digest(source)
        controls.append(dict(case=name, probes=len(points), maximum_error=float(abs(values-native).max()),
            frozen_receipt_verified=True))
    runtime = json.loads((w.ROOT/"handoff/macbook/tpms-tet4-runtime-20261111.json").read_text())
    assert all(w.digest(Path(path)) == digest for path,digest in runtime["binaries"].items())
    combined = w.ROOT/"handoff/macbook/implicit-tet4-comparison-20261111.json"
    report = json.loads(combined.read_text())
    assert report["specimen_count"] == 11 and report["case_count"] == 19
    assert all(w.digest(Path(path)) == digest for path,digest in report["sources"].items())
    for case in report["cases"]:
        expected = not (case["case"].startswith("primitive-") or case["case"] in
            ("graded-gyroid-linear", "graded-gyroid-density-linear"))
        assert all(arm["periodic_required"] == expected for arm in case["arms"].values())
        if case["case"].startswith("primitive-"):
            assert case["reference"]["kind"] == "analytic solid volume"
    w.save(a.run/"verification.json", dict(passed=True, controls=controls, inputs=inputs,
        historical_sources_resolved=resolved_sources,
        original_spectral_transfer_receipt=spectral_receipt, transfer_manifest_sha256=w.digest(transfer),
        groups={key:{k:v for k,v in data.items() if k not in ("artifacts", "distinct_meshes")}
            for key,data in outputs.items()}, runtime_binaries_verified=True,
        combined_report_sha256=w.digest(combined), driver_sha256=w.digest(Path(__file__)),
        scope="Replayed actual arrays and metrics, frozen receipts, public Microgen field identity and unchanged native/MMG binaries. Existing geometry and periodic admission remain unchanged. No FEA certification."))
    print("Verified eleven frozen field identities and both Tet4 evidence groups")


if __name__ == "__main__":
    main()
