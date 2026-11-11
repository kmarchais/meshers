"""One bounded Tet4 Gyroid count control using the existing comparison worker."""
import argparse
import json
from pathlib import Path
import shutil
import sys

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--adaptive", action="store_true", help="One geometry-driven successor to the rejected coarse control")
    parser.add_argument("--count-correction", action="store_true", help="One count correction from the measured13-cell output")
    args = parser.parse_args()
    args.run = args.run.resolve()
    args.run.mkdir(parents=True, exist_ok=False)
    _, _, meshers, *_ = w.imports()
    from comparison_benchmark import configs
    from job_parts import worker_part
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    findings = w.ROOT/"handoff/macbook/comparison-findings-20261026.json"
    rows = json.loads(findings.read_text())["rows"]
    native = next(r for r in rows if r["case"] == "gyroid" and r["name"] == "meshers-tet4")
    actual = next(r for r in rows if r["case"] == "gyroid" and r["name"] == "actual-default-tet4")
    seed = round(24*(actual["elements"]["median"]/native["elements"]["median"])**(1/3))
    assert seed == 13
    correction_source = None
    if args.count_correction:
        assert not args.adaptive
        correction_source = w.HERE/"runs/macbook-gyroid-tet4-budget-v1/result.json"
        measured = json.loads(correction_source.read_text())
        assert measured["complete"] and measured["candidate"]["elements"] == 15074
        seed = round(seed*(actual["elements"]["median"]/measured["candidate"]["elements"])**(1/3))
        assert seed == 12
    config = dict(next(c for c in configs("gyroid") if c["name"] == "meshers-tet4"),
                  cells=seed, periodic_trace=True, declared_unit_cell=True)
    if args.adaptive:
        config.update(name="meshers-tet4-coarse-adaptive", adapt=True,
                      tolerance=.002, adapt_tolerance=.002, min_length=.002)
    fixture_root = w.HERE/"runs/macbook-comparison-fixtures-v1"
    candidate = args.run/"candidate"
    candidate.mkdir()
    paths = [Path(__file__).resolve(), findings, Path(meshers._meshers.__file__)]
    if correction_source is not None:
        paths.append(correction_source)
    paths += [w.HERE/n for n in ("comparison_benchmark.py", "comparison_fixtures.py",
        "comparison_mesh_audit.py", "microgen_workflow.py", "scorecard.py", "periodic_trace.py",
        "box_projection.py", "facet_incidence_audit.py", "quadratic_edge_removal_probe.py", "job_parts.py")]
    for suffix in ("-samples.npz", "-specimen.json"):
        source = fixture_root/"gyroid"/("gyroid"+suffix)
        shutil.copyfile(source, candidate/source.name)
        paths += [source, candidate/source.name]
    hashes = {str(p): w.digest(p) for p in paths}
    shutil.copyfile(__file__, args.run/"driver.py")
    w.save(args.run/"manifest.json", dict(inputs=hashes, config=config,
        target_tets=actual["elements"]["median"], target_nodes=actual["nodes"]["median"],
        selection="Round24*(12250/72757)^(1/3); one coarse seed, optional geometry-driven successor after coarse rejection",
        count_correction_source=str(correction_source) if correction_source else None,
        count_correction_rule="Round13*(12250/15074)^(1/3)=12" if correction_source else None,
        adaptation_scope="Adaptive mode includes native geometry tolerance0.002, which also constrains surface motion; not an isolated switch ablation",
        numerical_threads=w.NUMERICAL_THREAD_POLICY, memory_bytes=w.MEMORY, timeout_s=180,
        actual_comparator_admitted=False, gates_changed=False,
        scope="Tet4 only. Count target does not qualify the rejected actual mesh. No FEA or speed-ratio claim."))
    with Job(parts=["candidate"]) as job:
        with worker_part(job, "candidate") as part:
            record = w.bounded([sys.executable, str(w.HERE/"comparison_benchmark.py"),
                "--run", str(candidate), "--fixtures", str(fixture_root), "--case", "gyroid",
                "--worker", json.dumps(config)], candidate/"worker.log", part, timeout=180)
        w.save(args.run/"worker.json", record)
        assert record["returncode"] == 0, "Read candidate/worker.log"
        row = json.loads((candidate/"result.json").read_text())
        assert all(w.digest(Path(p)) == h for p, h in hashes.items())
        w.save(args.run/"result.json", dict(complete=True, input_hashes_verified=True, candidate=row,
            target_tets=actual["elements"]["median"], target_nodes=actual["nodes"]["median"],
            admitted=row["qualified"], resources=record))
        job.result("Gyroid Tet4 count control", tets=row.get("elements"), nodes=row.get("nodes"),
                   admitted=row["qualified"], returncode=record["returncode"])


if __name__ == "__main__":
    main()
