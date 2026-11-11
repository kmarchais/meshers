"""Localize the saved spinodoid volume-correction distance regression.

No mesh moves and no fixtures are regenerated. Each metric is evaluated at the
same Cartesian probes within a saved mesh; corresponding barycentric locations
compare the three saved meshes. This is a sampling/metric diagnosis, not a
new admission result or a nearest-distance certificate.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w
import numpy as np
import comparison_fixtures as fixtures
from box_projection import constrained_foot
import scorecard
import volume_fidelity_control as volume


NAME = "spinodoid"
RUNS = w.HERE/"runs"
SOURCES = {
    "baseline": RUNS/"macbook-native-face-fit-nine-v1/spinodoid-repeat-0-fit",
    "global": RUNS/"macbook-volume-bias-seven-v1/spinodoid",
    "facewise": RUNS/"macbook-volume-facewise-pilot-v1/spinodoid",
}


def point_metrics(target, probes, *, project=False):
    raw = np.asarray(target.function(*probes.T))
    gradient = scorecard.gradient(target.function, probes)
    norm = np.linalg.norm(gradient, axis=1)
    assert np.all(norm > 1e-10)
    locked = np.isclose(abs(probes), .5, rtol=0, atol=1e-12)
    allowed = gradient.copy()
    allowed[locked] = 0
    allowed_norm = np.linalg.norm(allowed, axis=1)
    stationary = allowed_norm < 1e-14
    assert np.all(abs(raw[stationary]) <= 1e-12), "Stationary off-wall constrained probe"
    cap_residual = np.divide(raw, allowed_norm, out=np.zeros_like(raw), where=~stationary)
    result = dict(physical=raw/norm, cap_constrained=cap_residual,
                  locked=locked, gradient_norm=norm, allowed_gradient_norm=allowed_norm,
                  raw=raw)
    if project:
        foot = probes.copy()
        for _ in range(40):
            g = scorecard.gradient(target.function, foot)
            g[locked] = 0
            residual = target.function(*foot.T)
            step = residual[:, None]*g/np.maximum(np.sum(g*g, axis=1), 1e-30)[:, None]
            foot -= step
            if np.max(np.linalg.norm(step, axis=1), initial=0.) < 1e-12:
                break
        remaining = abs(target.function(*foot.T))/np.maximum(
            np.linalg.norm(scorecard.gradient(target.function, foot), axis=1), 1e-30)
        failed = (remaining > 1e-9) | np.any(abs(foot) > .5+1e-9, axis=1)
        fallback = []
        for i in np.flatnonzero(failed):
            candidate, info = constrained_foot(target.function, probes[i], 0., locked[i], foot[i])
            fallback.append(dict(index=int(i), **info))
            if info["feasible_target"]:
                foot[i] = candidate
        remaining = abs(target.function(*foot.T))/np.maximum(
            np.linalg.norm(scorecard.gradient(target.function, foot), axis=1), 1e-30)
        failed = (remaining > 1e-9) | np.any(abs(foot) > .5+1e-9, axis=1)
        result.update(projected=np.linalg.norm(foot-probes, axis=1), foot=foot,
                      failed=failed, fallback=fallback)
    return result


def controls():
    plane = scorecard.Case("slanted-plane", lambda x, y, z: x+2*y-.1, (-np.inf, 0.))
    # Equal field residuals. The first point is locked to x=.5, so its
    # permissible nearest foot is farther away than the unconstrained foot.
    probes = np.array([[.5, -.199, .1], [.2, -.049, .1]])
    values = point_metrics(plane, probes, project=True)
    physical = .002/np.sqrt(5)
    np.testing.assert_allclose(values["physical"], physical, atol=1e-12, rtol=0)
    np.testing.assert_allclose(values["cap_constrained"], [.001, physical], atol=1e-12, rtol=0)
    np.testing.assert_allclose(values["projected"], [.001, physical], atol=1e-12, rtol=0)
    assert np.array_equal(values["locked"], [[True, False, False], [False, False, False]])
    assert not values["failed"].any()
    return dict(slanted_plane_physical_residual=values["physical"].tolist(),
                slanted_plane_cap_residual=values["cap_constrained"].tolist(),
                slanted_plane_projected_distance=values["projected"].tolist(),
                cap_vs_interior_distinction_passed=True)


def scan(target, points, wall, surface_ids, bary, *, project, update):
    basis = volume.tri_basis(bary)[0]
    fields = ("physical", "cap_constrained", "projected") if project else ("physical", "cap_constrained")
    maxima = {k: np.empty(len(wall)) for k in fields}
    argmax = {k: np.empty(len(wall), dtype=int) for k in fields}
    categories = {k: dict(cap=None, interior=None) for k in fields}
    failures, fallbacks = 0, 0
    for start in range(0, len(wall), 192):
        end = min(start+192, len(wall))
        probes = np.einsum("sn,fnd->fsd", basis, points[wall[start:end]]).reshape(-1, 3)
        metrics = point_metrics(target, probes, project=project)
        on_cap = metrics["locked"].any(1)
        for k in fields:
            values = abs(metrics[k]).reshape(end-start, len(bary))
            maxima[k][start:end] = values.max(1)
            argmax[k][start:end] = values.argmax(1)
            for category, mask in (("cap", on_cap), ("interior", ~on_cap)):
                if not mask.any():
                    continue
                selected = np.flatnonzero(mask)
                index = int(selected[np.argmax(abs(metrics[k][selected]))])
                value = float(abs(metrics[k][index]))
                prior = categories[k][category]
                if prior is None or value > prior["maximum"]:
                    face, probe = divmod(index, len(bary))
                    categories[k][category] = dict(maximum=value, wall_index=start+face,
                        surface_face=int(surface_ids[start+face]), probe_index=probe,
                        barycentric=bary[probe].tolist(), point=probes[index].tolist())
        if project:
            failures += int(metrics["failed"].sum())
            fallbacks += len(metrics["fallback"])
        update(end/len(wall), end*len(bary))
    summaries = {}
    for k in fields:
        top = np.argsort(maxima[k])[-4:][::-1]
        summaries[k] = dict(maximum=float(maxima[k].max()),
            categories=categories[k],
            worst=[dict(wall_index=int(i), surface_face=int(surface_ids[i]),
                        value=float(maxima[k][i]), probe_index=int(argmax[k][i]),
                        barycentric=bary[argmax[k][i]].tolist()) for i in top])
    return dict(metrics=summaries, projection_failures=failures,
                projection_fallbacks=fallbacks), maxima, argmax


def local_comparison(target, meshes, wall, surface_ids, grids, scans, output):
    # Select only from the already-computed frozen observations. No new sampling
    # density, parameter sweep, correction or witness-fitted mesh is introduced.
    selected_faces = set()
    locations = set()
    for mode in scans:
        for grid in scans[mode]:
            for metric in scans[mode][grid]["metrics"].values():
                winner = metric["worst"][0]
                selected_faces.add(winner["wall_index"])
                locations.add((winner["wall_index"], grid, winner["probe_index"]))
    selected_faces = sorted(selected_faces)
    fields = ("physical", "cap_constrained", "projected")
    patches, arrays = [], {}
    for face in selected_faces:
        rows = {}
        for grid, bary in grids.items():
            basis = volume.tri_basis(bary)[0]
            state = {}
            for mode, points in meshes.items():
                probes = basis@points[wall[face]]
                values = point_metrics(target, probes, project=True)
                assert not values["failed"].any()
                state[mode] = values
                arrays[f"face{face}-{grid}-{mode}-points"] = probes
                for k in fields:
                    arrays[f"face{face}-{grid}-{mode}-{k}"] = values[k]
            rows[grid] = {}
            for mode, values in state.items():
                rows[grid][mode] = {
                    k: dict(maximum=float(abs(values[k]).max()),
                            argmax=int(abs(values[k]).argmax()),
                            barycentric=bary[abs(values[k]).argmax()].tolist(),
                            signed_value=float(values[k][abs(values[k]).argmax()]))
                    for k in fields}
            rows[grid]["same_barycentric_maximum_change"] = {
                mode: {k: float(np.max(abs(state[mode][k])-abs(state["baseline"][k]))) for k in fields}
                for mode in ("global", "facewise")}
        patches.append(dict(wall_index=face, surface_face=int(surface_ids[face]),
                            node_ids=wall[face].tolist(), grids=rows))
    witnesses = []
    for face, grid, index in sorted(locations):
        bary = grids[grid][index]
        basis = volume.tri_basis(bary[None])[0]
        states = {}
        for mode, points in meshes.items():
            probe = (basis@points[wall[face]])[0]
            values = point_metrics(target, probe[None], project=True)
            locked = values["locked"][0]
            # Independent local minimization at standard/dense winners checks
            # whether an iterative foot-path artifact explains the regression.
            closest = None
            if grid in ("standard", "dense13"):
                foot, info = constrained_foot(target.function, probe, 0., locked, values["foot"][0])
                closest = dict(foot=foot.tolist(), **info)
            states[mode] = dict(point=probe.tolist(), locked=locked.tolist(),
                raw=float(values["raw"][0]), physical_signed=float(values["physical"][0]),
                cap_constrained_signed=float(values["cap_constrained"][0]),
                gradient_norm=float(values["gradient_norm"][0]),
                allowed_gradient_norm=float(values["allowed_gradient_norm"][0]),
                projected=float(values["projected"][0]), foot=values["foot"][0].tolist(),
                independent_local_minimum=closest)
        witnesses.append(dict(wall_index=face, surface_face=int(surface_ids[face]),
                              grid=grid, probe_index=index, barycentric=bary.tolist(), states=states))
    np.savez_compressed(output/"patch-metrics.npz", **arrays)
    return patches, witnesses


def worker(a):
    w.imports()
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    fixtures_root = RUNS/"macbook-comparison-fixtures-v1"
    target = fixtures.case(NAME, fixtures_root)
    assert target.levels == (-np.inf, 0.)
    control = controls()
    paths = [Path(__file__), Path(volume.__file__), Path(fixtures.__file__),
             Path(w.__file__), w.HERE/"scorecard.py", w.HERE/"box_projection.py",
             fixtures.coefficients_path(NAME, fixtures_root)]
    meshes, records = {}, {}
    cells = surface = None
    for mode, source in SOURCES.items():
        record = json.loads((source/"result.json").read_text())
        records[mode] = record
        if mode != "baseline":
            coefficient_path = fixtures.coefficients_path(NAME, fixtures_root).resolve()
            assert w.digest(coefficient_path) == record["inputs"][str(coefficient_path)], \
                "Spectral coefficients differ from the frozen correction experiment"
        assert w.digest(source/"mesh.npz") == record["mesh_sha256"]
        paths += [source/"result.json", source/"mesh.npz", source/(NAME+"-samples.npz"),
                  source/(NAME+"-specimen.json")]
        with np.load(source/"mesh.npz") as mesh:
            meshes[mode] = mesh["points"].copy()
            if cells is None:
                cells, surface = mesh["tetrahedra"].copy(), mesh["surface"].copy()
            else:
                np.testing.assert_array_equal(cells, mesh["tetrahedra"])
                np.testing.assert_array_equal(surface, mesh["surface"])
    for mode in ("global", "facewise"):
        np.testing.assert_array_equal(meshes[mode][np.unique(cells[:, :4])],
                                      meshes["baseline"][np.unique(cells[:, :4])])
        for suffix in ("-samples.npz", "-specimen.json"):
            assert w.digest(SOURCES[mode]/(NAME+suffix)) == w.digest(SOURCES["baseline"]/(NAME+suffix))
    hashes = {str(p.resolve()): w.digest(p) for p in paths}
    wall, _, _ = volume.wall_data(target, meshes["baseline"], surface)
    lookup = {tuple(face): i for i, face in enumerate(surface)}
    surface_ids = np.array([lookup[tuple(face)] for face in wall])
    grids = dict(degree6=volume.barycentric(6), gauss6=volume.quadrature(6)[0],
                 gauss8=volume.quadrature(8)[0], standard=scorecard.BARY.copy(),
                 dense13=np.unique(np.vstack((scorecard.BARY, volume.barycentric(13))), axis=0))
    w.save(a.run/"manifest.json", dict(inputs=hashes, cases=[NAME], modes=list(meshes),
           grids={k: len(v) for k, v in grids.items()}, worker_timeout_s=300,
           memory_bytes=w.MEMORY, numerical_threads=w.NUMERICAL_THREAD_POLICY,
           no_mesh_changes=True, no_regenerated_specimen_samples=True,
           scope="Compare three saved meshes and metrics at identical surface probes; no admission or nearest-distance certificate"))
    for path in paths:
        if path.suffix == ".py":
            shutil.copyfile(path, a.run/("input-"+path.name))
    scans, arrays = {}, {}
    parts = {f"{mode}-{grid}": 3 if grid == "dense13" else 1 for mode in meshes for grid in grids}
    parts["local-patches"] = 1
    with Job(parts=parts) as job:
        for mode, points in meshes.items():
            scans[mode] = {}
            for grid, bary in grids.items():
                with job.part(f"{mode}-{grid}") as part:
                    part.update(0., "scanning frozen wall probes", force=True, faces=len(wall), probes=0)
                    summary, maxima, argmax = scan(target, points, wall, surface_ids, bary,
                        project=grid in ("standard", "dense13"),
                        update=lambda fraction, count: part.update(fraction, "scanning frozen wall probes", probes=count, faces=len(wall)))
                    scans[mode][grid] = summary
                    for k in maxima:
                        arrays[f"{mode}-{grid}-{k}-max"] = maxima[k]
                        arrays[f"{mode}-{grid}-{k}-argmax"] = argmax[k]
                    if grid in ("standard", "dense13"):
                        expected = records["global"]["before"]["geometry"] if mode == "baseline" else records[mode]["after"]["geometry"]
                        if grid == "dense13":
                            expected = records["global"]["dense_before"] if mode == "baseline" else records[mode]["dense_after"]
                        np.testing.assert_allclose(summary["metrics"]["projected"]["maximum"],
                                                   expected["sampled_forward_max"], atol=1e-12, rtol=0)
                        assert summary["projection_failures"] == 0
                    job.result(f"{mode}-{grid}", physical_max=summary["metrics"]["physical"]["maximum"],
                               cap_max=summary["metrics"]["cap_constrained"]["maximum"],
                               projected_max=summary["metrics"].get("projected", {}).get("maximum"))
        with job.part("local-patches") as part:
            part.update(0., "comparing metrics on selected frozen patches", force=True, meshes=3, grids=5)
            patches, witnesses = local_comparison(target, meshes, wall, surface_ids, grids, scans, a.run)
            part.update(1., "patch comparison complete", faces=len(patches), witnesses=len(witnesses))
            job.result("local-patches", faces=len(patches), witnesses=len(witnesses))
    np.savez_compressed(a.run/"face-maxima.npz", **arrays)
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    w.save(a.run/"result.json", dict(inputs=hashes, complete=True, controls=control, scans=scans,
        patches=patches, witnesses=witnesses, saved_forward_maxima_reproduced=True,
        arrays={name: w.digest(a.run/name) for name in ("patch-metrics.npz", "face-maxima.npz")},
        source_meshes_and_fixtures_unchanged=True,
        limits="Projected feet and local minimizers are feasible numerical targets, not global nearest-distance certificates. Finite face grids do not certify continuous maxima."))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--worker", action="store_true")
    a = p.parse_args()
    a.run = a.run.resolve()
    if a.worker:
        return worker(a)
    a.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, a.run/"input-volume_residual_probe.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), "--worker", "--run", str(a.run)],
        a.run/"worker.log", SimpleNamespace(update=lambda *args, **kwargs: None), timeout=300)
    w.save(a.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import result
    saved = json.loads((a.run/"result.json").read_text())
    for mode in SOURCES:
        result(f"{mode} saved forward maxima",
               standard=saved["scans"][mode]["standard"]["metrics"]["projected"]["maximum"],
               dense13=saved["scans"][mode]["dense13"]["metrics"]["projected"]["maximum"],
               jsonl=str(a.run/"reported-results.jsonl"))


if __name__ == "__main__":
    main()
