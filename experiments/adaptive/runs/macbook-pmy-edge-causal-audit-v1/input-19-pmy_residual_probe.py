"""Localize remaining PMY gaps without generating or modifying a mesh.

Worst fixed reverse samples select at most eight separated neighborhoods. The
complete final Tri6 surface is then tested by Bernstein ball exclusion. Original,
refined, and analytically snapped affine backgrounds are compared through actual
local level-set intersections, not field residual alone. Bounds use floating-point
padding, not interval arithmetic. The parent submits this bounded job to qrun.
"""

import argparse
from itertools import combinations, product
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w


NAME = "tpms-pmy-sheet"
GATE = .002
RADIUS = .04
WITNESSES = 8
TIMEOUT = 240


class AffineBackground:
    def __init__(self, points, cells, values):
        import numpy as np

        self.points, self.cells, self.values = points, cells, values
        self.coordinates = points[cells]
        self.low, self.high = self.coordinates.min(1), self.coordinates.max(1)
        self.node_values = values[cells]
        assert np.isfinite(self.coordinates).all() and np.isfinite(values).all()

    def candidates(self, query, radius=0.):
        import numpy as np

        delta = np.maximum(np.maximum(self.low-query, query-self.high), 0.)
        return np.flatnonzero(np.linalg.norm(delta, axis=1) <= radius+1e-12)

    def interpolation(self, query, level):
        import numpy as np
        from pmy_neck_control import barycentric

        candidates = self.candidates(query)
        bary = barycentric(self.points, self.cells[candidates], query)
        inside = bary.min(1) >= -1e-9
        indices, bary = candidates[inside], bary[inside]
        assert len(indices), "Query escaped the complete background"
        estimates = np.einsum("ij,ij->i", bary, self.node_values[indices])
        assert np.ptp(estimates) < 1e-8, "Discontinuous affine field at a shared face"
        minimum, maximum = self.node_values[indices].min(1), self.node_values[indices].max(1)
        strict = (minimum < level-1e-12) & (maximum > level+1e-12)
        contact = (minimum <= level+1e-12) & (maximum >= level-1e-12)
        return dict(value=float(estimates.mean()), containing_tets=indices.tolist(),
                    vertex_minimum=minimum.tolist(), vertex_maximum=maximum.tolist(),
                    strict_wall_bracket=bool(strict.any()), contact_wall_bracket=bool(contact.any()))


def clipped_cap(points, values, levels):
    """Clip one physical cap triangle by its affine band values."""
    polygon = list(zip(points, values))
    for level, sign in ((levels[0], -1), (levels[1], 1)):
        output = []
        for i, (point, value) in enumerate(polygon):
            other, next_value = polygon[(i+1) % len(polygon)]
            inside, next_inside = sign*(value-level) <= 0, sign*(next_value-level) <= 0
            if inside:
                output.append((point, value))
            if inside != next_inside:
                fraction = (level-value)/(next_value-value)
                output.append((point+fraction*(other-point), level))
        polygon = output
        if not polygon:
            break
    return [point for point, _ in polygon]


def level_polygon(points, values, level):
    import numpy as np

    if (values == level).all():
        return [], True
    cuts = [p for p, v in zip(points, values, strict=True) if v == level]
    for i, j in combinations(range(4), 2):
        if (values[i]-level)*(values[j]-level) < 0:
            fraction = (level-values[i])/(values[j]-values[i])
            cuts.append(points[i]+fraction*(points[j]-points[i]))
    if not cuts:
        return [], False
    cuts = np.unique(np.array(cuts), axis=0)
    if len(cuts) < 3:
        return list(cuts), False
    center = cuts.mean(0)
    first = cuts[0]-center
    first /= np.linalg.norm(first)
    normal = np.linalg.solve(points[1:]-points[:1], values[1:]-values[0])
    second = np.cross(normal, first)
    second /= np.linalg.norm(second)
    angles = np.arctan2((cuts-center) @ second, (cuts-center) @ first)
    return list(cuts[np.argsort(angles)]), False


def affine_boundary_distance(background, query, levels):
    """Cover every possible affine boundary point within the search ball.

    All tets with an AABB meeting the ball are included, plus their material cap
    polygons. A nearest primitive inside the search radius is globally nearest
    for this affine representation. Constant-level tets make exclusion ambiguous.
    """
    import numpy as np
    import scorecard

    triangles, segments, isolated, plateaus = [], [], [], []

    def add(polygon):
        # Exact level contacts can repeat a cap vertex during clipping.
        polygon = list({tuple(point): point for point in polygon}.values())
        if len(polygon) == 1:
            isolated.append(polygon[0])
        elif len(polygon) == 2:
            segments.append(polygon)
        elif len(polygon) >= 3:
            for i in range(1, len(polygon)-1):
                tri = [polygon[0], polygon[i], polygon[i+1]]
                if np.linalg.norm(np.cross(tri[1]-tri[0], tri[2]-tri[0])) > 0:
                    triangles.append(tri)
                else:
                    segments.extend(([tri[0], tri[1]], [tri[1], tri[2]]))

    candidates = background.candidates(query, RADIUS)
    for index in candidates:
        points, values = background.coordinates[index], background.node_values[index]
        for level in levels:
            if values.min() <= level <= values.max():
                polygon, plateau = level_polygon(points, values, level)
                add(polygon)
                if plateau:
                    plateaus.append(int(index))
        for face in combinations(range(4), 3):
            xyz = points[list(face)]
            if any((xyz[:, axis] == side).all() for axis in range(3) for side in (-.5, .5)):
                add(clipped_cap(xyz, values[list(face)], levels))
    distances = []
    if triangles:
        vertices = np.array(triangles).reshape(-1, 3)
        faces = np.arange(len(vertices)).reshape(-1, 3)
        distances.append(float(scorecard.distance_to_mesh(vertices, faces, query[None])[0]))
    for a, b in segments:
        edge = b-a
        fraction = np.clip(float((query-a) @ edge)/float(edge @ edge), 0., 1.) if (edge @ edge) > 0 else 0.
        distances.append(float(np.linalg.norm(query-(a+fraction*edge))))
    distances.extend(float(np.linalg.norm(query-point)) for point in isolated)
    nearest = min(distances) if distances else None
    in_radius = nearest is not None and nearest < RADIUS-1e-12
    lower = nearest if in_radius else RADIUS
    return dict(search_radius=RADIUS, candidate_tets=len(candidates), surface_triangles=len(triangles),
                lower_dimensional_contacts=len(segments)+len(isolated), constant_level_tets=plateaus,
                nearest_distance=nearest if in_radius else None,
                distance_lower_bound=0. if plateaus else lower,
                nearest_is_global_for_affine_representation=bool(in_radius and not plateaus),
                excludes_gate_ball=bool(not plateaus and lower > GATE+1e-12),
                scope="Piecewise-affine band and physical cap polygons; floating-point local enumeration, not the physical implicit boundary")


def controls():
    import numpy as np
    import pmy_failure_probe as probe

    result = probe.controls()
    points = np.array([[-.5, -.5, -.5], [.5, -.5, -.5], [-.5, .5, -.5], [-.5, -.5, .5]])
    cells = np.array([[0, 1, 2, 3]])
    query = np.array([-.01, -.3, -.3])
    background = AffineBackground(points, cells, points[:, 0])
    interpolation = background.interpolation(query, 0.)
    np.testing.assert_allclose(interpolation["value"], -.01, atol=1e-14)
    assert interpolation["strict_wall_bracket"]
    observation = affine_boundary_distance(background, query, (-2., 0.))
    np.testing.assert_allclose(observation["nearest_distance"], .01, atol=1e-12)
    assert observation["excludes_gate_ball"] and observation["nearest_is_global_for_affine_representation"]
    near = affine_boundary_distance(background, np.array([-.001, -.3, -.3]), (-2., 0.))
    assert not near["excludes_gate_ball"]
    empty = affine_boundary_distance(background, query, (2., 3.))
    assert empty["nearest_distance"] is None and empty["excludes_gate_ball"]
    plateau = AffineBackground(points, cells, np.zeros(4))
    assert not affine_boundary_distance(plateau, query, (-2., 0.))["excludes_gate_ball"]
    return dict(result, affine_plane_distance=True, affine_shared_interpolation=True,
                near_affine_wall_retained=True, empty_affine_band_excluded=True,
                constant_level_tet_is_inconclusive=True)


def classification(value, levels):
    return -1 if value < levels[0] else 1 if value > levels[1] else 0


def neighborhood(case, query, level, backgrounds):
    import numpy as np
    from pmy_neck_control import raw_gradient

    normal = raw_gradient(query[None])[0]
    normal /= np.linalg.norm(normal)
    offsets = np.array(list(product((-2*GATE, 0., 2*GATE), repeat=3)))
    offsets = np.vstack((offsets, np.array([-4, -2, -1, -.5, .5, 1, 2, 4])[:, None]*GATE*normal))
    queries = query+offsets
    queries = queries[(abs(queries) <= .5).all(1)]
    exact = np.asarray(case.function(*queries.T))
    reliable = np.min(abs(exact[:, None]-np.array(case.levels)), axis=1) > 1e-8
    exact_classes = np.array([classification(v, case.levels) for v in exact])
    results = {}
    for name, background in backgrounds.items():
        center = background.interpolation(query, level)
        values = [background.interpolation(p, level)["value"] for p in queries]
        classes = np.array([classification(v, case.levels) for v in values])
        mismatch = reliable & (classes != exact_classes)
        results[name] = dict(center=center, boundary=affine_boundary_distance(background, query, case.levels),
                            neighborhood_mismatches=int(mismatch.sum()),
                            false_material=int((mismatch & (classes == 0)).sum()),
                            missing_material=int((mismatch & (exact_classes == 0)).sum()),
                            affine_values=values, affine_classes=classes.tolist(), mismatch_indices=np.flatnonzero(mismatch).tolist())
    return dict(queries=queries.tolist(), exact_values=exact.tolist(), exact_classes=exact_classes.tolist(),
                reliable_classification_samples=int(reliable.sum()),
                ignored_near_level_samples=int((~reliable).sum()), backgrounds=results)


def diagnosis(exclusion, observations):
    if not exclusion["excluded"]:
        return "complete-surface gap not confirmed"
    refined, snapped = (observations["backgrounds"][name] for name in ("refined", "analytic-first-snap"))
    if refined["boundary"]["excludes_gate_ball"]:
        if not refined["center"]["contact_wall_bracket"]:
            return "refined sampling still lacks the wall in the containing tetrahedra"
        return "refined affine geometry is already outside tolerance; root-clipping contribution remains unresolved"
    if refined["boundary"]["nearest_distance"] is not None and refined["boundary"]["nearest_distance"] <= GATE:
        if snapped["boundary"]["excludes_gate_ball"]:
            return "loss appears in the analytic first-snap representation"
        if snapped["boundary"]["nearest_distance"] is not None and snapped["boundary"]["nearest_distance"] <= GATE:
            return "gap appears downstream of the analytic first-snap representation"
    return "stage attribution inconclusive"


def worker(a):
    import numpy as np
    import vtk
    from microgen.shape import surface_functions
    import pmy_failure_probe as probe
    import pmy_neck_control as neck
    import tpms_family_fixtures as fixtures

    vtk.vtkSMPTools.SetBackend("Sequential")
    assert vtk.vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    source_result = json.loads((a.source/"result.json").read_text())
    source_run = a.source.parent
    source_manifest = json.loads((source_run/"manifest.json").read_text())
    prepared = json.loads((source_run/"prepare-result.json").read_text())
    guard = Path(source_manifest["local_volume_guard_evidence"])
    guard_result = json.loads((guard/"result.json").read_text())
    assert source_result["complete"] and source_result["status"] == "admission-rejected"
    assert source_result["geometry"]["failed_gates"] == ["geometry"]
    assert w.digest(a.source/"mesh.npz") == source_result["mesh_sha256"]
    assert w.digest(a.source/"input.json") == prepared["generated_input_hashes"]["one-round"]
    assert w.digest(source_run/"baseline/input.json") == prepared["generated_input_hashes"]["baseline"]
    assert w.digest(guard/"one-round-determinants.npz") == guard_result["arms"]["one-round"]["arrays_sha256"]
    frozen = json.loads((a.fixtures/"result.json").read_text())
    known = {Path(path).resolve(): digest for path, digest in frozen["inputs"].items()}
    specimen = list((a.fixtures/NAME).iterdir())
    assert all(w.digest(path) == known[path.resolve()] for path in specimen)
    files = [Path(__file__), Path(probe.__file__), Path(neck.__file__), Path(fixtures.__file__),
             Path(surface_functions.__file__), a.source/"mesh.npz", a.source/"result.json", a.source/"input.json",
             source_run/"baseline/input.json", source_run/"manifest.json", source_run/"prepare-result.json",
             guard/"result.json", guard/"one-round-determinants.npz", a.fixtures/"result.json", *specimen]
    files += [w.HERE/name for name in ("microgen_workflow.py", "scorecard.py", "mmg_quality_localization.py")]
    hashes = {str(path.resolve()): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix == ".py":
            shutil.copyfile(path, a.run/f"input-{index}-{path.name}")
    w.save(a.run/"manifest.json", dict(inputs=hashes, timeout_s=TIMEOUT, memory_bytes=w.MEMORY,
        numerical_threads=w.NUMERICAL_THREAD_POLICY, witness_budget=WITNESSES, gate_distance=GATE,
        affine_search_radius=RADIUS, selection="Descending fixed-sample distance, with .025 minimum periodic separation",
        native_first_snap_scope="Archived analytic replay, not actual post-rollback native state",
        no_generation=True, no_mesh_edits=True, unchanged_admission=True))
    case = fixtures.case(NAME, a.fixtures)
    with np.load(a.source/"mesh.npz") as mesh:
        points, surface = mesh["points"], mesh["surface"]
    assert surface.shape[1] == 6
    with np.load(a.fixtures/NAME/(NAME+"-samples.npz")) as samples:
        queries = np.concatenate((samples["walls"], samples["caps"]))
        wall_count = len(samples["walls"])
    result = dict(inputs=hashes, source=str(a.source), source_geometry=source_result["geometry"],
                  candidate_generated=False, admission_claim=False, gates_unchanged=True)
    nodes = points[surface]
    with Job(parts={"controls": 1, "reverse-selection": 2, "complete-surface-bounds": 3, "affine-neighborhoods": 4}) as job:
        with job.part("controls") as part:
            result["controls"] = controls()
            part.update(1, "Analytic affine-distance and curved-exclusion controls passed", controls=len(result["controls"]))
        job.result("controls", passed=True)
        with job.part("reverse-selection") as part:
            tree = probe.locator(nodes, divisions=4)
            distances, face_ids = np.empty(len(queries)), np.empty(len(queries), dtype=np.int64)
            closest = [0., 0., 0.]
            cell, sub, squared = vtk.reference(0), vtk.reference(0), vtk.reference(0.)
            for i, query in enumerate(queries):
                tree.FindClosestPoint(query, closest, cell, sub, squared)
                distances[i], face_ids[i] = float(squared)**.5, int(cell)//16
                if i % 512 == 0:
                    part.update((i+1)/len(queries), "Locating fixed reverse samples", completed=i+1, samples=len(queries), maximum_distance=float(distances[:i+1].max()))
            selected = []
            for index in np.argsort(-distances, kind="stable"):
                if distances[index] <= GATE:
                    break
                delta = abs(queries[selected]-queries[index])
                delta = np.minimum(delta, 1-delta)
                if not len(selected) or np.linalg.norm(delta, axis=1).min() >= .025:
                    selected.append(int(index))
                if len(selected) == WITNESSES:
                    break
            assert selected
            result["selection"] = dict(samples=len(queries), locator_subdivisions=4,
                approximate_maximum=float(distances.max()), above_gate=int((distances > GATE).sum()), selected_indices=selected)
            np.savez_compressed(a.run/"reverse-selection.npz", queries=queries, distances=distances, nearest_faces=face_ids, selected=selected)
            part.update(1, "Selected separated worst neighborhoods", witnesses=len(selected), above_gate=result["selection"]["above_gate"])
        job.result("reverse-selection", witnesses=len(selected), above_gate=result["selection"]["above_gate"])
        rows = []
        with job.part("complete-surface-bounds") as part:
            all_controls = probe.bezier(nodes)
            for i, index in enumerate(selected):
                query, face = queries[index], int(face_ids[index])
                exclusion = probe.exclude_ball(all_controls, query, radius=GATE)
                nearest = probe.nearest_patch(nodes[face], query)
                level = min(case.levels, key=lambda x: abs(float(case.function(*query))-x))
                residual = abs(float(case.function(*query))-level)/max(float(np.linalg.norm(neck.raw_gradient(query[None]))), 1e-30)
                assert index >= wall_count or residual < 1e-9
                rows.append(dict(query_index=index, kind="wall" if index < wall_count else "cap", probe=query.tolist(),
                    wall_level=level, field_residual_distance=residual, locator_distance=float(distances[index]),
                    locator_face=face, exact_patch_feasible_upper_bound=nearest,
                    complete_surface_ball_exclusion=exclusion))
                part.update((i+1)/len(selected), "Bounding complete curved surface distance", completed=i+1, witnesses=len(selected), excluded=int(exclusion["excluded"]))
        result["witnesses"] = rows
        w.save(a.run/"partial.json", result)
        with job.part("affine-neighborhoods") as part:
            backgrounds = {}
            for name, directory in (("original", source_run/"baseline"), ("refined", a.source)):
                data = json.loads((directory/"input.json").read_text())
                p, t = np.array(data["points"]), np.array(data["tetrahedra"], dtype=np.int64)
                assert tuple(data["levels"]) == case.levels
                backgrounds[name] = AffineBackground(p, t, np.asarray(case.function(*p.T)))
            with np.load(guard/"one-round-determinants.npz") as replay:
                np.testing.assert_array_equal(replay["original_points"], backgrounds["refined"].points)
                np.testing.assert_array_equal(replay["tetrahedra"], backgrounds["refined"].cells)
                snapped = replay["snapped_points"]
            changed = np.any(snapped != backgrounds["refined"].points, axis=1)
            values = np.asarray(case.function(*np.where(snapped == .5, -.5, snapped).T))
            raw = np.asarray(case.function(*np.where(backgrounds["refined"].points == .5, -.5, backgrounds["refined"].points).T))
            values[changed] = np.where(raw[changed] >= 0, case.levels[1], case.levels[0])
            backgrounds["analytic-first-snap"] = AffineBackground(snapped, backgrounds["refined"].cells, values)
            for i, row in enumerate(rows):
                observations = neighborhood(case, np.array(row["probe"]), row["wall_level"], backgrounds)
                row["neighborhood"] = observations
                row["diagnosis"] = diagnosis(row["complete_surface_ball_exclusion"], observations)
                part.update((i+1)/len(rows), "Compared affine wall intersections and local field classes", completed=i+1, witnesses=len(rows),
                    refined_mismatches=observations["backgrounds"]["refined"]["neighborhood_mismatches"])
                w.save(a.run/"partial.json", result)
        summary = {}
        for row in rows:
            summary[row["diagnosis"]] = summary.get(row["diagnosis"], 0)+1
        result["summary"] = summary
        for name, count in summary.items():
            job.result("diagnosis", classification=name, witnesses=count)
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    result.update(complete=True, input_hashes_verified=True,
        reverse_selection_sha256=w.digest(a.run/"reverse-selection.npz"),
        interpretation="Bernstein exclusion covers every final Tri6 patch with floating-point padding. Affine comparisons explain sampled representations, not exact physical topology. The first-snap archive is an analytic replay, so attribution to a native stage remains conditional.")
    w.save(a.run/"result.json", result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-pmy-native-background-v4/one-round")
    p.add_argument("--fixtures", type=Path, default=w.HERE/"runs/macbook-tpms-family-fixtures-v1")
    p.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = p.parse_args()
    for name in ("run", "source", "fixtures"):
        setattr(a, name, getattr(a, name).resolve())
    if a.worker:
        return worker(a)
    a.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, a.run/"input-pmy_residual_probe.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        a.run/"worker.log", SimpleNamespace(update=lambda *args, **kwargs: None), timeout=TIMEOUT)
    w.save(a.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import result
    findings = json.loads((a.run/"result.json").read_text())
    result("PMY residual witnesses", selected=len(findings["witnesses"]),
           complete_surface_exclusions=sum(row["complete_surface_ball_exclusion"]["excluded"] for row in findings["witnesses"]))
    for name, count in findings["summary"].items():
        result("PMY residual diagnosis", classification=name, witnesses=count)


if __name__ == "__main__":
    main()
