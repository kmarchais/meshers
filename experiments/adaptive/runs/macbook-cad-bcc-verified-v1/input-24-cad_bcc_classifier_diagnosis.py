"""Diagnose frozen BCC CAD membership failures without changing any mesh gate."""

import argparse
import inspect
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w
import cad_lattice_fixture as fixture

SOURCE = w.HERE/"runs/macbook-cad-bcc-coverage-v1"


def state_name(state):
    from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_OUT
    return {TopAbs_IN: "inside", TopAbs_ON: "boundary", TopAbs_OUT: "outside"}.get(state, "unknown")


def inventory(cad):
    import numpy as np
    import microgen.cad
    from OCP.BRep import BRep_Tool, BRep_Builder
    from OCP.BRepAdaptor import BRepAdaptor_Surface, BRepAdaptor_Curve
    from OCP.BRepCheck import BRepCheck_Analyzer, BRepCheck_Shell, BRepCheck_NoError
    from OCP.BRepTools import BRepTools
    from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX, TopAbs_SHELL, TopAbs_INTERNAL, TopAbs_EXTERNAL
    from OCP.TopExp import TopExp
    from OCP.TopTools import TopTools_IndexedMapOfShape
    from OCP.TopoDS import TopoDS, TopoDS_Compound

    entities, records = {}, {}
    for name, kind, cast in (("faces", TopAbs_FACE, TopoDS.Face_s), ("edges", TopAbs_EDGE, TopoDS.Edge_s),
                            ("vertices", TopAbs_VERTEX, TopoDS.Vertex_s), ("shells", TopAbs_SHELL, TopoDS.Shell_s)):
        indexed = TopTools_IndexedMapOfShape(); TopExp.MapShapes_s(cad.wrapped, kind, indexed)
        entities[name] = [cast(indexed.FindKey(index)) for index in range(1, indexed.Extent()+1)]
        rows = []
        for index, entity in enumerate(entities[name]):
            row = dict(index=index, orientation=str(entity.Orientation()), closed_flag=bool(entity.Closed()))
            if name != "shells":
                row["tolerance"] = float(BRep_Tool.Tolerance_s(entity))
            else:
                closure = BRepCheck_Shell(entity).Closed()
                row.update(closure_check=str(closure), closure_passed=closure == BRepCheck_NoError)
            if name == "vertices":
                point = BRep_Tool.Pnt_s(entity); row["point"] = [point.X(), point.Y(), point.Z()]
            if name == "faces":
                adaptor = BRepAdaptor_Surface(entity, True)
                row.update(surface_type=str(adaptor.GetType()), uv_bounds=list(BRepTools.UVBounds_s(entity)))
                if "Plane" in row["surface_type"]:
                    plane = adaptor.Plane(); normal = plane.Axis().Direction(); point = plane.Location()
                    n, p = np.array([normal.X(), normal.Y(), normal.Z()]), np.array([point.X(), point.Y(), point.Z()])
                    row.update(plane_normal=n.tolist(), plane_point=p.tolist(), plane_distance_from_origin=float(abs(n@p)))
            if name == "edges":
                row["degenerate"] = bool(BRep_Tool.Degenerated_s(entity))
                if not row["degenerate"]:
                    adaptor = BRepAdaptor_Curve(entity)
                    row.update(curve_type=str(adaptor.GetType()), bounds=[adaptor.FirstParameter(), adaptor.LastParameter()])
            rows.append(row)
        records[name] = rows
    compound = TopoDS_Compound(); builder = BRep_Builder(); builder.MakeCompound(compound)
    boundary_face_ids = []
    for index, face in enumerate(entities["faces"]):
        if face.Orientation() not in (TopAbs_INTERNAL, TopAbs_EXTERNAL):
            builder.Add(compound, face); boundary_face_ids.append(index)
    summaries = {}
    for name in ("faces", "edges", "vertices"):
        values = np.array([row["tolerance"] for row in records[name]])
        summaries[name] = dict(count=len(values), minimum=float(values.min()), median=float(np.median(values)), maximum=float(values.max()))
    records.update(tolerances=summaries, solid_count=len(microgen.cad.enumerate_solids(cad)),
        valid=bool(BRepCheck_Analyzer(cad.wrapped).IsValid()), volume=cad.volume(),
        shells_closed=bool(records["shells"]) and all(row["closure_passed"] for row in records["shells"]),
        boundary_face_ids=boundary_face_ids, central_support_planes=[row for row in records["faces"]
            if row.get("plane_distance_from_origin", 1.) <= 1e-8])
    return SimpleNamespace(cad=cad, entities=entities, boundary=compound, boundary_face_ids=boundary_face_ids), records


def classify(cad, points, part):
    import microgen.cad
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.gp import gp_Pnt

    solids = microgen.cad.enumerate_solids(cad)
    reused = [BRepClass3d_SolidClassifier(solid) for solid in solids]
    rows = []
    combine = lambda states: "inside" if "inside" in states else "boundary" if "boundary" in states else "unknown" if "unknown" in states or not states else "outside"
    for index, point in enumerate(points):
        p = gp_Pnt(*point.tolist())
        for classifier in reused:
            classifier.Perform(p, 1e-9)
        repeated = combine([state_name(classifier.State()) for classifier in reused])
        fresh = combine([state_name(BRepClass3d_SolidClassifier(solid, p, 1e-9).State()) for solid in solids])
        rows.append(dict(reused=repeated, fresh=fresh))
        if index % 128 == 0:
            part.update(index/max(len(points), 1), "Comparing reused and fresh CAD classifiers", probes=index)
    part.update(1, "Frozen classifier queries completed", probes=len(points))
    return rows


def distances(context, query):
    import numpy as np
    from OCP.BRep import BRep_Tool
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.BRepTools import BRepTools
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
    from OCP.Extrema import Extrema_ExtPC
    from OCP.gp import gp_Pnt, gp_Pnt2d

    xyz = lambda point: np.array([point.X(), point.Y(), point.Z()])
    p = gp_Pnt(*query.tolist())
    point_shape = BRepBuilderAPI_MakeVertex(p).Vertex()
    distance = BRepExtrema_DistShapeShape()
    distance.SetMultiThread(False)
    distance.SetDeflection(1e-12)
    distance.LoadS1(point_shape); distance.LoadS2(context.boundary); distance.Perform()
    brep = dict(done=bool(distance.IsDone()))
    if distance.IsDone():
        solutions = []
        for index in range(1, distance.NbSolution()+1):
            foot = xyz(distance.PointOnShape2(index))
            solutions.append(dict(point=foot.tolist(), euclidean_distance=float(np.linalg.norm(query-foot)),
                support_type=str(distance.SupportTypeShape2(index))))
        brep.update(value=float(distance.Value()), inner_solution=bool(distance.InnerSolution()), solutions=solutions)
    surface_rows, failures = [], []
    for index in context.boundary_face_ids:
        face = context.entities["faces"][index]
        try:
            bounds = BRepTools.UVBounds_s(face)
            projection = GeomAPI_ProjectPointOnSurf(p, BRep_Tool.Surface_s(face), *bounds, 1e-12)
            for solution in range(1, projection.NbPoints()+1):
                u, v = projection.Parameters(solution)
                foot = xyz(projection.Point(solution))
                trim = state_name(BRepClass_FaceClassifier(face, gp_Pnt2d(u, v), 1e-9).State())
                surface_rows.append(dict(face=index, distance=float(np.linalg.norm(query-foot)), foot=foot.tolist(),
                    uv=[u, v], trim_state=trim, face_tolerance=float(BRep_Tool.Tolerance_s(face))))
        except Exception as error:
            failures.append(dict(face=index, error=str(error)))
    surface_rows.sort(key=lambda row: row["distance"])
    edge_rows, edge_failures = [], []
    for index, edge in enumerate(context.entities["edges"]):
        if BRep_Tool.Degenerated_s(edge):
            continue
        adaptor = BRepAdaptor_Curve(edge)
        a, b = adaptor.FirstParameter(), adaptor.LastParameter()
        candidates = [(xyz(adaptor.Value(a)), a), (xyz(adaptor.Value(b)), b)]
        try:
            extrema = Extrema_ExtPC(p, adaptor, a, b, 1e-12)
            if not extrema.IsDone():
                edge_failures.append(dict(edge=index, error="Extrema not done"))
            else:
                for solution in range(1, extrema.NbExt()+1):
                    found = extrema.Point(solution)
                    candidates.append((xyz(found.Value()), found.Parameter()))
        except Exception as error:
            edge_failures.append(dict(edge=index, error=str(error)))
        foot, parameter = min(candidates, key=lambda candidate: np.linalg.norm(candidate[0]-query))
        edge_rows.append(dict(edge=index, distance=float(np.linalg.norm(foot-query)), foot=foot.tolist(),
            parameter=float(parameter), edge_tolerance=float(BRep_Tool.Tolerance_s(edge))))
    edge_rows.sort(key=lambda row: row["distance"])
    vertex_rows = [dict(vertex=index, distance=float(np.linalg.norm(xyz(BRep_Tool.Pnt_s(vertex))-query)),
        vertex_tolerance=float(BRep_Tool.Tolerance_s(vertex))) for index, vertex in enumerate(context.entities["vertices"])]
    vertex_rows.sort(key=lambda row: row["distance"])
    outside = np.maximum(abs(query)-.5, 0)
    return dict(point=query.tolist(), analytic_box_signed_max=float(abs(query).max()-.5),
        analytic_outside_box_distance=float(np.linalg.norm(outside)), boundary_brep=brep,
        nearest_support_surface=surface_rows[:1], nearest_trim_in_surface=[row for row in surface_rows if row["trim_state"] == "inside"][:1],
        nearest_trim_on_surface=[row for row in surface_rows if row["trim_state"] == "boundary"][:1],
        nearest_edge=edge_rows[:1], nearest_vertex=vertex_rows[:1], face_projection_failures=failures,
        edge_extrema_failures=edge_failures, support_surface_scope="Underlying bounded surface projections; trim ON can inherit topology tolerance, so it is not an independent membership certificate")


def worker(args):
    np, _, _, microgen, *_ = w.imports()
    import microgen.cad
    import OCP.OCP
    from vtkmodules.vtkCommonCore import vtkSMPTools
    assert vtkSMPTools.SetBackend("Sequential")
    vtkSMPTools.Initialize(1)
    assert vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    prior = json.loads((args.source/"result.json").read_text())
    assert not prior["fixture_verified"] and not prior["native_attempted"]
    assert w.digest(args.source/"actual-bcc.step") == prior["cad"]["step_sha256"]
    cad_inputs = {Path(path): digest for path, digest in prior["inputs"].items()
        if "/microgen/" in path or "/OCP/" in path or "dist-info/" in path
        or path == str(Path(fixture.__file__).resolve())}
    for path, digest in cad_inputs.items():
        assert w.digest(path) == digest, f"Changed CAD reconstruction source: {path}"
    with np.load(args.source/"membership-probes.npz") as archive:
        points, old_states, normalized, raw = (archive[key] for key in ("points", "cad_states", "normalized_field", "raw_field"))
    mismatch_ids = [row["index"] for row in prior["membership"]["mismatches"]]
    center = int(np.flatnonzero((points == 0).all(1))[0])
    selected = sorted(set(mismatch_ids+[center]))
    paths = {Path(__file__).resolve(), Path(fixture.__file__).resolve(), Path(w.__file__).resolve(),
        Path(microgen.cad.__file__), Path(inspect.getfile(microgen.BodyCenteredCubic)),
        Path(inspect.getfile(microgen.BodyCenteredCubic.__mro__[1])), Path(OCP.OCP.__file__)}
    paths.update(args.source/name for name in ("result.json", "actual-bcc.step", "primitive-receipt.json",
        "raw-primitive-receipt.json", "membership-probes.npz", "membership-probes-input.npz", "manifest.json"))
    paths.update(cad_inputs)
    hashes = {str(path): w.digest(path) for path in sorted(paths)}
    for index, path in enumerate(sorted(paths)):
        if path.suffix == ".py":
            shutil.copyfile(path, args.run/f"input-{index}-{path.name}")
    w.save(args.run/"manifest.json", dict(inputs=hashes, worker_timeout_s=240, memory_bytes=w.MEMORY,
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY, mismatch_ids=mismatch_ids, center_id=center,
        vtk_backend=vtkSMPTools.GetBackend(),
        queries_unchanged=True, classifier_tolerance=1e-9, distance_deflection=1e-12,
        scope="Diagnosis only; no mesh, new offsets, tolerance repair, admission relaxation or fixture promotion"))
    result = dict(inputs=hashes, complete=False, geometry_verification_remains_blocked=True, native_attempted=False,
        mismatch_ids=mismatch_ids, center_id=center, center=dict(point=points[center].tolist(),
            saved_cad_state=str(old_states[center]), normalized_field=float(normalized[center]), raw_field=float(raw[center])),
        comparisons={}, controls={})
    with Job(parts={"load": 1, "controls": 1, "step-classification": 2, "step-distances": 3,
                    "fresh-classification": 2, "fresh-distances": 3}) as job:
        with job.part("load") as part:
            part.update(0, "Loading saved CAD and reconstructing the same recipe", frozen_queries=len(points))
            imported = microgen.cad.import_step(args.source/"actual-bcc.step")
            rebuilt = microgen.BodyCenteredCubic(strut_radius=.1, cell_size=1., strut_joints=False).generate_cad()
            contexts = {}
            for name, cad in (("step", imported), ("fresh", rebuilt)):
                contexts[name], data = inventory(cad)
                result["comparisons"][name] = dict(inventory=data)
                w.save(args.run/(name+"-inventory.json"), data)
            part.update(1, "Stored CAD entity tolerances and topology", step_faces=len(contexts["step"].entities["faces"]))
        with job.part("controls") as part:
            part.update(0, "Checking a unit box and a finite cylinder", controls=0)
            control_cases = {
                "box": (microgen.cad.make_box([1., 1., 1.], [0., 0., 0.]),
                    np.array([[0., 0., 0.], [.49999, 0., 0.], [.5, 0., 0.], [.50001, 0., 0.],
                              [.49999, .49999, .49999], [.50001, .50001, .50001]]),
                    ["inside", "inside", "boundary", "outside", "inside", "outside"],
                    [.5, 1e-5, 0., 1e-5, 1e-5, float(np.sqrt(3)*1e-5)]),
                "cylinder": (microgen.cad.make_cylinder(.1, .5, [0., 0., 0.]),
                    np.array([[0., 0., 0.], [.24999, 0., 0.], [.25, 0., 0.], [.25001, 0., 0.],
                              [0., .09999, 0.], [0., .1, 0.], [0., .10001, 0.]]),
                    ["inside", "inside", "boundary", "outside", "inside", "boundary", "outside"],
                    [.1, 1e-5, 0., 1e-5, 1e-5, 0., 1e-5])}
            for name, (cad, queries, expected, expected_distances) in control_cases.items():
                context, data = inventory(cad)
                states = classify(cad, queries, SimpleNamespace(update=lambda *a, **k: None))
                distance_rows = [distances(context, point) for point in queries]
                distance_passed = all(row["boundary_brep"]["done"] and
                    abs(row["boundary_brep"]["value"]-truth) <= 1e-9
                    for row, truth in zip(distance_rows, expected_distances, strict=True))
                result["controls"][name] = dict(inventory=data, expected=expected, classifications=states,
                    expected_boundary_distances=expected_distances, boundary_distances_passed=distance_passed,
                    passed=distance_passed and all(row["reused"] == truth and row["fresh"] == truth for row, truth in zip(states, expected, strict=True)),
                    distances=distance_rows)
            part.update(1, "Known-geometry controls recorded", controls=len(control_cases))
        for name in ("step", "fresh"):
            with job.part(name+"-classification") as part:
                rows = classify(contexts[name].cad, points, part)
                result["comparisons"][name].update(classifications=rows,
                    reuse_disagreements=sum(row["fresh"] != row["reused"] for row in rows),
                    differences_from_saved=sum(row["reused"] != old for row, old in zip(rows, old_states, strict=True)))
                w.save(args.run/(name+"-classifications.json"), rows)
            with job.part(name+"-distances") as part:
                rows = []
                for offset, index in enumerate(selected):
                    row = dict(index=index, saved_state=str(old_states[index]), normalized_field=float(normalized[index]),
                        raw_field=float(raw[index]), **distances(contexts[name], points[index]))
                    rows.append(row)
                    part.update((offset+1)/len(selected), "Measuring frozen mismatch distances to CAD geometry",
                        completed=offset+1, total=len(selected), center=index == center)
                w.save(args.run/(name+"-distances.json"), rows)
                result["comparisons"][name].update(distance_queries=len(rows), center_distance=next(row for row in rows if row["index"] == center),
                    surface_projection_failures=sum(len(row["face_projection_failures"]) for row in rows),
                    curve_extrema_failures=sum(len(row["edge_extrema_failures"]) for row in rows))
            job.result(name, reuse_disagreements=result["comparisons"][name]["reuse_disagreements"],
                differences_from_saved=result["comparisons"][name]["differences_from_saved"],
                maximum_face_tolerance=result["comparisons"][name]["inventory"]["tolerances"]["faces"]["maximum"],
                maximum_vertex_tolerance=result["comparisons"][name]["inventory"]["tolerances"]["vertices"]["maximum"])
            w.save(args.run/"result.json", result)
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    result.update(complete=True, input_hashes_verified=True)
    w.save(args.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args(); args.source, args.run = args.source.resolve(), args.run.resolve()
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import report, result
    report(0, "Starting frozen BCC CAD classifier diagnosis", eta_s=240, memory_ceiling_gib=8)
    receipt = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **k: None), timeout=240)
    w.save(args.run/"worker.json", receipt)
    assert receipt["returncode"] == 0, "Read worker.log and worker.json"
    row = json.loads((args.run/"result.json").read_text())
    result("bcc-classifier-diagnosis", complete=row["complete"], geometry_verification_remains_blocked=True,
        controls_passed=all(value["passed"] for value in row["controls"].values()),
        step_reuse_disagreements=row["comparisons"]["step"]["reuse_disagreements"],
        fresh_reuse_disagreements=row["comparisons"]["fresh"]["reuse_disagreements"])


if __name__ == "__main__":
    main()
