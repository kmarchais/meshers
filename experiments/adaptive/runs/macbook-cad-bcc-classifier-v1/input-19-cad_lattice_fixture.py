"""Experimental finite-cylinder field adapter for Microgen CAD fixtures.

Original primitives remain in the receipt. Periodic images are bounded by
their analytic support boxes. Only exactly collinear anchored cylinder chains
are coalesced, removing internal cap zeros without rounding a positive gap.
The union identity applies to the endpoint-anchored representation. Its axes
can differ from source CAD rotations at recorded arithmetic roundoff, so
independent CAD checks remain required. There is no all-class coverage claim.
"""

from itertools import product

import microgen_workflow as w
import numpy as np


ARITHMETIC = 64*np.finfo(float).eps


def primitives(shape):
    assert not shape.strut_joints, "Spherical joint adaptation is outside this control"
    records = []
    for index, (center, rotation, height, pair) in enumerate(zip(shape.strut_centers,
            shape.strut_rotations, shape.strut_heights, shape.strut_vertex_pairs, strict=True)):
        endpoints = np.asarray(shape.vertices)[pair]
        direction = endpoints[1]-endpoints[0]
        length = float(np.linalg.norm(direction))
        axis = rotation.as_matrix()[:, 0]
        orientation_error = float(min(np.linalg.norm(axis-direction/length), np.linalg.norm(axis+direction/length)))
        assert orientation_error <= ARITHMETIC
        assert float(height) == length, "Cannot certify source endpoint anchors for nonstandard strut heights"
        assert np.array_equal(center, endpoints.mean(axis=0))
        nominal = np.stack((center-.5*height*axis, center+.5*height*axis))
        endpoint_error = min(float(abs(nominal-endpoints).max()), float(abs(nominal-endpoints[::-1]).max()))
        assert endpoint_error <= ARITHMETIC*max(1., length)
        records.append(dict(id=index, center=np.asarray(center).tolist(), axis=axis.tolist(),
            rotation=rotation.as_matrix().tolist(), height=float(height), radius=float(shape.strut_radius),
            endpoints=endpoints.tolist(), endpoint_roundoff_bound=ARITHMETIC*max(1., length),
            endpoint_reconstruction_error=endpoint_error, rotation_direction_error=orientation_error))
    return records


def periodic_images(records, bounds):
    lo, hi = np.asarray(bounds)
    periods = hi-lo
    images, proof = [], []
    for original in records:
        center, axis = np.asarray(original["center"]), np.asarray(original["axis"])
        extent = .5*original["height"]*abs(axis)+original["radius"]*np.sqrt(np.maximum(1-axis*axis, 0))
        error = ARITHMETIC*max(1., float(abs(center).max()), original["height"])
        lower, upper = center-extent-error, center+extent+error
        # The installed Microgen route deliberately skips an axis when a
        # primitive crosses both opposite planes. Do not replace that semantics.
        assert not np.any((lower < lo-error*2) & (upper > hi+error*2)), "Two-sided primitive wrapping is unsupported"
        first = np.ceil((lo-upper)/periods).astype(int)
        last = np.floor((hi-lower)/periods).astype(int)
        excluded_below = lo-(upper+(first-1)*periods)
        excluded_above = lower+(last+1)*periods-hi
        assert (excluded_below > 0).all() and (excluded_above > 0).all()
        for shift in product(*(range(a, b+1) for a, b in zip(first, last, strict=True))):
            translation = np.asarray(shift)*periods
            images.append(dict(original, source_id=original["id"], image=list(shift),
                center=(center+translation).tolist(),
                endpoints=(np.asarray(original["endpoints"])+translation).tolist()))
        proof.append(dict(source_id=original["id"], support_lower=lower.tolist(), support_upper=upper.tolist(),
            arithmetic_padding=error, first_shift=first.tolist(), last_shift=last.tolist(),
            nearest_excluded_lower_margin=excluded_below.tolist(), nearest_excluded_upper_margin=excluded_above.tolist()))
    return images, dict(bounds=[lo.tolist(), hi.tolist()], periods=periods.tolist(), primitives=proof,
        argument="Any omitted translated support box is strictly separated from the cell in at least one axis",
        limitation="Image-range bound, not a formal CAD/field equivalence proof")


def coalesce(records):
    groups = []
    for record in records:
        endpoints = np.asarray(record["endpoints"])
        direction = endpoints[1]-endpoints[0]
        for group in groups:
            # No angular tolerance or proximity threshold can join two lines.
            if record["radius"] == group[0]["radius"] and np.array_equal(np.cross(direction, group[1]), np.zeros(3)) and np.array_equal(
                    np.cross(endpoints[0]-np.asarray(group[0]["endpoints"])[0], group[1]), np.zeros(3)):
                group[2].append(record)
                break
        else:
            groups.append((record, direction, [record]))
    merged, reports, former_caps = [], [], []
    for _, direction, group in groups:
        axis_index = int(abs(direction).argmax())
        intervals = []
        for record in group:
            ends = np.asarray(record["endpoints"])
            ends = ends[np.argsort(ends[:, axis_index])]
            intervals.append((float(ends[0, axis_index]), float(ends[1, axis_index]), ends, record))
        intervals.sort(key=lambda value: (value[0], value[1]))
        clusters = []
        for item in intervals:
            if not clusters:
                clusters.append([item])
                continue
            last_end = max(previous[1] for previous in clusters[-1])
            gap = item[0]-last_end
            error = ARITHMETIC*max(1., abs(last_end), abs(item[0]))
            if 0 < gap <= error:
                raise ValueError(f"Ambiguous positive axial gap {gap}; no merging across it")
            if gap > 0:
                clusters.append([item])
            else:
                if gap == 0:
                    end = next(previous[2][1] for previous in clusters[-1] if previous[1] == last_end)
                    assert np.array_equal(end, item[2][0]), "Unverified cap-touch relation"
                clusters[-1].append(item)
        for cluster in clusters:
            begin = min(cluster, key=lambda item: item[0])[2][0]
            end = max(cluster, key=lambda item: item[1])[2][1]
            axis = (end-begin)/np.linalg.norm(end-begin)
            result = dict(center=((begin+end)/2).tolist(), axis=axis.tolist(), height=float(np.linalg.norm(end-begin)),
                radius=cluster[0][3]["radius"], endpoints=np.stack((begin, end)).tolist())
            sources = [dict(source_id=item[3].get("source_id", item[3]["id"]), image=item[3].get("image", [0, 0, 0]),
                interval=[item[0], item[1]], original_height=item[3]["height"],
                endpoint_reconstruction_error=item[3].get("endpoint_reconstruction_error", 0)) for item in cluster]
            gaps = [float(cluster[i][0]-max(item[1] for item in cluster[:i])) for i in range(1, len(cluster))]
            assert all(gap <= 0 for gap in gaps)
            reports.append(dict(output=len(merged), sources=sources, axial_gaps=gaps,
                merged_endpoints=result["endpoints"], collinearity_cross_products_exactly_zero=True,
                radius_bitwise_equal=True, no_positive_gap_bridged=True,
                equivalence_scope="Anchored collinear cylinder union identity; original CAD rotation differences remain separately recorded and checked"))
            for item in cluster:
                for point in item[2]:
                    if not np.array_equal(point, begin) and not np.array_equal(point, end):
                        former_caps.append(dict(center=point.tolist(), axis=axis.tolist(), radius=result["radius"]))
            merged.append(result)
    return merged, reports, former_caps


def field(records, *, preserve_rotations=False):
    terms = []
    for record in records:
        axis = np.asarray(record["axis"])
        if preserve_rotations:
            matrix = np.asarray(record["rotation"]).T
        else:
            coordinate = np.eye(3)[int(abs(axis).argmin())]
            transverse = np.cross(axis, coordinate)
            transverse /= np.linalg.norm(transverse)
            matrix = np.stack((axis, transverse, np.cross(axis, transverse)))
        terms.append((record["center"], matrix.tolist(), record["radius"], record["height"]/2))

    def evaluate(x, y, z):
        result = None
        for center, matrix, radius, half in terms:
            dx, dy, dz = x-center[0], y-center[1], z-center[2]
            local = [row[0]*dx+row[1]*dy+row[2]*dz for row in matrix]
            radial = np.sqrt(local[1]*local[1]+local[2]*local[2])-radius
            axial = np.abs(local[0])-half
            outside = np.sqrt(np.maximum(radial, 0)**2+np.maximum(axial, 0)**2)
            value = outside+np.minimum(np.maximum(radial, axial), 0)
            result = value if result is None else np.minimum(result, value)
        return result
    return evaluate


def cad_samples(cad, part):
    """Evaluate OCCT surfaces/curves directly, retaining every trimmed face."""
    from OCP.BRep import BRep_Tool
    from OCP.BRepAdaptor import BRepAdaptor_Surface, BRepAdaptor_Curve
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.BRepTools import BRepTools
    from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX, TopAbs_IN, TopAbs_ON
    from OCP.TopExp import TopExp
    from OCP.TopTools import TopTools_IndexedMapOfShape
    from OCP.TopoDS import TopoDS
    from OCP.gp import gp_Pnt2d

    coordinates = lambda p: [p.X(), p.Y(), p.Z()]
    faces, edges, vertices = (TopTools_IndexedMapOfShape() for _ in range(3))
    for kind, indexed in ((TopAbs_FACE, faces), (TopAbs_EDGE, edges), (TopAbs_VERTEX, vertices)):
        TopExp.MapShapes_s(cad.wrapped, kind, indexed)
    samples, rows = [], []
    for index in range(1, faces.Extent()+1):
        face = TopoDS.Face_s(faces.FindKey(index))
        bounds = BRepTools.UVBounds_s(face)
        assert np.isfinite(bounds).all()
        u0, u1, v0, v1 = bounds
        adaptor = BRepAdaptor_Surface(face, True)
        accepted = []
        for u in u0+(np.arange(32)+.5)/32*(u1-u0):
            for v in v0+(np.arange(32)+.5)/32*(v1-v0):
                classifier = BRepClass_FaceClassifier(face, gp_Pnt2d(float(u), float(v)), 1e-10)
                if classifier.State() in (TopAbs_IN, TopAbs_ON):
                    accepted.append(coordinates(adaptor.Value(float(u), float(v))))
        assert accepted, f"CAD face {index} has no accepted exact surface samples"
        samples.extend(accepted)
        rows.append(dict(face=index, uv_bounds=list(bounds), accepted=len(accepted), surface_type=str(adaptor.GetType())))
        assert len(samples) <= 200000, "Fixed CAD query sample resource cap exceeded"
        part.update(index/max(faces.Extent(), 1)*.8, "Sampling exact trimmed CAD surfaces", faces=index, samples=len(samples))
    curve_samples, edge_rows = [], []
    for index in range(1, edges.Extent()+1):
        edge = TopoDS.Edge_s(edges.FindKey(index))
        if BRep_Tool.Degenerated_s(edge):
            edge_rows.append(dict(edge=index, degenerate=True))
            continue
        adaptor = BRepAdaptor_Curve(edge)
        a, b = adaptor.FirstParameter(), adaptor.LastParameter()
        assert np.isfinite([a, b]).all()
        curve_samples.extend(coordinates(adaptor.Value(float(t))) for t in np.linspace(a, b, 129))
        edge_rows.append(dict(edge=index, degenerate=False, bounds=[a, b], samples=129, curve_type=str(adaptor.GetType())))
    vertex_samples = [coordinates(BRep_Tool.Pnt_s(TopoDS.Vertex_s(vertices.FindKey(i)))) for i in range(1, vertices.Extent()+1)]
    all_points = np.vstack((samples, curve_samples, vertex_samples))
    assert len(all_points) <= 300000
    on_cap = np.any(abs(abs(all_points)-.5) <= 1e-8, axis=1)
    part.update(1, "Exact CAD face/edge/vertex samples frozen", samples=len(all_points), edges=edges.Extent())
    return dict(walls=all_points[~on_cap], caps=all_points[on_cap], edges=np.asarray(curve_samples), vertices=np.asarray(vertex_samples)), dict(
        faces=rows, edges=edge_rows, vertices=len(vertex_samples), exact_surface_points=len(samples),
        method="32x32 UV midpoint grid inside trimmed faces; 129 exact curve points per nondegenerate edge; all CAD vertices",
        scope="Points lie on OCCT geometry; finite sample coverage is not a Hausdorff certificate")
