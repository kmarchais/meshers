"""Conforming Tet10 edge-star bisection that retains each parent curved map."""

import numpy as np

EDGES = np.array(((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)))
FACES = np.array(((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)))


def child_maps(offset):
    first, second = EDGES[offset]
    maps = np.repeat(np.eye(4)[None], 2, axis=0)
    maps[0, first] = maps[1, second] = (np.eye(4)[first] + np.eye(4)[second]) / 2
    return maps


def positive_quality(coordinates):
    determinant = np.linalg.det(coordinates[:, 1:] - coordinates[:, :1])
    squared = sum(((coordinates[:, first] - coordinates[:, second])**2).sum(axis=1)
                  for first, second in EDGES)
    quality = (12 * (abs(determinant) / 2)**(2 / 3) / squared)**1.5
    return np.where(determinant > 0, quality, 0)


def refine(points, cells, surface, fields, *, periodic=True, minimum_quality=.1,
           maximum_nodes=250000, maximum_tets=250000):
    from fedoo.lib_elements.tetrahedron import Tet10
    from facet_incidence_audit import audit

    assert cells.shape[1] == 10 and surface.shape[1] == 6
    assert fields.shape == (6, len(points), 3) and np.isfinite(fields).all()
    assert len(points) <= maximum_nodes and len(cells) <= maximum_tets
    assert 0 < minimum_quality <= 1
    original = audit(points, cells, surface)
    assert original["oriented_facet_consistency_passed"] and not original["nonpositive_corner_determinants"]
    classes = np.arange(len(points))
    if periodic:
        from periodic import periodic_pairs
        _, classes, _ = periodic_pairs(points, surface)
    pairs = np.sort(cells[:, EDGES].reshape(-1, 2), axis=1)
    edges, first_rows, inverse = np.unique(pairs, axis=0, return_index=True, return_inverse=True)
    mids = cells[:, 4:].ravel()[first_rows]
    assert np.array_equal(mids[inverse], cells[:, 4:].ravel())
    assert not np.isin(mids, np.unique(cells[:, :4])).any()
    incidence = [[] for _ in edges]
    for row, edge in enumerate(inverse):
        incidence[edge].append(divmod(row, 6))
    groups = {}
    for index, edge in enumerate(edges):
        groups.setdefault(tuple(sorted(classes[edge])), []).append(index)
    lengths = np.linalg.norm(points[edges[:, 1]] - points[edges[:, 0]], axis=1)
    ordered = sorted(groups.values(), key=lambda group: (-float(lengths[group].max()), group[0]))
    allowed = np.ones((len(cells), 6), dtype=bool)
    for offset, endpoints in enumerate(EDGES):
        for endpoint in endpoints:
            coordinates = points[cells[:, :4]].copy()
            coordinates[:, endpoint] = points[cells[:, offset + 4]]
            allowed[:, offset] &= positive_quality(coordinates) >= minimum_quality + 1e-12
    selected = np.full(len(cells), -1, dtype=int)
    projected_nodes, projected_tets = len(points), len(cells)
    counts = {"accepted_edge_orbits": 0, "quality_rejected_orbits": 0,
              "overlapping_star_orbits": 0, "budget_rejected_orbits": 0}
    for group in ordered:
        incident = [item for edge in group for item in incidence[edge]]
        parents, offsets = np.asarray(incident).T
        if len(set(parents)) != len(parents) or (selected[parents] >= 0).any():
            counts["overlapping_star_orbits"] += 1
            continue
        if not allowed[parents, offsets].all():
            counts["quality_rejected_orbits"] += 1
            continue
        added = sum(2 + len(set(cells[[parent for parent, _ in incidence[edge]], :4].ravel())
                            - set(edges[edge])) for edge in group)
        if projected_nodes + added > maximum_nodes or projected_tets + len(parents) > maximum_tets:
            counts["budget_rejected_orbits"] += 1
            continue
        selected[parents] = offsets
        projected_nodes += added
        projected_tets += len(parents)
        counts["accepted_edge_orbits"] += 1
    element = Tet10()
    maps = [child_maps(offset) for offset in range(6)]
    weights = []
    for pair in maps:
        child_weights = []
        for mapping in pair:
            barycentric = np.concatenate((mapping, mapping[EDGES].mean(axis=1)))
            child_weights.append(element.shape_function(barycentric @ element.xi_nd[:4]))
        weights.append(child_weights)
    nodes, values = list(points.copy()), list(fields.transpose(1, 0, 2).copy())
    edge_nodes = {tuple(edge): int(mid) for edge, mid in zip(edges, mids, strict=True)}
    refined, parents, mappings = [], [], []
    shared_position_error, shared_field_error = 0., 0.
    for parent, cell in enumerate(cells):
        offset = selected[parent]
        if offset < 0:
            refined.append(cell.tolist())
            parents.append(parent)
            mappings.append(np.eye(4))
            continue
        for child, endpoint in enumerate(EDGES[offset]):
            corners = cell[:4].copy()
            corners[endpoint] = cell[4 + offset]
            positions = weights[offset][child] @ points[cell]
            displacement = np.einsum("an,lnc->alc", weights[offset][child], fields[:, cell])
            connectivity = corners.tolist()
            for local, (first, second) in enumerate(EDGES, start=4):
                key = tuple(sorted((int(corners[first]), int(corners[second]))))
                if key not in edge_nodes:
                    edge_nodes[key] = len(nodes)
                    nodes.append(positions[local])
                    values.append(displacement[local])
                node = edge_nodes[key]
                shared_position_error = max(shared_position_error, float(abs(nodes[node] - positions[local]).max()))
                shared_field_error = max(shared_field_error, float(abs(values[node] - displacement[local]).max()))
                connectivity.append(node)
            refined.append(connectivity)
            parents.append(parent)
            mappings.append(maps[offset][child])
    assert shared_position_error < 1e-12 and shared_field_error < 1e-11
    points_out, cells_out = np.asarray(nodes), np.asarray(refined)
    faces = cells_out[:, FACES].reshape(-1, 3)
    _, indices, multiplicities = np.unique(np.sort(faces, axis=1), axis=0, return_index=True, return_counts=True)
    boundary = faces[indices[multiplicities == 1]]
    surface_out = np.array([list(face) + [edge_nodes[tuple(sorted((int(face[first]), int(face[second]))))]
        for first, second in ((0, 1), (1, 2), (2, 0))] for face in boundary])
    assert len(points_out) <= projected_nodes <= maximum_nodes
    assert len(cells_out) == projected_tets <= maximum_tets
    assert len(np.unique(cells_out)) == len(points_out)
    consistency = audit(points_out, cells_out, surface_out)
    assert consistency["oriented_facet_consistency_passed"] and not consistency["nonpositive_corner_determinants"]
    counts.update(source_nodes=len(points), source_elements=len(cells), nodes=len(points_out),
        elements=len(cells_out), bisected_parent_elements=int((selected >= 0).sum()),
        shared_position_error=shared_position_error, shared_field_error=shared_field_error,
        projected_node_upper_bound=projected_nodes, facet_consistency=consistency)
    if periodic:
        _, _, counts["periodic"] = periodic_pairs(points_out, surface_out)
    return {"points": points_out, "tetrahedra": cells_out, "surface": surface_out,
            "displacement": np.asarray(values).transpose(1, 0, 2), "parents": np.asarray(parents),
            "parent_barycentric_maps": np.asarray(mappings), "diagnostics": counts}


def interpolation_check(source_points, source_cells, source_fields, refined):
    from fedoo.lib_elements.tetrahedron import Tet10

    element = Tet10()
    samples = np.random.default_rng(20261023).dirichlet(np.ones(4), 12)
    values = element.shape_function(samples @ element.xi_nd[:4])
    parents = refined["parents"]
    maximum_geometry_error, maximum_field_error = 0., 0.
    for start in range(0, len(parents), 512):
        selection = slice(start, start + 512)
        mapped = np.einsum("pa,cab->cpb", samples, refined["parent_barycentric_maps"][selection])
        parent_values = element.shape_function(mapped.reshape(-1, 4) @ element.xi_nd[:4]).reshape(-1, 12, 10)
        old = source_cells[parents[selection]]
        new = refined["tetrahedra"][selection]
        expected = np.einsum("cpn,cnd->cpd", parent_values, source_points[old])
        actual = np.einsum("pn,cnd->cpd", values, refined["points"][new])
        maximum_geometry_error = max(maximum_geometry_error, float(abs(actual - expected).max()))
        expected = np.einsum("cpn,lcnd->lcpd", parent_values, source_fields[:, old])
        actual = np.einsum("pn,lcnd->lcpd", values, refined["displacement"][:, new])
        maximum_field_error = max(maximum_field_error, float(abs(actual - expected).max()))
    assert maximum_geometry_error < 1e-12 and maximum_field_error < 1e-11
    return {"points_per_element": len(samples), "elements_checked": len(parents),
            "maximum_geometry_error": maximum_geometry_error, "maximum_field_error": maximum_field_error}
