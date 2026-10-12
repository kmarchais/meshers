"""Experimental HXT volume construction from a frozen periodic native surface."""

import time

import numpy as np


def generate(case, cells, config, directory, checkpoint):
    import gmsh
    from scipy import sparse
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
    from microgen_benchmark.backends import generate as native_generate
    from microgen_benchmark.geometry import boundary, canonical_surface

    points, tets, seed = native_generate(case, "meshers-pass1", cells, config, checkpoint)
    started = time.perf_counter()
    surface = boundary(points, tets)
    surface_points, faces, _, _ = canonical_surface(points, surface)
    # A periodic material can be disconnected in the cut cube. Each connected
    # solid must be declared as a volume, including small corner components.
    left, right = tets[:, [0, 0, 0]].ravel(), tets[:, [1, 2, 3]].ravel()
    _, components = connected_components(
        sparse.coo_matrix((np.ones(len(left)), (left, right)), shape=(len(points), len(points))),
        directed=False,
    )
    surface_components = components[np.unique(surface)]
    face_components = surface_components[faces[:, 0]]
    assert np.all(surface_components[faces] == face_components[:, None])
    # HXT distinguishes the two Gyroid void networks through boundary-surface
    # identifiers. One tag for every triangle gives both voids the same BRep.
    coordinates = surface_points[faces]
    patches = np.full(len(faces), -1, dtype=np.int64)
    for axis in range(3):
        for side, value in enumerate([-case.scale / 2, case.scale / 2]):
            cap = np.all(np.abs(coordinates[:, :, axis] - value) < 1e-12 * case.scale, axis=1)
            patches[cap] = 2 * axis + side
    wall = patches < 0
    if wall.any():
        values = case.field(*coordinates[wall].reshape(-1, 3).T).reshape(-1, 3)
        if not np.all(np.abs(np.abs(values) - 1) < 1e-8):
            raise ValueError("Unexpected native wall vertices away from band isosurfaces")
        if not np.all(np.sign(values) == np.sign(values[:, :1])):
            raise ValueError("A native wall face crosses both band isosurfaces")
        patches[wall] = 6 + (values.mean(axis=1) > 0)
    options = {
        "General.NumThreads": 1,
        "Mesh.MaxNumThreads1D": 1,
        "Mesh.MaxNumThreads2D": 1,
        "Mesh.MaxNumThreads3D": 1,
        "Mesh.Algorithm3D": 10,
        "Mesh.MeshOnlyEmpty": 1,
        "Mesh.MeshSizeMin": case.scale / cells,
        "Mesh.MeshSizeMax": case.scale / cells,
        "Mesh.MeshSizeFromPoints": 0,
        "Mesh.MeshSizeFromCurvature": 0,
        "Mesh.MeshSizeExtendFromBoundary": 0,
        "Mesh.Optimize": 1,
        "Mesh.OptimizeNetgen": 0,
        "Mesh.ElementOrder": 1,
    }
    np.savez_compressed(directory / "native-surface.npz", points=surface_points, triangles=faces,
                        patches=patches, components=face_components)
    checkpoint(stage="HXT volume construction", surface_triangles=len(faces))
    gmsh.initialize([], readConfigFiles=False)
    log = []
    try:
        gmsh.logger.start()
        gmsh.model.add("periodic-native-surface")
        for name, value in options.items():
            gmsh.option.setNumber(name, value)
        shells, registered = [], np.zeros(len(surface_points), dtype=bool)
        for component in np.unique(face_components):
            surface_tags = []
            for patch in np.unique(patches[face_components == component]):
                surface_tag = gmsh.model.addDiscreteEntity(2)
                surface_tags.append(surface_tag)
                triangles = faces[(face_components == component) & (patches == patch)]
                used = np.unique(triangles)
                new = used[~registered[used]]
                gmsh.model.mesh.addNodes(2, surface_tag, new + 1, surface_points[new].ravel())
                registered[new] = True
                gmsh.model.mesh.addElementsByType(surface_tag, 2, [], (triangles + 1).ravel())
            shells.append(surface_tags)
        assert registered.all()
        # Material-wall/cap interfaces must be explicit curves, so boundary
        # recovery treats their vertices as part of the prescribed topology.
        gmsh.model.mesh.createTopology(makeSimplyConnected=False, exportDiscrete=True)
        gmsh.model.mesh.reclassifyNodes()
        for surface_tags in shells:
            shell = gmsh.model.geo.addSurfaceLoop(surface_tags)
            gmsh.model.geo.addVolume([shell])
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(3)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        tags, xyz = np.asarray(tags), np.asarray(xyz).reshape(-1, 3)
        order = np.argsort(tags)
        tags, xyz = tags[order], xyz[order]
        types, _, connectivity = gmsh.model.mesh.getElements(3)
        if list(types) != [4]:
            raise ValueError(f"Expected only linear tetrahedra, received Gmsh types {types}")
        element_tags = np.asarray(connectivity[0])
        mapped = np.searchsorted(tags, element_tags)
        if not np.array_equal(tags[mapped], element_tags):
            raise ValueError("Gmsh tetrahedra reference missing nodes")
        result_tets = mapped.reshape(-1, 4)
        if len(result_tets) > config["max_tetrahedra"]:
            raise ValueError("HXT output exceeds declared tetrahedron budget")
        # Verify the whole boundary, not just the periodic caps. HXT is free to
        # reorder nodes, but may not change the input surface in this control.
        from microgen_benchmark.geometry import canonical_mesh
        xyz, result_tets, _, _ = canonical_mesh(xyz, result_tets)
        np.savez_compressed(directory / "hxt-volume.npz", points=xyz, tetrahedra=result_tets)
        output_surface = boundary(xyz, result_tets)
        output_points, output_faces, _, _ = canonical_surface(xyz, output_surface)
        distance, match = cKDTree(surface_points).query(output_points)
        if (len(output_points) != len(surface_points) or len(np.unique(match)) != len(match)
                or distance.max() > 1e-12 * case.scale):
            raise ValueError(f"HXT changed boundary vertices: {len(surface_points)} input, "
                             f"{len(output_points)} output, maximum nearest-input distance {distance.max():.6g}")
        expected = np.sort(faces, axis=1)
        actual = np.sort(match[output_faces], axis=1)
        expected = expected[np.lexsort(expected.T[::-1])]
        actual = actual[np.lexsort(actual.T[::-1])]
        if not np.array_equal(expected, actual):
            raise ValueError("HXT changed boundary triangulation")
        construction_s = time.perf_counter() - started
        return xyz, result_tets, {
            **seed,
            "mesh_s": seed["mesh_s"] + construction_s,
            "generation_s": seed["generation_s"] + construction_s,
            "native_seed_mesh_s": seed["mesh_s"],
            "native_seed_tetrahedra": len(tets),
            "surface_triangles": len(faces),
            "material_components": len(shells),
            "surface_patch_triangle_counts": {str(p): int((patches == p).sum()) for p in np.unique(patches)},
            "hxt_construction_and_boundary_checks_s": construction_s,
            "hxt_options": options,
            "hxt_boundary_unchanged": True,
            "hxt_scope": "Full native seed cost included; no surface-only speed claim",
        }
    finally:
        log = gmsh.logger.get()
        gmsh.finalize()
        (directory / "gmsh.log").write_text("\n".join(log) + "\n")
