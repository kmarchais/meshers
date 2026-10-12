"""Experimental HXT volume construction from a frozen periodic native surface."""

import time

import numpy as np


def generate(case, cells, config, directory, checkpoint):
    import gmsh
    from scipy.spatial import cKDTree
    from microgen_benchmark.backends import generate as native_generate
    from microgen_benchmark.geometry import boundary, canonical_surface

    points, tets, seed = native_generate(case, "meshers-pass1", cells, config, checkpoint)
    started = time.perf_counter()
    surface = boundary(points, tets)
    surface_points, faces, _, _ = canonical_surface(points, surface)
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
    checkpoint(stage="HXT volume construction", surface_triangles=len(faces))
    gmsh.initialize([], readConfigFiles=False)
    log = []
    try:
        gmsh.logger.start()
        gmsh.model.add("periodic-native-surface")
        for name, value in options.items():
            gmsh.option.setNumber(name, value)
        surface_tag = gmsh.model.addDiscreteEntity(2)
        gmsh.model.mesh.addNodes(2, surface_tag, np.arange(1, len(surface_points) + 1), surface_points.ravel())
        gmsh.model.mesh.addElementsByType(surface_tag, 2, [], (faces + 1).ravel())
        shell = gmsh.model.geo.addSurfaceLoop([surface_tag])
        volume = gmsh.model.geo.addVolume([shell])
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(3)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        tags, xyz = np.asarray(tags), np.asarray(xyz).reshape(-1, 3)
        order = np.argsort(tags)
        tags, xyz = tags[order], xyz[order]
        types, _, connectivity = gmsh.model.mesh.getElements(3, volume)
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
        output_surface = boundary(xyz, result_tets)
        output_points, output_faces, _, _ = canonical_surface(xyz, output_surface)
        distance, match = cKDTree(surface_points).query(output_points)
        if (len(output_points) != len(surface_points) or len(np.unique(match)) != len(match)
                or distance.max() > 1e-12 * case.scale):
            raise ValueError("HXT changed boundary vertices")
        expected = np.sort(faces, axis=1)
        actual = np.sort(match[output_faces], axis=1)
        expected = expected[np.lexsort(expected.T[::-1])]
        actual = actual[np.lexsort(actual.T[::-1])]
        if not np.array_equal(expected, actual):
            raise ValueError("HXT changed boundary triangulation")
        construction_s = time.perf_counter() - started
        np.savez_compressed(directory / "native-surface.npz", points=surface_points, triangles=faces)
        return xyz, result_tets, {
            **seed,
            "mesh_s": seed["mesh_s"] + construction_s,
            "generation_s": seed["generation_s"] + construction_s,
            "native_seed_mesh_s": seed["mesh_s"],
            "native_seed_tetrahedra": len(tets),
            "surface_triangles": len(faces),
            "hxt_construction_and_boundary_checks_s": construction_s,
            "hxt_options": options,
            "hxt_boundary_unchanged": True,
            "hxt_scope": "Full native seed cost included; no surface-only speed claim",
        }
    finally:
        log = gmsh.logger.get()
        gmsh.finalize()
        (directory / "gmsh.log").write_text("\n".join(log) + "\n")
