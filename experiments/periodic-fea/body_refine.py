"""Prototype nested interior-edge refinement using a coarse six-load guide.

This reads frozen, already audited artifacts. Adaptation time is measured, but
the complete direct recipe is only estimated from archived seed/guide timings.
No generation speed claim is valid until the whole recipe is run live.
"""

import json
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from homogenization import element_geometry
from residual_analysis import faces, stress_and_strain
from run import digest

TET_EDGES = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
TRI_EDGES = np.array([[0, 1], [1, 2], [2, 0]])


def replay_mesh(case, cells):
    """Reuse an admitted mesh to compare solvers without meshing cost."""
    directory = Path(__file__).parent / "body-seed"
    row = json.loads((directory / "result.json").read_text())
    definition = json.loads((directory / "input-provenance.json").read_text())["specimen"]
    if row["case"] != case.id or row["cells"] != cells or definition["offset"] != case.offset or definition["function"] != case.function:
        raise ValueError("Replay input differs from the requested specimen")
    if digest(directory / "mesh.npz") != row["mesh_sha256"]:
        raise ValueError("Replay mesh changed")
    with np.load(directory / "mesh.npz") as mesh:
        return mesh["points"], mesh["tetrahedra"], {
            "generation_s": 0., "mesh_s": 0., "setup_s": 0., "prototype_only": True,
            "generation_timing_scope": "Saved mesh replay; no meshing performed",
            "input_mesh_sha256": row["mesh_sha256"],
        }


def verify_preserved_boundary(seed_points, seed_surface, points, surface):
    """Compare retained coordinates and oriented triangles, independent of row order."""
    if len(points) < len(seed_points) or not np.array_equal(points[:len(seed_points)], seed_points):
        raise ValueError("Body refinement changed an original vertex")

    def ordered(triangles):
        indices = (triangles.argmin(axis=1)[:, None] + np.arange(3)) % 3
        rotated = np.take_along_axis(triangles, indices, axis=1)
        return rotated[np.lexsort(rotated.T[::-1])]

    if not np.array_equal(ordered(seed_surface), ordered(surface)):
        raise ValueError("Body refinement changed an oriented boundary triangle")
    return {"retained_vertices": len(seed_points), "unchanged_oriented_boundary_triangles": len(surface)}


def guide_score(directory):
    """Worst mixed-load traction residual per cell, normalized by guide energy."""
    from types import SimpleNamespace

    row = json.loads((directory / "result.json").read_text())
    archive = np.load(directory / "mesh.npz")
    if digest(directory / "mesh.npz") != row["mesh_sha256"]:
        raise ValueError("Guide mesh changed")
    points, tets = archive["points"], archive["tetrahedra"]
    classes = np.load(directory / "periodic.npz")["classes"]
    topology = faces(points, tets, classes)
    left, right, interior, normal, weight, _ = topology
    volume, gradients = element_geometry(SimpleNamespace(points=points, tetrahedra=tets))
    jumps = np.empty((len(left), 3, 6))
    columns = []
    for load, record in enumerate(row["homogenization"]["displacement_fields"]):
        path = directory / record["artifact"]
        if digest(path) != record["sha256"]:
            raise ValueError("Guide field changed")
        field = np.load(path)
        if not np.array_equal(field["engineering"], np.eye(6)[load]):
            raise ValueError("Wrong guide load convention")
        stress, _ = stress_and_strain(field["displacement"], tets, gradients)
        mean = np.einsum("t,tij->ij", volume, stress)
        columns.append(mean[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]])
        jumps[:, :, load] = np.einsum("fij,fj->fi", stress[left], normal)
        jumps[interior, :, load] -= np.einsum("fij,fj->fi", stress[right[interior]], normal[interior])
    stiffness = np.asarray(row["homogenization"]["stiffness_engineering_voigt"])
    if np.linalg.norm(np.asarray(columns).T - stiffness) / np.linalg.norm(stiffness) > 1e-10:
        raise ValueError("Guide fields do not reproduce stiffness with the expected constitutive law")
    factor = np.linalg.cholesky((stiffness + stiffness.T) / 2)
    normalized = jumps @ np.linalg.inv(factor).T
    face_gram = weight[:, None, None] * np.einsum("fci,fcj->fij", normalized, normalized)
    cell_gram = np.zeros((len(tets), 6, 6))
    np.add.at(cell_gram, left, face_gram * np.where(interior, .5, 1)[:, None, None])
    np.add.at(cell_gram, right[interior], face_gram[interior] / 2)
    score = np.maximum(np.linalg.eigvalsh(cell_gram)[:, -1], 0)
    # Strip the leading h^5 scale of a smooth P1 traction-jump indicator before
    # transferring its spatial density to another mesh. This remains a heuristic.
    density = score / volume ** (5 / 3)
    return points[tets].mean(axis=1), density, row


def split_interior_edges(points, tets, surface, density, max_new_nodes, minimum_quality=.1):
    """Split complete edge stars in disjoint cavities; preserve every boundary face."""
    count = len(points)
    all_edges = np.sort(tets[:, TET_EDGES].reshape(-1, 2), axis=1)
    edges, inverse, counts = np.unique(all_edges, axis=0, return_inverse=True, return_counts=True)
    boundary_edges = np.sort(surface[:, TRI_EDGES].reshape(-1, 2), axis=1)
    boundary_codes = np.unique(boundary_edges[:, 0] * count + boundary_edges[:, 1])
    eligible = ~np.isin(edges[:, 0] * count + edges[:, 1], boundary_codes)
    midpoint = points[edges].mean(axis=1)
    eligible &= np.all((midpoint > -.5 + 1e-9) & (midpoint < .5 - 1e-9), axis=1)
    length2 = np.sum((points[edges[:, 0]] - points[edges[:, 1]])**2, axis=1)
    priority = density * length2
    ordering = np.argsort(inverse, kind="stable")
    starts = np.r_[0, np.cumsum(counts)]
    cells = np.repeat(np.arange(len(tets)), 6)
    removed = np.zeros(len(tets), dtype=bool)
    new_points, additions = [], []
    rejected_quality = 0
    candidates = np.flatnonzero(eligible & (priority > 1e-20))
    for ei in candidates[np.argsort(-priority[candidates], kind="stable")]:
        if len(new_points) >= max_new_nodes:
            break
        incident = cells[ordering[starts[ei]:starts[ei + 1]]]
        if removed[incident].any():
            continue
        a, b = edges[ei]
        children = np.tile(tets[incident], (2, 1))
        positions = np.r_[np.argmax(tets[incident] == a, axis=1), np.argmax(tets[incident] == b, axis=1)]
        xyz = points[children].copy()
        xyz[np.arange(len(children)), positions] = midpoint[ei]
        determinant = np.linalg.det(xyz[:, 1:] - xyz[:, :1])
        edge_sum = np.sum((xyz[:, TET_EDGES[:, 0]] - xyz[:, TET_EDGES[:, 1]])**2, axis=(1, 2))
        quality = np.sqrt(432 * determinant**2 / edge_sum**3)
        if determinant.min() <= 0 or quality.min() < minimum_quality:
            rejected_quality += 1
            continue
        children[np.arange(len(children)), positions] = count + len(new_points)
        new_points.append(midpoint[ei])
        additions.append(children)
        removed[incident] = True
    refined_points = np.vstack([points, np.asarray(new_points).reshape(-1, 3)])
    refined_tets = np.vstack([tets[~removed], *additions]) if additions else tets.copy()
    return refined_points, refined_tets, {"inserted_body_nodes": len(new_points),
                                          "bisected_parent_tetrahedra": int(removed.sum()),
                                          "quality_rejected_edge_stars": rejected_quality}


def generate(case, cells, config, directory, checkpoint):
    source = Path(__file__).parent
    seed_dir, guide_dir = source / "body-seed", source / "body-guide"
    started = time.perf_counter()
    checkpoint(stage="Coarse mixed-load residual map")
    centers, density, guide = guide_score(guide_dir)
    seed = json.loads((seed_dir / "result.json").read_text())
    if seed["case"] != case.id or guide["case"] != case.id:
        raise ValueError("Seed/guide geometry differs from requested case")
    definitions = [json.loads((folder / "input-provenance.json").read_text())["specimen"] for folder in [seed_dir, guide_dir]]
    if definitions[0] != definitions[1] or definitions[0]["offset"] != case.offset or definitions[0]["function"] != case.function:
        raise ValueError("Seed/guide field definitions differ from the requested specimen")
    mesh = np.load(seed_dir / "mesh.npz")
    if digest(seed_dir / "mesh.npz") != seed["mesh_sha256"]:
        raise ValueError("Seed mesh changed")
    points, tets, surface = mesh["points"], mesh["tetrahedra"], mesh["surface"]
    all_edges = np.unique(np.sort(tets[:, TET_EDGES].reshape(-1, 2), axis=1), axis=0)
    query = (points[all_edges].mean(axis=1) + .5) % 1
    tree = cKDTree((centers + .5) % 1, boxsize=1)
    _, nearest = tree.query(query, workers=1)
    checkpoint(stage="Nested body-edge bisection")
    points, tets, metadata = split_interior_edges(points, tets, surface, density[nearest],
                                                 config["body_refine_nodes"], config["minimum_quality"])
    if len(tets) > config["max_tetrahedra"]:
        raise ValueError("Body-refinement tetrahedron budget exceeded")
    elapsed = time.perf_counter() - started
    return points, tets, {
        **metadata, "generation_s": elapsed, "mesh_s": elapsed, "setup_s": 0.,
        "generation_timing_scope": "Only Python adaptation of saved artifacts; not direct generation",
        "prototype_only": True, "guide_cells": guide["cells"], "seed_cells": seed["cells"],
        "archived_guide_generation_s": guide["generation_s"],
        "archived_guide_fea_s": guide["homogenization"]["seconds"],
        "archived_seed_generation_s": seed["generation_s"],
        "estimated_live_recipe_s": elapsed + guide["generation_s"] + guide["homogenization"]["seconds"] + seed["generation_s"],
        "seed_tensor": seed["homogenization"]["stiffness_engineering_voigt"],
        "seed_volume": seed["solid_fraction"],
        "body_marker": "Coarse worst mixed-load traction residual, transferred by periodic nearest centroid",
    }
