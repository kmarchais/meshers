"""Read-only derivative certificate for the archived CAD restoration control.

No optimizer starts and no mesh displacement is applied. The linear systems
diagnose the archived solution and one adjacent complete-star expansion.
"""

import argparse
import json
from pathlib import Path
import shutil

import cad_ridge_quality_restore as restore
import numpy as np


ACTIVE_TOLERANCE = 1e-9
DUAL_TOLERANCE = 1e-8


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=restore.w.HERE/"runs/macbook-cad-ridge-quality-restore-v1")
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    args.source, args.run = args.source.resolve(), args.run.resolve()
    args.run.mkdir(parents=True, exist_ok=False)
    saved = json.loads((args.source/"result.json").read_text())
    assert saved["complete"] and saved["ordinary_admitted"]
    assert restore.w.digest(Path(restore.__file__)) == saved["inputs"][str(Path(restore.__file__).resolve())]
    original_file = Path(saved["source"])/"mesh.npz"
    proposal_file = Path(saved["proposals"])/"frozen-proposals.json"
    for path in (original_file, proposal_file):
        assert restore.w.digest(path) == saved["inputs"][str(path)]
    assert restore.w.digest(args.source/"candidate.npz") == saved["candidate_sha256"]
    sources = {Path(__file__).resolve(), Path(restore.__file__).resolve(), Path(restore.w.__file__).resolve(),
               Path(restore.cavity.__file__).resolve()}
    inputs = {original_file, proposal_file, args.source/"result.json", args.source/"candidate.npz",
              args.source/"repaired-before-restoration.npz", args.source/"quality-vectors.npz"}
    hashes = {str(path): restore.w.digest(path) for path in sorted(sources | inputs)}
    for path in sources:
        shutil.copyfile(path, args.run/path.name)
    restore.w.save(args.run/"manifest.json", dict(inputs=hashes, read_only=True,
        optimization_performed=False, candidate_mesh_written=False, active_tolerance=ACTIVE_TOLERANCE,
        positive_dual_tolerance=DUAL_TOLERANCE, numerical_threads=1,
        scope="Small least-squares derivative audit of saved arrays; no new optimizer or geometry candidate"))
    with np.load(original_file) as archive:
        points, cells, surface = (archive[k] for k in ("points", "tetrahedra", "surface"))
    with np.load(args.source/"repaired-before-restoration.npz") as archive:
        repaired = archive["points"]
    with np.load(args.source/"candidate.npz") as archive:
        moved = archive["points"]
    proposal = json.loads(proposal_file.read_text())["frozen"][0]
    plan = restore.inventory(points, cells, surface, proposal)
    model = restore.Restoration(repaired, points, cells, plan)
    variables = np.asarray(saved["optimizer"]["variables"])
    residual, jacobian, quality, _, direction = model.evaluate(variables)
    ordering = np.argsort(quality, kind="stable")
    original_quality = restore.quality_gradient(points[model.cells])[0]
    original_ordering = np.argsort(original_quality, kind="stable")
    active = np.flatnonzero(residual[:len(quality)] < ACTIVE_TOLERANCE)
    objective_gradient = np.r_[np.zeros(model.dimension), 1.]
    dual = np.linalg.lstsq(jacobian[active].T, objective_gradient, rcond=None)[0]
    stationarity = float(np.linalg.norm(jacobian[active].T@dual-objective_gradient))
    assert stationarity < 1e-12 and dual.min() > -1e-12
    records = []
    for rank, multiplier in zip(active, dual, strict=True):
        index = ordering[rank]
        records.append(dict(rank=int(rank), tetrahedron=int(plan["expanded"][index]),
            target_tetrahedron=int(plan["expanded"][original_ordering[rank]]), quality=float(quality[index]),
            target=float(model.target[rank]), residual=float(residual[rank]), dual_weight=float(multiplier),
            gradient_norm=float(np.linalg.norm(jacobian[rank, :-1])), corners=model.cells[index].tolist(),
            free_corners=[int(n) for n in model.cells[index] if n in model.free],
            quality_gap_below=float(quality[index]-quality[ordering[rank-1]]) if rank else None,
            quality_gap_above=float(quality[ordering[rank+1]]-quality[index]) if rank+1 < len(quality) else None))
    selected = dual > DUAL_TOLERANCE
    support, weights = active[selected], dual[selected]
    supported_cells = model.cells[ordering[support]]
    _, gradients, _ = restore.quality_gradient(moved[supported_cells])
    reactions = {}
    for weight, cell, gradient in zip(weights, supported_cells, gradients, strict=True):
        for node, vector in zip(cell, gradient, strict=True):
            reactions[int(node)] = reactions.get(int(node), np.zeros(3))+weight*vector
    physical_nodes, free_nodes = set(map(int, surface.ravel())), set(map(int, model.free))
    reaction_rows = []
    for node, vector in reactions.items():
        star = np.flatnonzero((cells[:, :4] == node).any(1))
        neighbors = np.setdiff1d(np.unique(cells[star, :4]), [node])
        length = float(np.linalg.norm(points[neighbors]-points[node], axis=1).min())
        reaction_rows.append(dict(node=node,
            role="free-interior" if node in free_nodes else "fixed-wall" if node in physical_nodes else "fixed-interior",
            reaction=vector.tolist(), reaction_norm=float(np.linalg.norm(vector)),
            scaled_norm=float(np.linalg.norm(vector)*.1*length), original_shortest_edge=length,
            star_tetrahedra=len(star), expanded_tetrahedra=len(np.union1d(plan["expanded"], star)),
            additional_tetrahedra=len(np.setdiff1d(star, plan["expanded"]))))
    reaction_rows.sort(key=lambda row: (-row["scaled_norm"], row["node"]))
    chosen = next(row for row in reaction_rows if row["role"] == "fixed-interior")
    node = chosen["node"]
    star = np.flatnonzero((cells[:, :4] == node).any(1))
    expanded = dict(plan, interior=np.r_[plan["interior"], node],
        expanded=np.union1d(plan["expanded"], star), shortest=np.r_[plan["shortest"], chosen["original_shortest_edge"]])
    boundary = restore.cavity.boundary_keys(restore.cavity.boundary(cells[expanded["expanded"]]))
    physical = restore.cavity.boundary_keys(surface)
    interface_nodes = np.unique([n for key, face in boundary.items() if key not in physical for n in face])
    assert not np.isin(expanded["interior"], interface_nodes).any()
    assert not np.isin(expanded["interior"], surface).any()
    expanded_model = restore.Restoration(repaired, points, cells, expanded)
    expanded_variables = np.r_[variables[:-1], [0., 0., 0.], variables[-1]]
    er, ej, eq, _, _ = expanded_model.evaluate(expanded_variables)
    expanded_variables[-1] = float((expanded_model.target-np.sort(eq)).max())
    er, ej, eq, _, _ = expanded_model.evaluate(expanded_variables)
    ea = np.flatnonzero(er[:len(eq)] < ACTIVE_TOLERANCE)
    matrix = ej[ea, :-1]
    descent = np.linalg.lstsq(matrix, np.ones(len(ea)), rcond=None)[0]
    equation_error = float(abs(matrix@descent-1).max())
    descent /= np.linalg.norm(descent.reshape(-1, 3), axis=1).max()
    derivatives = matrix@descent
    assert equation_error < 1e-10 and derivatives.min() > 0
    quality_rows = {name: saved["checks"][name] for name in ("original_quality", "repaired_quality", "restored_quality")}
    np.savez_compressed(args.run/"derivative-certificate.npz", active_ranks=active, dual_weights=dual,
        active_jacobian=jacobian[active], expanded_active_ranks=ea, expanded_quality_jacobian=matrix,
        normalized_direction_not_applied=descent, expanded_free_corners=expanded["interior"],
        expanded_tetrahedra=expanded["expanded"])
    result = dict(inputs=hashes, complete=True, active_constraints=records,
        positive_dual_ranks=support.tolist(), positive_dual_tetrahedra=plan["expanded"][ordering[support]].tolist(),
        positive_dual_weights=weights.tolist(), positive_dual_sum=float(weights.sum()),
        stationarity_residual=stationarity,
        smallest_sorted_quality_gap=float(np.diff(np.sort(quality)).min()),
        maximum_radius_fraction=float(np.linalg.norm(direction, axis=1).max()),
        support_gradient_rank=int(np.linalg.matrix_rank(jacobian[support, :-1])),
        fixed_corner_reactions=reaction_rows, quality=quality_rows,
        proposal0_all_reported_quantiles_nonworse=all(quality_rows["repaired_quality"][key] >= value-1e-13
            for key, value in quality_rows["original_quality"].items()),
        ordinary_admitted=saved["ordinary_admitted"], research_sorted_quality_feasible=saved["research_sorted_quality_feasible"],
        interpretation="Local first-order minimax obstruction in the current cavity, not a sorting tie or active radius limit; no global impossibility claim",
        next_control=dict(selected_corner=node, selection="Largest original-edge-scaled dual reaction among fixed interior corners",
            free_corners=expanded["interior"].tolist(), tetrahedra=len(expanded["expanded"]),
            additional_tetrahedra=chosen["additional_tetrahedra"], maximum_displacement_radius=float(expanded_model.radii[-1]),
            cap_clearance=float(.5-abs(points[np.unique(cells[expanded["expanded"]])]).max()),
            interface_contains_no_free_corners=True, at_archived_point_slack=float(expanded_variables[-1]),
            active_ranks=ea.tolist(), active_tetrahedra=expanded["expanded"][np.argsort(eq, kind="stable")[ea]].tolist(),
            active_gradient_rank=int(np.linalg.matrix_rank(matrix)), active_count=len(ea),
            directional_quality_derivative_minimum=float(derivatives.min()),
            directional_quality_derivative_maximum=float(derivatives.max()),
            least_squares_equation_error=equation_error,
            smallest_inactive_margin=float(er[:len(eq)][er[:len(eq)] >= ACTIVE_TOLERANCE].min()),
            direction_applied=False, geometry_or_quality_feasibility_proven=False,
            scope="One 224-tet cavity with corner 3419 newly free; physical repaired boundary stays fixed; no solve has been run"))
    assert all(restore.w.digest(Path(path)) == digest for path, digest in hashes.items())
    result["input_hashes_verified"] = True
    restore.w.save(args.run/"result.json", result)
    print(json.dumps(dict(stationarity_residual=stationarity, positive_dual_constraints=len(support),
        selected_corner=node, expanded_tetrahedra=len(expanded["expanded"]),
        minimum_directional_derivative=float(derivatives.min()), mesh_modified=False)))


if __name__ == "__main__":
    main()
