"""Feasible target feet on the fixed level surface inside the study box."""

import numpy as np
from scipy.optimize import minimize


def constrained_foot(function, probe, level, locked, initial):
    from scorecard import gradient

    assert np.isfinite(probe).all()
    initial = initial if np.isfinite(initial).all() else probe.copy()
    scale = max(float(np.linalg.norm(initial-probe)), 1e-3)
    lower, upper = np.full(3, -.5), np.full(3, .5)
    lower[locked] = upper[locked] = np.sign(probe[locked])*.5
    start = np.clip(initial, lower, upper)
    if np.all(locked):
        foot, success, message, iterations = start, True, "all coordinates fixed", 0
    else:
        def residual(displacement):
            point = probe+scale*displacement
            return float(function(*point[None].T)[0]-level)/scale
        result = minimize(lambda displacement: float(displacement @ displacement), (start-probe)/scale,
            method="SLSQP", bounds=list(zip((lower-probe)/scale, (upper-probe)/scale)),
            constraints={"type": "eq", "fun": residual}, options={"maxiter": 300, "ftol": 1e-13})
        foot = probe+scale*result.x
        success, message, iterations = bool(result.success), str(result.message), int(result.get("nit", 0))
    remaining = abs(float(function(*foot[None].T)[0]-level))
    normalized = remaining/max(float(np.linalg.norm(gradient(function, foot[None])[0])), 1e-30)
    feasible = bool(np.isfinite(foot).all() and np.isfinite(normalized) and normalized <= 1e-9
        and abs(foot).max() <= .5+1e-12 and np.all(abs(foot[locked]-lower[locked]) <= 1e-12))
    return foot, dict(solver_success=success, message=message, iterations=iterations,
        absolute_field_residual=remaining, normalized_field_residual=normalized,
        feasible_target=feasible, distance=float(np.linalg.norm(foot-probe)))
