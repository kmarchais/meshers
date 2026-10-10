"""Common P1 elasticity test with exact quadratic fields and natural tractions.

The primary score is the worst relative energy error over every linear
combination of 18 quadratic displacement fields. It is a finite test space,
not an assertion about every engineering load. Coordinates are nondimensional.
"""

import time

import numpy as np
from scipy import linalg, sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu

from .geometry import EDGES

LAM = 0.3 / (1.3 * 0.4)
MU = 1 / 2.6
C = np.zeros((6, 6))
C[:3, :3] = LAM
np.fill_diagonal(C[:3, :3], LAM + 2 * MU)
C[3:, 3:] = np.eye(3) * MU
QUADRATURE = np.full((4, 4), 0.1381966011250105)
np.fill_diagonal(QUADRATURE, 0.5854101966249685)


def modes():
    h = np.zeros((24, 3, 3, 3))
    a = np.zeros((24, 3, 3))
    names = []
    for component in range(3):
        for j, k in [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (2, 0)]:
            m = len(names)
            h[m, component, j, k] = 2 if j == k else 1
            h[m, component, k, j] = 2 if j == k else 1
            names.append(f"u{'xyz'[component]}={'xyz'[j]}{'xyz'[k]}")
    for i in range(3):
        a[18 + i, i, i] = 1
    for m, (i, j) in enumerate([(1, 2), (0, 2), (0, 1)], 21):
        a[m, i, j] = a[m, j, i] = 0.5
    return h, a, names + ["patch-xx", "patch-yy", "patch-zz", "patch-yz", "patch-xz", "patch-xy"]


H, A, NAMES = modes()


def exact(x):
    return 0.5 * np.einsum("miab,...a,...b->...mi", H, x, x) + np.einsum("mia,...a->...mi", A, x)


def strain(g):
    return np.stack(
        [
            g[..., 0, 0],
            g[..., 1, 1],
            g[..., 2, 2],
            g[..., 1, 2] + g[..., 2, 1],
            g[..., 0, 2] + g[..., 2, 0],
            g[..., 0, 1] + g[..., 1, 0],
        ],
        axis=-1,
    )


def rigid_pins(x, t):
    edges = np.concatenate([t[:, [i, j]] for i, j in EDGES])
    graph = sparse.coo_matrix(
        (np.ones(len(edges)), (edges[:, 0], edges[:, 1])), shape=(len(x), len(x))
    )
    components, labels = connected_components(graph, directed=False)
    pins = []
    for c in range(components):
        ids = np.flatnonzero(labels == c)
        r = x[ids] - x[ids].mean(axis=0)
        basis = np.zeros((len(ids), 3, 6))
        basis[:, :, :3] = np.eye(3)
        for axis in range(3):
            basis[:, :, 3 + axis] = np.cross(np.eye(3)[axis], r)
        matrix = basis.reshape(-1, 6)
        _, qr, pivot = linalg.qr(matrix.T, pivoting=True, mode="economic")
        if np.min(np.abs(np.diag(qr))) < 1e-12:
            raise ValueError("A component cannot support six independent rigid-motion constraints")
        local = pivot[:6]
        pins.extend((3 * ids[local // 3] + local % 3).tolist())
    return np.asarray(sorted(pins)), components


def worst_error(error, reference):
    """Generalized Rayleigh quotient; invariant to an invertible load-basis change."""
    if np.linalg.eigvalsh(reference)[0] <= 1e-12 * np.linalg.norm(reference, 2):
        raise ValueError("Exact energy matrix is singular or too ill-conditioned")
    eigenvalues, eigenvectors = linalg.eigh((error + error.T) / 2, (reference + reference.T) / 2)
    if eigenvalues[0] < -1e-10:
        raise ValueError("Negative error energy exceeds roundoff tolerance")
    return float(np.sqrt(max(0, eigenvalues[-1]))), eigenvectors[:, -1]


def solve(points, tets, surface, scale):
    start = time.perf_counter()
    x = points / scale
    vertices = x[tets]
    jac = vertices[:, 1:] - vertices[:, :1]
    volumes = np.linalg.det(jac) / 6
    grad = np.zeros((len(tets), 4, 3))
    grad[:, 1:] = np.linalg.inv(jac).transpose(0, 2, 1)
    grad[:, 0] = -grad[:, 1:].sum(axis=1)
    b = np.zeros((len(tets), 6, 12))
    for node in range(4):
        gx, gy, gz = grad[:, node].T
        j = 3 * node
        b[:, 0, j] = gx
        b[:, 1, j + 1] = gy
        b[:, 2, j + 2] = gz
        b[:, 3, j + 1] = gz
        b[:, 3, j + 2] = gy
        b[:, 4, j] = gz
        b[:, 4, j + 2] = gx
        b[:, 5, j] = gy
        b[:, 5, j + 1] = gx
    ke = np.einsum("tki,kl,tlj,t->tij", b, C, b, volumes)
    dofs = (3 * tets[..., None] + np.arange(3)).reshape(-1, 12)
    k = sparse.coo_matrix(
        (
            ke.ravel(),
            (
                np.broadcast_to(dofs[:, :, None], ke.shape).ravel(),
                np.broadcast_to(dofs[:, None, :], ke.shape).ravel(),
            ),
        ),
        shape=(3 * len(x), 3 * len(x)),
    ).tocsc()
    del b, ke
    body = -MU * np.einsum("mijj->mi", H) - (LAM + MU) * np.einsum("mjji->mi", H)
    force = np.zeros((len(x), 24, 3))
    for node in range(4):
        np.add.at(force, tets[:, node], volumes[:, None, None] * body[None] / 4)
    tri = x[surface]
    normal = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    area = np.linalg.norm(normal, axis=1) / 2
    normal /= 2 * area[:, None]
    bary = np.full((3, 3), 1 / 6)
    np.fill_diagonal(bary, 2 / 3)
    for weights in bary:
        sample = np.einsum("a,taj->tj", weights, tri)
        g = np.einsum("miab,tb->tmia", H, sample) + A
        stress = MU * (g + g.transpose(0, 1, 3, 2))
        stress += LAM * np.trace(g, axis1=2, axis2=3)[..., None, None] * np.eye(3)
        traction = np.einsum("tmij,tj->tmi", stress, normal)
        for node in range(3):
            np.add.at(force, surface[:, node], area[:, None, None] * weights[node] * traction / 3)
    f = force.transpose(0, 2, 1).reshape(3 * len(x), 24)
    pins, components = rigid_pins(x, tets)
    free = np.setdiff1d(np.arange(3 * len(x)), pins)
    if not len(free):
        raise ValueError("No unconstrained degrees of freedom")
    u = exact(x).transpose(0, 2, 1).reshape(3 * len(x), 24)
    rhs = f[free] - k[free][:, pins] @ u[pins]
    assembly_s = time.perf_counter() - start
    start = time.perf_counter()
    factor = splu(k[free][:, free])
    u[free] = factor.solve(rhs)
    residual = k[free] @ u - f[free]
    relative_residual = np.linalg.norm(residual, axis=0) / np.maximum(
        np.linalg.norm(rhs, axis=0), 1e-30
    )
    if relative_residual.max() > 1e-8:
        raise ValueError(f"Free equilibrium residual {relative_residual.max():.3g} exceeds 1e-8")
    solve_s = time.perf_counter() - start
    fields = u.reshape(len(x), 3, 24).transpose(0, 2, 1)
    start = time.perf_counter()
    error, reference = np.zeros((24, 24)), np.zeros((24, 24))
    for first in range(0, len(tets), 128):
        sl = slice(first, first + 128)
        discrete = strain(np.einsum("tami,taj->tmij", fields[tets[sl]], grad[sl]))
        for weights in QUADRATURE:
            sample = np.einsum("a,taj->tj", weights, vertices[sl])
            analytical = strain(np.einsum("miab,tb->tmia", H, sample) + A)
            delta = discrete - analytical
            error += np.einsum("tmi,ij,tnj,t->mn", delta, C, delta, volumes[sl] / 4)
            reference += np.einsum("tmi,ij,tnj,t->mn", analytical, C, analytical, volumes[sl] / 4)
    relative = np.sqrt(np.maximum(0, error.diagonal()) / reference.diagonal())
    if relative[18:].max() > 1e-8:
        raise ValueError(f"Affine patch energy error {relative[18:].max():.3g} exceeds 1e-8")
    worst, combination = worst_error(error[:18, :18], reference[:18, :18])
    metrics = {
        "worst_energy_error": worst,
        "per_load_errors": dict(zip(NAMES[:18], relative[:18].tolist(), strict=True)),
        "patch_max_error": float(relative[18:].max()),
        "free_residual": float(relative_residual.max()),
        "components": components,
        "dofs": 3 * len(x),
        "free_dofs": len(free),
        "assembly_s": assembly_s,
        "solve_s": solve_s,
        "evaluation_s": time.perf_counter() - start,
        "worst_combination": combination.tolist(),
    }
    return metrics, {
        "displacements": fields,
        "error_energy": error,
        "exact_energy": reference,
        "pins": pins,
    }
