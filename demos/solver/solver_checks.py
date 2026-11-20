# Copyright (c) 2026 Kevin Marchais

"""Numerical acceptance checks shared by the elasticity studies."""

MATRIX_SYMMETRY_ATOL = 1e-10
MAX_RELATIVE_RESIDUAL = 2e-9
AFFINE_DISPLACEMENT_ATOL = 1e-8
AFFINE_ENERGY_ATOL = 1e-7


def require(condition: object, message: str) -> None:
    """Enforce a numerical acceptance check in every Python execution mode.

    Raises:
        ValueError: The numerical check failed.
    """
    if not condition:
        raise ValueError(message)
