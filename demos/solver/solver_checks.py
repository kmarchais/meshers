"""Numerical acceptance checks shared by the elasticity studies."""


def require(condition: object, message: str) -> None:
    """Enforce a numerical acceptance check in every Python execution mode.

    Raises:
        ValueError: The numerical check failed.
    """
    if not condition:
        raise ValueError(message)
