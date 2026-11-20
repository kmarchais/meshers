# Copyright (c) 2026 Kevin Marchais

"""Keyword handling shared by Python meshing entry points."""

from collections.abc import Mapping
from typing import TypedDict

from ._meshers import CancellationToken


class MeshingOptions(TypedDict, total=False):
    """Callback execution and adaptive edge sizing keywords."""

    batch_size: int
    cancel: CancellationToken | None
    min_length: float | None
    max_length: float | None
    gradation: float


def check_keywords(options: Mapping[str, object], schema: type) -> None:
    """Reject unsupported keywords before field preparation or native meshing.

    Raises:
        TypeError: A keyword does not belong to the requested operation.
    """
    unknown = options.keys() - schema.__annotations__.keys()
    if unknown:
        msg = f"Unexpected keyword argument: {min(unknown)!r}"
        raise TypeError(msg)
