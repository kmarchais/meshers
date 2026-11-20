# Copyright (c) 2026 Kevin Marchais

"""Input fields consumed from the frozen release-wheel baseline."""

from typing import NotRequired, TypedDict


class NativeConfig(TypedDict):
    """Native settings read by the periodic mesh generator."""

    cells: int
    adapt: bool
    tolerance: float
    min_length: float
    order: NotRequired[int]
    phase_x: NotRequired[float]
    optimization_passes: NotRequired[int]
    max_length: NotRequired[float]
    gradation: NotRequired[float]
    adapt_tolerance: NotRequired[float]


class Limits(TypedDict):
    """Numerical gates shared by generation and the independent audit."""

    minimum_quality: float
    max_nodes_or_elements: int
    relative_volume_error: float
    two_sided_sampled_distance: float


class Route(TypedDict):
    """Settings and optional archived stiffness for one route."""

    config: NativeConfig
    archived_stiffness: NotRequired[list[list[float]]]


class Baseline(TypedDict):
    """Portion of the release record consumed by the solver studies."""

    limits: Limits
    routes: dict[str, Route]
    reference_volumes: dict[str, float]
