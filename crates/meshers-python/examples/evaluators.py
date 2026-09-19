"""Automatic compilation, reusable compiled fields and batched callbacks."""

import meshers
import numpy as np


def sphere(x, y, z):
    return (x - 0.47) ** 2 + (y - 0.51) ** 2 + (z - 0.49) ** 2 - 0.29**2


def gradient(x, y, z):
    return 2 * (x - 0.47), 2 * (y - 0.51), 2 * (z - 0.49)


options = {"cells": 12, "geometry_tolerance": 0.05, "optimize_passes": 1}
automatic = meshers.generate(sphere, **options)
compiled = meshers.compile_field(sphere, gradient=gradient)
print(compiled(0.5, 0.5, 0.5), compiled.gradient(0.5, 0.5, 0.5))
reused = meshers.generate(compiled, **options)


def callback_gradient(x, y, z):
    # Callback gradients have shape (N, 3); compiled gradients return a tuple.
    return np.column_stack(gradient(x, y, z))


callback = meshers.generate(
    sphere, gradient=callback_gradient, compile=False, batch_size=128, **options
)
finite_difference = meshers.generate(sphere, compile=False, batch_size=128, **options)
for result in (automatic, reused, callback, finite_difference):
    print(result.diagnostics["evaluator"], result.diagnostics["volume"])
