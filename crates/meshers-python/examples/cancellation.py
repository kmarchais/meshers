"""Handle cancellation, invalid options and a resource-budget failure."""

import meshers

token = meshers.CancellationToken()
# A GUI or another Python thread can call cancel() while a mesh is running.
token.cancel()
try:
    meshers.generate("gyroid", cells=16, cancel=token)
except meshers.CancelledError:
    print("Cancelled; use a new token for a new request.")

try:
    meshers.generate("gyroid", cells=0)
except ValueError as error:
    print("Invalid option:", error)

try:
    meshers.generate("gyroid", cells=16, max_tetrahedra=1)
except meshers.MeshingError as error:
    print("Meshing failed:", error)
