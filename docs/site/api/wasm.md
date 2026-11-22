---
icon: simple/webassembly
---

# WebAssembly for the mesh explorer

The [mesh explorer](../index.md) runs the Rust volume mesher in a browser worker. Its WebAssembly interface ships with the documentation and can change with the explorer. Use the [Python](python.md) or [Rust](rust.md) API for custom fields and volume connectivity.

## Generate a periodic sheet

Instantiate the [explorer module](../assets/playground/meshers.wasm) in a worker. Pass its `instance.exports` to this function.

```javascript
function generateGyroid(core) {
  const ptr = core.generate(0, 10, 0.5);
  const bytes = new Uint8Array(core.memory.buffer, ptr, core.response_length());
  const mesh = JSON.parse(new TextDecoder().decode(bytes));
  if (mesh.error) throw new Error(mesh.error);
  return mesh;
}
```

| Argument | Accepted values |
| --- | --- |
| `field` | `0` for gyroid, `1` for Schwarz P |
| `cells` | Integer from 6 to 16 on each axis |
| `width` | Band half-width from 0.3 to 0.8 in field units |

All axes are periodic. The explorer uses a geometry tolerance of 0.05 and two optimization passes.

## Read the result

| JSON fields | Meaning |
| --- | --- |
| `points`, `faces`, `tags` | 3D node coordinates, zero-based boundary triangles, and one boundary tag per triangle |
| `tets` | Number of volume tetrahedra |
| `minQuality`, `minAngle` | Minimum [tetrahedron shape quality](../quality.md) and smallest boundary triangle angle in degrees |
| `error` | On rejection, the response contains an error message in place of mesh data |

The module owns the response buffer. Decode it before the next call to `generate` reuses it. The explorer terminates the worker after each result and enforces a five-second generation timeout. The module caps linear memory at 32 MiB.
