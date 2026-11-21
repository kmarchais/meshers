const ready = fetch("playground/meshers.wasm")
  .then(response => {
    if (!response.ok) throw new Error("Could not load the mesher.");
    return response.arrayBuffer();
  })
  .then(bytes => WebAssembly.instantiate(bytes))
  .then(({ instance }) => instance.exports);

self.onmessage = async ({ data: { field, cells, width, id } }) => {
  try {
    const core = await ready;
    const address = core.generate(field, cells, width);
    const bytes = new Uint8Array(core.memory.buffer, address, core.response_length());
    const mesh = JSON.parse(new TextDecoder().decode(bytes));
    if (mesh.error) throw new Error(mesh.error);
    // Share wall normals between adjacent triangles and keep the six cap planes sharp.
    const vertices = new Map();
    const positions = [];
    const wallIndices = [], capIndices = [];
    mesh.faces.forEach((face, i) => {
      const tag = mesh.tags[i];
      const group = tag === 0 ? wallIndices : capIndices;
      for (const node of face) {
        const key = node * 7 + tag;
        if (!vertices.has(key)) {
          vertices.set(key, positions.length / 3);
          positions.push(...mesh.points[node]);
        }
        group.push(vertices.get(key));
      }
    });
    const coords = new Float32Array(positions);
    const indices = new Uint32Array([...wallIndices, ...capIndices]);
    self.postMessage({
      id, coords, indices, wallCount: wallIndices.length,
      tets: mesh.tets, nodes: mesh.points.length,
    }, [coords.buffer, indices.buffer]);
  } catch (error) {
    self.postMessage({ id, error: error.message });
  }
};
