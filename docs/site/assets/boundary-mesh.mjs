export function boundaryMesh(mesh) {
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
  return {
    coords, indices, wallCount: wallIndices.length, tets: mesh.tets, nodes: mesh.points.length,
    minQuality: mesh.minQuality, minAngle: mesh.minAngle,
  };
}
