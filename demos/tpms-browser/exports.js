export function stlBytes(mesh, scale) {
  const bytes = new ArrayBuffer(84 + 50 * mesh.triangles.length),
    view = new DataView(bytes);
  view.setUint32(80, mesh.triangles.length, true);
  mesh.triangles.forEach((face, index) => {
    const p = face.map((i) => mesh.points[i].map((x) => x * scale));
    const u = p[1].map((x, i) => x - p[0][i]),
      v = p[2].map((x, i) => x - p[0][i]);
    const n = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
      ],
      length = Math.hypot(...n);
    let offset = 84 + index * 50;
    for (const vector of [n.map((x) => x / length), ...p])
      for (const x of vector) {
        view.setFloat32(offset, x, true);
        offset += 4;
      }
  });
  return bytes;
}
export function vtuText(mesh, scale) {
  return `<?xml version="1.0"?><VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian"><UnstructuredGrid><Piece NumberOfPoints="${mesh.points.length}" NumberOfCells="${mesh.tetrahedra.length}"><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">${mesh.points.map((p) => p.map((x) => x * scale).join(" ")).join("\n")}</DataArray></Points><Cells><DataArray type="Int32" Name="connectivity" format="ascii">${mesh.tetrahedra.map((t) => t.join(" ")).join("\n")}</DataArray><DataArray type="Int32" Name="offsets" format="ascii">${mesh.tetrahedra.map((_, i) => 4 * (i + 1)).join(" ")}</DataArray><DataArray type="UInt8" Name="types" format="ascii">${mesh.tetrahedra.map(() => 10).join(" ")}</DataArray></Cells></Piece></UnstructuredGrid></VTKFile>`;
}
