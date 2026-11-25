// Orthographic teaching figures. All coordinates and scores use the displayed geometry.
const project = ([x, y, z], cx, cy, scale) => [
  cx + scale * (0.86 * x - 0.5 * y),
  cy + scale * (0.22 * x + 0.38 * y - 0.9 * z),
];
const points = (vertices) => vertices.map((p) => p.map((x) => x.toFixed(2)).join(",")).join(" ");
const text = (x, y, value) => `<text x="${x}" y="${y}" text-anchor="middle">${value}</text>`;
const line = (a, b, style = "") =>
  `<line x1="${a[0]}" y1="${a[1]}" x2="${b[0]}" y2="${b[1]}" ${style}/>`;
const svg = (label, content) =>
  `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 580 190" role="img" aria-label="${label}"><title>${label}</title>${content}</svg>`;

export function qualityFigure(percent) {
  const h = percent / 100;
  // Equilateral unit base, apex above its centroid at h times the regular height.
  // Meshers MMG = (12 (|det|/2)^(2/3) / sum(edge²))^(3/2).
  const score = h * (6 / (4 + 2 * h * h)) ** 1.5;
  function cell(cx, height, q) {
    const v = [
      [-0.5, -0.29, 0],
      [0.5, -0.29, 0],
      [0, Math.sqrt(3) / 2 - 0.29, 0],
      [0, Math.sqrt(3) / 6 - 0.29, Math.sqrt(2 / 3) * height],
    ];
    const p = v.map((a) => project(a, cx, 126, 138));
    const color = [0, 1, 2].map((i) => Math.round([199, 100, 60][i] * (1 - q) + [45, 159, 143][i] * q));
    // Back-to-front by centroid depth along the viewing direction.
    return [
      [0, 1, 2],
      [0, 1, 3],
      [1, 2, 3],
      [2, 0, 3],
    ]
      .map((face, i) => ({
        depth: face.reduce((s, j) => s + 0.45 * v[j][0] + 0.78 * v[j][1] + 0.44 * v[j][2], 0),
        face,
        shade: [0.73, 0.87, 1, 0.78][i],
      }))
      .sort((a, b) => a.depth - b.depth)
      .map(
        ({ face, shade }) =>
          `<polygon points="${points(face.map((j) => p[j]))}" fill="rgb(${color.map((c) => Math.round(c * shade)).join(" ")})" stroke="var(--figure-edge)" stroke-width="1.4"/>`,
      )
      .join("");
  }
  return {
    graphic: svg(
      "Regular and adjustable tetrahedra colored by tetrahedron shape quality",
      cell(145, 1, 1) +
        cell(435, h, score) +
        text(145, 176, "Regular · 1.000") +
        text(435, 176, "Adjusted · " + score.toFixed(3)),
    ),
    value: percent + "% of regular height",
  };
}

export function pairsFigure() {
  const base = [
    [0, 0],
    [1, 0],
    [1, 1],
    [0, 1],
  ];
  const nodes = [
    [0.32, 0.34],
    [0.68, 0.68],
    [0.68, 0.22],
  ];
  const at = (x, [y, z]) => project([x * 2, y, z], 152, 121, 110);
  let body = base.map((p) => line(at(0, p), at(1, p), 'class="guide-grid"')).join("");
  for (const x of [0, 1]) {
    const outline = base.map((p) => at(x, p));
    body += `<polygon points="${points(outline)}" class="pair-face"/>`;
  }
  for (const node of nodes) {
    const a = at(0, node), b = at(1, node);
    body += line(a, b, 'class="pair-link" stroke-dasharray="5 5"');
    for (const [cx, cy] of [a, b]) body += `<circle cx="${cx}" cy="${cy}" r="4" class="pair-point"/>`;
  }
  return {
    graphic: svg(
      "Corresponding nodes on opposite x caps",
      body + text(127, 172, "x low · tag 1") + text(365, 172, "x high · tag 2"),
    ),
  };
}

export function bandFigure(percent) {
  const width = percent / 100,
    px = (x) => 62 + 221 * x,
    py = (y) => 88 - 140 * y;
  const curve = (sign) =>
    Array.from({ length: 61 }, (_, i) => [px(i / 30), py((sign * width) / (1 + i / 30))]);
  const body =
    `<polygon points="${points([...curve(1), ...curve(-1).reverse()])}" class="pair-face"/>` +
    line([px(0), py(0)], [px(2), py(0)], 'class="guide-grid"') +
    text(px(0) + 30, 176, "x = 0") +
    text(px(2) - 30, 176, "x = 2");
  return {
    graphic: svg("Section of the band around (1+x)y = 0", body),
    value:
      "±" +
      width.toFixed(2) +
      " · thickness " +
      (2 * width).toFixed(2) +
      " → " +
      ((2 * width) / 3).toFixed(2),
  };
}
