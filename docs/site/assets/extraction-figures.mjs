// Analytic level sets used only in the teaching diagrams.
const surface = "var(--mint, var(--teal))";
const root = "var(--chart-copper, var(--copper))";
const stroke = (p, q, color = "currentColor", extra = "") =>
  `<path d="M${p} L${q}" fill="none" stroke="${color}" ${extra}/>`;
const dot = (p, color, r = 4) => `<circle cx="${p[0]}" cy="${p[1]}" r="${r}" fill="${color}"/>`;
const face = (p, color, opacity) =>
  `<polygon points="${p.join(" ")}" fill="${color}" fill-opacity="${opacity}" stroke="${color}" stroke-width="2"/>`;
const text = ([x, y], value) => `<text x="${x}" y="${y}" text-anchor="middle">${value}</text>`;
const svg = (title, body) =>
  `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 330 200" role="img" aria-label="${title}"><title>${title}</title>${body}</svg>`;

const tetra = [
  [0, 0, 0],
  [1, 0, 0],
  [1, 1, 0],
  [1, 1, 1],
];
const edges = [
  [0, 1],
  [0, 2],
  [0, 3],
  [1, 2],
  [1, 3],
  [2, 3],
];
const projectTet = ([x, y, z]) => [44 + 236 * x - 169 * y + 45 * z, 122 + 7 * x + 41 * y - 144 * z];
const mix = (a, b, t) => a.map((v, i) => v + t * (b[i] - v));

export function triangleSection(level) {
  // f=x+y+z takes values 0,1,2,3 at the four vertices, in that order.
  const crossings = edges
    .filter(([a, b]) => a < level && level < b)
    .map(([a, b]) => mix(tetra[a], tetra[b], (level - a) / (b - a)));
  if (crossings.length === 4) [crossings[2], crossings[3]] = [crossings[3], crossings[2]];
  // An exact corner belongs to several edges but supplies just one section vertex.
  if (Number.isInteger(level) && tetra[level]) crossings.push(tetra[level]);
  return crossings;
}

export function triangleStage(level, step) {
  const section = triangleSection(level).map(projectTet);
  const corners = tetra.map(projectTet);
  let g = edges
    .map(([a, b]) => stroke(corners[a], corners[b], "currentColor", 'opacity=".45"'))
    .join("");
  if (step >= 2 && section.length >= 3) g += face(section, surface, step === 3 ? ".44" : ".08");
  // For this tetrahedron and plane family the AD-to-BC diagonal gives the
  // better minimum triangle shape score throughout the four-root interval.
  if (step === 3 && section.length === 4)
    g += stroke(section[1], section[3], surface, 'stroke-width="2"');
  corners.forEach((p, i) => {
    g +=
      dot(p, "currentColor", 3) +
      text(
        [p[0], p[1] + (i === 3 ? -7 : 18)],
        String.fromCharCode(65 + i) +
          " " +
          (step === 0 ? (i - level).toFixed(2) : i === level ? "0" : i < level ? "−" : "+"),
      );
  });
  if (step >= 1) g += section.map((p) => dot(p, root)).join("");
  const n = section.length;
  const captions = [
    `${tetra.filter((_, i) => i < level).length} of 4 corners lie below the plane. Labels show f − t.`,
    Number.isInteger(level) && n
      ? `The plane passes through corner ${"ABCD"[level]}. Count shared endpoints only once.`
      : `${n} of 6 edges have opposite signs at their ends. Each carries one root.`,
    n === 1
      ? "One corner touches the plane. There is no polygon to join."
      : n
        ? `Join the ${n} roots in order around the cutting plane.`
        : "No crossing edges. There is no contour to join.",
    n === 4
      ? "Four roots form a quadrilateral. One diagonal gives two triangles."
      : n === 3
        ? "Three roots form one triangle. Its vertices lie on the original edges."
        : n === 1
          ? "The plane touches one corner. This tetrahedron emits no face."
          : "All four corners lie below the plane. This tetrahedron emits no face.",
  ];
  return {
    graphic: svg("One Freudenthal tetrahedron cut by x + y + z = " + level, g),
    caption: captions[step],
  };
}

// Four unit cubes around the vertical edge at x=y=1. The tilted plane
// z=.25+.12x+.08y crosses all nine vertical edges, with one loop per cube.
const grid = Array.from({ length: 9 }, (_, i) => [i % 3, Math.floor(i / 3)]);
const cells = [
  [0, 1, 4, 3],
  [1, 2, 5, 4],
  [4, 5, 8, 7],
  [3, 4, 7, 6],
];
const projectGrid = ([x, y, z]) => [166 + 57 * x - 58 * y, 72 + 28 * x + 25 * y - 56 * z];
const roots = grid.map(([x, y]) => [x, y, 0.25 + 0.12 * x + 0.08 * y]);
const centers = cells.map((ids) =>
  [0, 1, 2].map((k) => ids.reduce((s, i) => s + roots[i][k], 0) / 4),
);

export function quadStage(step) {
  let g = "";
  for (const z of [0, 1]) {
    for (const i of [0, 1, 2]) {
      g += stroke(projectGrid([0, i, z]), projectGrid([2, i, z]), "currentColor", 'opacity=".25"');
      g += stroke(projectGrid([i, 0, z]), projectGrid([i, 2, z]), "currentColor", 'opacity=".25"');
    }
  }
  grid.forEach(([x, y], i) => {
    g += stroke(
      projectGrid([x, y, 0]),
      projectGrid([x, y, 1]),
      i === 4 ? surface : "currentColor",
      i === 4 ? 'stroke-width="3"' : 'opacity=".25"',
    );
  });
  if (step === 0) {
    for (const [z, sign] of [
      [0, "−"],
      [1, "+"],
    ]) {
      const [x, y] = projectGrid([1, 1, z]);
      g += text([x + 12, y], sign);
    }
  }
  if (step >= 2) {
    for (const ids of cells)
      g += face(
        ids.map((i) => projectGrid(roots[i])),
        root,
        ".035",
      );
  }
  if (step === 3) g += face(centers.map(projectGrid), surface, ".65");
  if (step >= 1) g += roots.map((p) => dot(projectGrid(p), root, 3)).join("");
  if (step >= 2) g += centers.map((p) => dot(projectGrid(p), surface, 5)).join("");
  const captions = [
    "Four cubes meet at the highlighted edge. Its bottom is inside (−) and its top is outside (+).",
    "The tilted plane cuts nine vertical edges. The central crossing belongs to all four cubes.",
    "Each cube has one loop of four roots. Average the four positions to place one teal vertex per loop.",
    "Walk around the central edge. Its four neighbouring vertices form one quad.",
  ];
  return { graphic: svg("Four cubes around a common vertical edge", g), caption: captions[step] };
}
