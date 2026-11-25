// Small analytic constructions for the lessons. Each function also renders its static poster.
const line = ([x, y], [u, v], attrs = "") =>
  `<line x1="${x}" y1="${y}" x2="${u}" y2="${v}" ${attrs}/>`;
const point = ([x, y], fill = "currentColor") =>
  `<circle cx="${x}" cy="${y}" r="4" fill="${fill}"/>`;
const label = (x, y, text) => `<text x="${x}" y="${y}" text-anchor="middle">${text}</text>`;
const copper = "var(--chart-copper, var(--copper))";
const mint = "var(--mint, var(--teal))";
const svg = (title, body) =>
  `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 540 192" role="img" aria-label="${title}"><title>${title}</title>${body}</svg>`;

export function cutFigure(percent) {
  // The reference tetrahedron has f=x+y+z. Its edge roots all lie a fraction t from O.
  const t = percent / 100;
  const origin = [76, 124],
    corners = [
      [308, 155],
      [192, 28],
      [412, 75],
    ];
  const roots = corners.map((p) => p.map((v, i) => origin[i] + t * (v - origin[i])));
  const outline = corners.map((p, i) => line(p, corners[(i + 1) % 3]) + line(origin, p)).join("");
  const face = `<polygon points="${roots.map((p) => p.join(",")).join(" ")}" fill="${mint}" fill-opacity=".45" stroke="${mint}" stroke-width="3"/>`;
  return {
    graphic: svg(
      "Three shared edge roots form the cut face",
      `<g fill="none" stroke="currentColor" opacity=".6">${outline}</g>` +
        face +
        roots.map((p) => point(p, copper)).join("") +
        point(origin) +
        label(65, 149, "O") +
        label(285, 181, "Kept volume fraction = " + t.toFixed(2) + "³"),
    ),
    value: t.toFixed(2) + " · " + (100 * t ** 3).toFixed(1) + "% volume",
  };
}

export function sagFigure(percent) {
  const fraction = percent / 100; // chord length as a fraction of the diameter
  const radius = 112,
    cx = 262,
    cy = 147;
  const half = radius * fraction,
    offset = Math.sqrt(1 - fraction ** 2);
  const y = cy - radius * offset;
  return {
    graphic: svg(
      "Sag between a circular arc and its chord",
      `<path d="M${cx - radius} ${cy} A${radius} ${radius} 0 0 1 ${cx + radius} ${cy}" fill="none" stroke="${mint}" stroke-width="3"/>` +
        line([cx - half, y], [cx + half, y], `stroke="${copper}" stroke-width="3"`) +
        line([cx, cy - radius], [cx, y], 'stroke="currentColor" stroke-dasharray="3 2"') +
        point([cx - half, y], copper) +
        point([cx + half, y], copper) +
        label(cx, y + 21, "chord") +
        label(cx, 180, "Radius = 1"),
    ),
    value: "h = " + (2 * fraction).toFixed(2) + " · sag = " + (1 - offset).toFixed(3),
  };
}

export function foldFigure(percent) {
  const d = percent / 100;
  const position = (t) => t + 4 * d * t * (1 - t);
  const axis = (t, x) => [82 + 336 * t, 150 - 108 * x];
  const curve = Array.from({ length: 65 }, (_, i) => axis(i / 64, position(i / 64)).join(",")).join(
    " ",
  );
  const minSlope = 1 - 4 * d;
  return {
    graphic: svg(
      "Physical position along the quadratic reference edge",
      line(axis(0, 0), axis(1, 0), 'stroke="currentColor"') +
        line(axis(0, 0), axis(0, 1.2), 'stroke="currentColor"') +
        line(axis(0, 0), axis(1, 1), 'stroke="currentColor" stroke-dasharray="4 4" opacity=".45"') +
        `<polyline points="${curve}" fill="none" stroke="${mint}" stroke-width="3"/>` +
        [0, 0.5, 1].map((t) => point(axis(t, position(t)), copper)).join("") +
        label(51, 42, "x = 1") +
        label(80, 174, "0") +
        label(250, 174, "Reference position t") +
        label(418, 174, "1"),
    ),
    value:
      "min dx/dt = " +
      minSlope.toFixed(2) +
      (minSlope < 0 ? " · folded" : minSlope === 0 ? " · singular endpoint" : " · increasing"),
  };
}
