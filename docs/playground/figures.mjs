// Build transparent, self-contained plots from the archived FEA measurements.
import { readFile, writeFile } from "node:fs/promises";
import { qualityFigure, bandFigure, pairsFigure } from "../site/assets/guide-figures.mjs";
import { cutFigure, foldFigure, sagFigure } from "../site/assets/course-figures.mjs";

const root = new URL("../../", import.meta.url);
const assets = new URL("docs/site/assets/playground/", root);
const read = async (path) => JSON.parse(await readFile(new URL(path, root), "utf8"));
const style = `<style>
:root{color:#243e40;--teal:#00796d;--copper:#a75027;--blue:#426e9d;--grid:#c4d2d0;--wash:#d6eeea;--figure-edge:#29483f}
@media(prefers-color-scheme:dark){:root{color:#dfece9;--teal:#76d9c1;--copper:#f1ac7e;--blue:#98bbf0;--grid:#43575c;--wash:#284541;--figure-edge:#172f34}}
text{font:14px system-ui,sans-serif;fill:currentColor}.axis{font-size:12px}.grid,.guide-grid{stroke:var(--grid);stroke-dasharray:3 5}.pair-face{fill:var(--wash);stroke:var(--teal)}.pair-link,.pair-point{stroke:var(--copper);fill:var(--copper);stroke-width:2}
</style>`;
const save = (name, body) =>
  writeFile(new URL(name + ".svg", assets), body.replace("</title>", "</title>" + style));
const wrap = (title, w, h, body) =>
  `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" role="img"><title>${title}</title>${body}</svg>`;
const label = (x, y, t, anchor = "start", cls = "") =>
  `<text x="${x}" y="${y}" text-anchor="${anchor}" class="${cls}">${t}</text>`;
const colors = ["var(--teal)", "var(--copper)", "var(--blue)"];
function linePlot(title, series, { xRange, yRange, xTicks, yTicks, xLabel, log = false, reference }) {
  const scale = (n, [a, b]) => (log ? Math.log(n / a) / Math.log(b / a) : (n - a) / (b - a));
  const px = (x) => 66 + scale(x, xRange) * 345,
    py = (y) => 195 - scale(y, yRange) * 133;
  let body = label(12, 22, title);
  for (const y of yTicks)
    body += `<path d="M60 ${py(y)}H423" class="grid"/>` + label(52, py(y) + 4, y, "end", "axis");
  for (const x of xTicks) body += label(px(x), 219, x, "middle", "axis");
  body += label(236, 244, xLabel, "middle", "axis");
  if (reference) {
    const [[x1, y1], [x2, y2]] = reference.ends;
    body +=
      `<path d="M${px(x1)} ${py(y1)} L${px(x2)} ${py(y2)}" stroke="currentColor" stroke-dasharray="5 4" fill="none"/>` +
      label(px(x1) + 2, py(y1) - 7, reference.name, "start", "axis");
  }
  for (const [i, { name, xy }] of series.entries()) {
    const color = colors[i];
    body +=
      `<line x1="${14 + i * 129}" x2="${30 + i * 129}" y1="43" y2="43" stroke="${color}" stroke-width="3"/>` +
      label(35 + i * 129, 47, name, "start", "axis");
    body += `<polyline points="${xy.map(([x, y]) => [px(x), py(y)].join(",")).join(" ")}" fill="none" stroke="${color}" stroke-width="2.6" stroke-dasharray="${i === 2 ? "4 4" : "none"}"/>`;
    for (const [x, y] of xy)
      body += `<circle cx="${px(x)}" cy="${py(y)}" r="3.5" fill="${color}"><title>${name} · ${x.toPrecision(4)}, ${y.toPrecision(4)}</title></circle>`;
  }
  return wrap(title, 446, 258, body);
}
const convergence = await read("docs/solver-convergence-results.json");
for (const [name, key, title, yRange, yTicks, ref] of [
  [
    "displacement",
    "l2_displacement_error",
    "Displacement · L² error",
    [0.00085, 0.018],
    [0.001, 0.002, 0.005, 0.01],
    {
      ends: [
        [0.125, 0.004],
        [0.235, 0.004 * (0.235 / 0.125) ** 2],
      ],
      name: "Slope 2",
    },
  ],
  [
    "energy",
    "energy_error",
    "Elastic energy · error",
    [0.105, 0.66],
    [0.12, 0.2, 0.3, 0.5],
    {
      ends: [
        [0.125, 0.3],
        [0.235, (0.3 * 0.235) / 0.125],
      ],
      name: "Slope 1",
    },
  ],
]) {
  const series = convergence.cases.map((c) => ({
    name: { cube: "Cube", cylinder: "Cylinder", torus: "Torus" }[c.geometry],
    xy: c.refinement.map((r) => [r.maximum_edge_length, r[key]]),
  }));
  await save(
    name,
    linePlot(title, series, {
      xRange: [0.112, 0.245],
      yRange,
      xTicks: [0.12, 0.16, 0.2, 0.24],
      yTicks,
      xLabel: "Maximum physical edge length · log–log",
      log: true,
      reference: ref,
    }),
  );
}
const gyroid = (await read("docs/fedoo-homogenization-results.json")).cases.filter(
  (c) => c.geometry === "gyroid_sheet",
);
await save(
  "stiffness",
  linePlot(
    "Effective axial modulus",
    [{ name: "Gyroid · Eₓ", xy: gyroid.map((c) => [c.cells, c.young_modulus_xyz[0]]) }],
    {
      xRange: [11, 33],
      yRange: [0.378, 0.425],
      xTicks: gyroid.map((c) => c.cells),
      yTicks: [0.38, 0.39, 0.4, 0.41, 0.42],
      xLabel: "Background cells per axis",
    },
  ),
);
const costs = await read("demos/solver/release-wheel/schwarz-matched-cost-20261118/summary.json");
let bars =
  label(16, 26, "Median workflow time · seconds") + label(16, 52, "Native Tet4 / Microgen + repair");
for (const [i, [key, name]] of [
  ["apple-m4", "M4 · macOS"],
  ["core-i7-10750h", "i7-10750H · Linux"],
  ["core-i7-12700h", "i7-12700H · Windows"],
  ["core-ultra-9-285k", "Ultra 9 285K · Windows"],
  ["ryzen-7-7730u", "Ryzen 7 7730U · Linux"],
].entries()) {
  const y = 80 + i * 46;
  bars += label(16, y + 9, name);
  for (const [j, route] of ["near-target", "h050"].entries()) {
    const v = costs[key].bracket.medians[route].wall,
      extent = v * 2.65;
    bars +=
      `<rect x="198" y="${y + j * 13 - 5}" width="${extent}" height="6" rx="3" fill="${colors[j]}"/>` +
      label(205 + extent, y + j * 13 + 2, v.toFixed(1), "start", "axis");
  }
}
await save("matched-cost", wrap("Time at the one percent stiffness threshold", 564, 303, bars));
for (const [name, figure] of [
  ["quality", qualityFigure(5)],
  ["band", bandFigure(30)],
  ["pairs", pairsFigure()],
  ["cut", cutFigure(50)],
  ["sag", sagFigure(60)],
  ["fold", foldFigure(10)],
])
  await save(name, figure.graphic);
