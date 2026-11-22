import { fillSlider } from "./slider.js";
import { bandFigure, qualityFigure } from "./guide-figures.mjs";
import { cutFigure, foldFigure, sagFigure } from "./course-figures.mjs";

const figures = { band: bandFigure, quality: qualityFigure, cut: cutFigure, sag: sagFigure, fold: foldFigure };
for (const figure of document.querySelectorAll("[data-guide]")) {
  const input = figure.querySelector("input");
  const output = figure.querySelector("output");
  const diagram = figure.querySelector(".guide-diagram");
  function draw() {
    const { graphic, value } = figures[figure.dataset.guide](input.valueAsNumber);
    diagram.innerHTML = graphic;
    output.value = value;
    fillSlider(input);
  }
  input.addEventListener("input", draw);
  draw();
  figure.dataset.ready = "";
}

const layers = document.querySelector("#cut-layers");
if (layers)
  layers.addEventListener("change", (event) => {
    const image = layers.querySelector("img");
    const url = new URL(image.src);
    url.hash = event.target.value;
    image.src = url.href;
  });

const section = document.querySelector("#quality-section");
if (section) {
  import("./shape-range.js")
    .then(({ mountShapeRange }) => mountShapeRange(section))
    .catch(() => {
      section.querySelector('[role="status"]').textContent =
        "The interactive filter is unavailable. Showing the 46 cells up to 0.10.";
    });
}
