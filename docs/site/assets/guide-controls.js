import { bandFigure, pairsFigure, qualityFigure } from "./guide-figures.mjs";

const figures = { band: bandFigure, pairs: pairsFigure, quality: qualityFigure };
for (const figure of document.querySelectorAll("[data-guide]")) {
  const input = figure.querySelector("input");
  const output = figure.querySelector("output");
  const diagram = figure.querySelector(".guide-diagram");
  function draw() {
    const { graphic, value } = figures[figure.dataset.guide](input.valueAsNumber);
    diagram.innerHTML = graphic;
    output.value = value;
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
if (section)
  section.addEventListener("change", (event) => {
    const isolated = event.target.value === "poor";
    const image = section.querySelector("img");
    image.src = new URL(
      isolated ? "gyroid-quality-poor.webp" : "gyroid-quality-section.webp",
      import.meta.url,
    );
    image.alt = isolated
      ? "The 46 cells below MMG quality 0.1, in their original positions"
      : "Cutaway through the gyroid colored by tetrahedron quality";
  });
