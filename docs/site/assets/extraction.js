import { quadStage, triangleStage } from "./extraction-figures.mjs";

const lab = document.querySelector("#surface-extraction");
const steps = [...lab.querySelectorAll('input[name="stage"]')];
const level = lab.querySelector("select");
function draw() {
  const step = steps.findIndex((input) => input.checked);
  for (const [id, result] of [
    ["triangle", triangleStage(Number(level.value), step)],
    ["quad", quadStage(step)],
  ]) {
    const panel = lab.querySelector("#" + id);
    panel.querySelector(".guide-diagram").innerHTML = result.graphic;
    panel.querySelector("figcaption").textContent = result.caption;
  }
}
lab.addEventListener("change", draw);
draw();
lab.dataset.ready = "";
