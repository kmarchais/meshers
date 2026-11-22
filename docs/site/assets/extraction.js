import { quadStage, triangleStage } from "./extraction-figures.mjs";

const lab = document.querySelector("#surface-extraction");
function draw() {
  const step = Number(lab.querySelector('input[name="stage"]:checked').value);
  const level = Number(lab.querySelector('input[name="level"]:checked').value);
  for (const [id, result] of [
    ["triangle", triangleStage(level, step)],
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
