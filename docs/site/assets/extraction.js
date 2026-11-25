import { quadStage, triangleStage } from "./extraction-figures.mjs";
import { fillSlider } from "./slider.js";

const lab = document.querySelector("#surface-extraction");
const level = lab.querySelector("#cut-plane");
const stage = () => Number(lab.querySelector('input[name="stage"]:checked').value);
function show(id, result) {
  const panel = lab.querySelector(id);
  panel.querySelector(".guide-diagram").innerHTML = result.graphic;
  const caption = panel.querySelector("figcaption");
  if (caption.textContent !== result.caption) caption.textContent = result.caption;
}
function movePlane() {
  lab.querySelector("output").value = level.valueAsNumber.toFixed(2);
  fillSlider(level);
  show("#triangle", triangleStage(level.valueAsNumber, stage()));
}
function changeStage() {
  movePlane();
  show("#quad", quadStage(stage()));
}
level.addEventListener("input", movePlane);
lab.querySelector(".mesh-fields").addEventListener("change", changeStage);
changeStage();
lab.dataset.ready = "";
