const lab = document.querySelector(".mesh-lab");
const form = lab.querySelector(".mesh-controls");
const stage = lab.querySelector(".mesh-stage");
const status = lab.querySelector(".mesh-status");
const start = form.querySelector("#mesh-start");
const reset = form.querySelector("#mesh-reset");
const lines = form.querySelector("#mesh-lines");
const cells = form.querySelector("#mesh-cells");
const width = form.querySelector("#mesh-width");
const field = form.elements.namedItem("field");
const code = document.querySelector("#mesh-code");

function updateParameters() {
  for (const slider of [cells, width]) {
    const fraction = (slider.valueAsNumber - Number(slider.min)) / (Number(slider.max) - Number(slider.min));
    slider.style.setProperty("--fill", fraction * 100 + "%");
  }
  form.querySelector('output[for="mesh-cells"]').value = cells.value;
  form.querySelector('output[for="mesh-width"]').value = width.valueAsNumber.toFixed(2);
  code.querySelector('[data-param="field"]').textContent = field.value === "0" ? '"gyroid"' : '"schwarz_p"';
  code.querySelector('[data-param="cells"]').textContent = cells.value;
  for (const span of code.querySelectorAll('[data-param="width"]')) span.textContent = width.value;
}

let viewer;
let active = false;
let session = 0;
let worker;
let deadline;
let pending;

function cancelJob() {
  clearTimeout(pending);
  pending = undefined;
  clearTimeout(deadline);
  worker?.terminate();
  worker = undefined;
}

function stop(message = "Explorer stopped. Start again to generate a mesh.") {
  active = false;
  session++;
  cancelJob();
  viewer?.dispose();
  viewer = undefined;
  reset.disabled = true;
  start.textContent = "Start 3D explorer";
  status.textContent = message;
}

function generate() {
  cancelJob();
  if (!active || !viewer || document.hidden) return;
  status.textContent = "Generating mesh…";
  try {
    const job = new Worker(new URL("./mesh-worker.js", import.meta.url));
    worker = job;
    deadline = setTimeout(() => stop("The mesher took too long to load. Try again."), 45_000);
    job.onerror = () => { if (worker === job) stop("The mesher could not run. Try again."); };
    job.onmessage = ({ data }) => {
      if (worker !== job) return;
      if (data.ready) {
        clearTimeout(deadline);
        deadline = setTimeout(() => stop("Meshing timed out. Try fewer cells."), 5_000);
        return;
      }
      cancelJob(); // Drop the worker and its WebAssembly heap after every result.
      if (data.error) { status.textContent = "Mesh rejected. " + data.error; return; }
      viewer.showMesh(data, lines.checked);
      reset.disabled = false;
      status.textContent = [
        data.tets.toLocaleString() + " tetrahedra", data.nodes.toLocaleString() + " nodes",
        "Min. MMG quality " + data.minQuality.toFixed(3),
        "Min. boundary angle " + data.minAngle.toFixed(1) + "\u00b0",
      ].join(" \u00b7 ");
    };
    job.postMessage({ field: Number(field.value), cells: cells.valueAsNumber, width: width.valueAsNumber });
  } catch {
    stop("The mesher could not start. Try the Python example below.");
  }
}

form.addEventListener("input", updateParameters);
form.addEventListener("change", ({ target }) => {
  if (target === lines) { viewer?.setEdges(lines.checked); return; }
  if (!viewer) return;
  cancelJob();
  status.textContent = "Generating mesh…";
  pending = setTimeout(generate, 200);
});
reset.onclick = () => viewer?.reset();
start.onclick = async () => {
  if (active) { stop(); return; }
  active = true;
  const current = ++session;
  start.textContent = "Stop explorer";
  status.textContent = "Starting 3D view…";
  try {
    const { createViewer } = await import("./playground/viewer.js");
    if (session !== current) return;
    viewer = createViewer(stage);
    generate();
  } catch {
    stop("The 3D view needs WebGL 2. You can run the Python example below.");
  }
};
document.addEventListener("visibilitychange", () => {
  if (document.hidden && active && (!viewer || worker || status.textContent === "Generating mesh…")) {
    stop("Meshing paused while the page was hidden. Start again to resume.");
  }
});
window.addEventListener("pagehide", () => { if (active) stop(); });
updateParameters();
start.disabled = false;
