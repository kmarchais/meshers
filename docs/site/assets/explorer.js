const lab = document.querySelector(".mesh-lab");
const form = lab.querySelector(".mesh-controls");
const stage = lab.querySelector(".mesh-stage");
const status = lab.querySelector(".mesh-status");
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
let loading = false;
let starter;
let inView = false;
let dirty = true;
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

function show(data) {
  viewer.showMesh(data, lines.checked);
  reset.disabled = false;
  status.textContent = [
    data.tets.toLocaleString() + " tetrahedra", data.nodes.toLocaleString() + " nodes",
    "Min. MMG quality " + data.minQuality.toFixed(3),
    "Min. boundary angle " + data.minAngle.toFixed(1) + "\u00b0",
  ].join(" \u00b7 ");
}

function fail(message) {
  cancelJob();
  status.textContent = message;
}

function generate() {
  cancelJob();
  if (!viewer || !inView || document.hidden) return;
  dirty = false;
  if (field.value === "0" && cells.valueAsNumber === 10 && width.valueAsNumber === 0.5) {
    show(starter);
    return;
  }
  status.textContent = "Generating mesh…";
  try {
    const job = new Worker(new URL("./mesh-worker.js", import.meta.url), { type: "module" });
    worker = job;
    deadline = setTimeout(() => fail("The mesher took too long to load. Try again."), 45_000);
    job.onerror = () => { if (worker === job) fail("The mesher could not run. Try again."); };
    job.onmessage = ({ data }) => {
      if (worker !== job) return;
      if (data.ready) {
        clearTimeout(deadline);
        deadline = setTimeout(() => fail("Meshing timed out. Try fewer cells."), 5_000);
        return;
      }
      cancelJob(); // Drop the worker and its WebAssembly heap after every result.
      if (data.error) { status.textContent = "Mesh rejected. " + data.error; return; }
      show(data);
    };
    job.postMessage({ field: Number(field.value), cells: cells.valueAsNumber, width: width.valueAsNumber });
  } catch {
    fail("The mesher could not start. Try the Python example below.");
  }
}

async function prepare() {
  if (!inView || document.hidden || loading) return;
  if (viewer) {
    if (dirty) generate();
    return;
  }
  const current = session;
  loading = true;
  try {
    const [module, data] = await Promise.all([
      import("./playground/viewer.js"),
      fetch(new URL("./playground/gyroid.json", import.meta.url)).then(response => {
        if (!response.ok) throw new Error("Initial mesh unavailable");
        return response.json();
      }),
    ]);
    if (current !== session || !inView || document.hidden) return;
    starter = { ...data, coords: new Float32Array(data.coords), indices: new Uint32Array(data.indices) };
    viewer = module.createViewer(stage);
    generate();
  } catch {
    if (current === session) fail("The 3D view could not load. You can run the Python example below.");
  } finally {
    if (current === session) loading = false;
  }
}

function updateVisibility() {
  if (inView && !document.hidden) {
    prepare();
  } else if (worker || pending !== undefined) {
    cancelJob();
    dirty = true;
  }
}

form.addEventListener("input", updateParameters);
form.addEventListener("change", ({ target }) => {
  if (target === lines) { viewer?.setEdges(lines.checked); return; }
  dirty = true;
  cancelJob();
  if (inView && !document.hidden) pending = setTimeout(prepare, 200);
});
reset.onclick = () => viewer?.reset();
new IntersectionObserver(([entry]) => {
  inView = entry.isIntersecting;
  updateVisibility();
}).observe(stage);
document.addEventListener("visibilitychange", updateVisibility);
window.addEventListener("pagehide", () => {
  session++;
  cancelJob();
  viewer?.dispose();
  viewer = undefined;
  starter = undefined;
  loading = false;
  dirty = true;
  reset.disabled = true;
});
window.addEventListener("pageshow", updateVisibility);
updateParameters();
