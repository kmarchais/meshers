import * as THREE from "three";
import { stlBytes, vtuText } from "./exports.js";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

const $ = (id) => document.getElementById(id);
const names = ["Gyroid", "Split-P", "Schwarz P"],
  presets = ["Quick view", "3D printing", "FEA"],
  defaults = [12, 20, 16];
let preset = 0,
  worker = null,
  current = null,
  history = [],
  busy = false,
  queue = [],
  timer = null,
  start = 0;
let wire = false,
  quality = false,
  tetWire = false;
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(38, 1, 0.01, 100);
camera.position.set(1.9, 1.45, 2.1);
const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
renderer.localClippingEnabled = true;
renderer.setClearColor(0x000000, 0);
renderer.outputColorSpace = THREE.SRGBColorSpace;
$("viewport").append(renderer.domElement);
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
scene.add(new THREE.HemisphereLight(0xebf6dc, 0x263831, 2.5));
const light = new THREE.DirectionalLight(0xffe4bc, 3);
light.position.set(2, 4, 3);
scene.add(light);
const rim = new THREE.DirectionalLight(0x8cbcb2, 2);
rim.position.set(-3, 1, -2);
scene.add(rim);
const group = new THREE.Group();
scene.add(group);
const plane = new THREE.Plane(new THREE.Vector3(-1, 0, 0), 1);
const material = new THREE.MeshStandardMaterial({
  color: 0xc3dba6,
  roughness: 0.6,
  metalness: 0.12,
  side: THREE.DoubleSide,
  clippingPlanes: [plane],
});
let mesh = null,
  wireMesh = null,
  tetMesh = null,
  box = null;
const fmt = (n) => n.toLocaleString(undefined, { maximumFractionDigits: 0 });
const num = (n, d = 3) => Number(n).toFixed(d);
const small = (n) => (n < 0.001 ? n.toExponential(2) : num(n));
new ResizeObserver(() => {
  const { width, height } = $("viewport").getBoundingClientRect();
  renderer.setSize(width, height, false);
  camera.aspect = width / height;
  // Keep the whole cell visible in the narrow in-app browser panel.
  camera.zoom = Math.min(1, camera.aspect);
  camera.updateProjectionMatrix();
}).observe($("viewport"));
function frame() {
  requestAnimationFrame(frame);
  controls.update();
  renderer.render(scene, camera);
}
frame();

function status(text, error = false, running = false) {
  $("status").hidden = !text;
  $("status").classList.toggle("error", error);
  $("status-text").textContent = text;
  $("stop").hidden = !running;
  if (!running) $("elapsed").textContent = "";
}
function config() {
  return {
    shape: +$("shape").value,
    preset,
    repeat: +$("repeat").value,
    resolution: +$("resolution").value,
    thickness: +$("thickness").value,
    grade: +$("grade").value,
    cellSize: +$("cell-size").value,
  };
}
function labels() {
  const c = config();
  $("repeat-value").textContent = `${c.repeat} × ${c.repeat} × ${c.repeat}`;
  $("thickness-value").textContent = num(c.thickness, 2);
  $("grade-value").textContent = c.grade
    ? `±${Math.round(100 * c.grade)}%`
    : "Uniform";
  $("resolution-value").textContent = c.resolution;
  $("cell-size-value").textContent = `${c.cellSize} mm`;
  $("budget").textContent =
    `${c.repeat * c.resolution - 1} grid divisions per domain axis.`;
}
function choosePreset(p, reset = true) {
  preset = p;
  document
    .querySelectorAll(".preset")
    .forEach((b) => b.classList.toggle("active", +b.dataset.preset === p));
  if (reset)
    $("resolution").value = Math.min(
      defaults[p],
      Math.floor(((p === 2 ? 40 : 64) + 1) / +$("repeat").value),
    );
  labels();
}
function dirty() {
  labels();
  if (current)
    $("stale").textContent =
      "Settings changed. Generate to update the displayed mesh.";
}
document.querySelectorAll(".preset").forEach(
  (b) =>
    (b.onclick = () => {
      choosePreset(+b.dataset.preset);
      dirty();
    }),
);
for (const id of [
  "shape",
  "repeat",
  "resolution",
  "thickness",
  "grade",
  "cell-size",
])
  $(id).addEventListener("input", dirty);
function newWorker() {
  const url = URL.createObjectURL(
    new Blob([WORKER_SOURCE], { type: "text/javascript" }),
  );
  const w = new Worker(url);
  URL.revokeObjectURL(url);
  w.onerror = (e) => {
    e.preventDefault();
    finishError(
      e.message || "The meshing worker stopped. Reduce sampling and try again.",
    );
  };
  w.onmessage = ({ data }) => {
    if (!data.ok) {
      finishError(data.error);
      return;
    }
    clearInterval(timer);
    data.id = crypto.randomUUID ? crypto.randomUUID() : String(Date.now());
    history.unshift(data);
    history = history.slice(0, 6);
    show(data);
    renderHistory();
    if (queue.length) {
      launch(queue.shift());
    } else {
      setBusy(false);
      status("");
    }
  };
  return w;
}
function setBusy(value) {
  busy = value;
  $("generate").disabled = value;
  $("compare").disabled = value;
  $("clear").disabled = value;
  document
    .querySelectorAll("aside input,aside select,.preset,#history button")
    .forEach((e) => (e.disabled = value));
}
function launch(c) {
  setBusy(true);
  worker ??= newWorker();
  start = performance.now();
  status(`Generating ${presets[c.preset].toLowerCase()} mesh… `, false, true);
  timer = setInterval(
    () =>
      ($("elapsed").textContent =
        `${((performance.now() - start) / 1000).toFixed(1)} s`),
    100,
  );
  worker.postMessage(c);
}
function finishError(message) {
  clearInterval(timer);
  queue = [];
  worker?.terminate();
  worker = null;
  setBusy(false);
  status(`Generation failed: ${message}`, true);
}
$("generate").onclick = () => {
  queue = [];
  launch(config());
};
$("compare").onclick = () => {
  const c = config();
  queue = [0, 1, 2].map((p) => ({
    ...c,
    preset: p,
    resolution: Math.min(
      defaults[p],
      Math.floor(((p === 2 ? 40 : 64) + 1) / c.repeat),
    ),
  }));
  launch(queue.shift());
};
$("stop").onclick = () => {
  clearInterval(timer);
  worker?.terminate();
  worker = null;
  queue = [];
  setBusy(false);
  status("Stopped. Your previous mesh is still available.");
};

function disposeView() {
  for (const obj of [...group.children]) {
    group.remove(obj);
    obj.geometry?.dispose();
    if (obj.material !== material) obj.material?.dispose();
  }
  mesh = wireMesh = tetMesh = box = null;
}
function colorFor(angle) {
  return new THREE.Color().setHSL(
    (Math.min(angle, 60) / 60) * 0.26 + 0.025,
    0.49,
    0.58,
  );
}
function show(result) {
  current = result;
  const { mesh: m, config: c, seconds } = result;
  $("shape").value = c.shape;
  $("repeat").value = c.repeat;
  $("resolution").value = c.resolution;
  $("thickness").value = c.thickness;
  $("grade").value = c.grade;
  $("cell-size").value = c.cellSize;
  choosePreset(c.preset, false);
  $("stale").textContent = "";
  $("scene-mode").textContent =
    `${presets[c.preset]} / ${c.preset === 2 ? "tetrahedral volume" : "surface mesh"}`;
  $("scene-title").textContent =
    `${names[c.shape]} · ${c.repeat === 1 ? "1 cell" : c.repeat + "³ cells"}${c.grade ? " · graded" : ""}`;
  disposeView();
  const positions = new Float32Array(m.triangles.length * 9),
    colors = new Float32Array(positions.length);
  m.triangles.forEach((f, i) => {
    const color = colorFor(m.angles[i]);
    f.forEach((v, j) => {
      positions.set(
        m.points[v].map((x) => x / c.repeat),
        i * 9 + j * 3,
      );
      colors.set([color.r, color.g, color.b], i * 9 + j * 3);
    });
  });
  const geo = new THREE.BufferGeometry();
  geo.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  geo.setAttribute("color", new THREE.BufferAttribute(colors, 3));
  geo.computeVertexNormals();
  mesh = new THREE.Mesh(geo, material);
  group.add(mesh);
  box = new THREE.Box3Helper(
    new THREE.Box3(
      new THREE.Vector3(-0.5, -0.5, -0.5),
      new THREE.Vector3(0.5, 0.5, 0.5),
    ),
    0x52715c,
  );
  group.add(box);
  $("volume-wire").disabled = !m.tetrahedra.length;
  tetWire = false;
  $("volume-wire").classList.remove("on");
  $("slice").value = 100;
  plane.constant = 0.51;
  updateView();
  $("time").textContent =
    seconds < 1 ? `${Math.round(seconds * 1000)} ms` : `${num(seconds, 2)} s`;
  $("triangles").textContent = fmt(m.triangles.length);
  $("tets").textContent = m.tetrahedra.length
    ? fmt(m.tetrahedra.length)
    : "Surface only";
  const q = m.metrics;
  $("closed").textContent =
    q.bad_edges === 0 ? "Pass" : `${q.bad_edges} bad edges`;
  $("periodic").textContent = q.periodic_matches
    .map((x) => (x === null ? "—" : x ? "✓" : "✕"))
    .join(" / ");
  $("angle").textContent = `${num(q.angle_p01, 1)}°`;
  $("area").textContent = `${num(q.area_cv, 2)} CV`;
  $("deviation").textContent =
    `${small(q.centroid_distance_p95 * c.cellSize)} mm`;
  $("tet-quality").textContent = q.volume ? num(q.volume.minimum_quality) : "—";
  const target = q.volume
    ? q.volume.target_met &&
      q.bad_edges === 0 &&
      !q.periodic_matches.includes(false)
    : q.bad_edges === 0 && !q.periodic_matches.includes(false);
  $("quality-status").textContent = q.volume
    ? target
      ? "FEA quality targets met"
      : "FEA quality target missed. Refine sampling."
    : target
      ? c.preset === 1
        ? "Closed print surface. Check deviation."
        : "Closed preview surface."
      : "Mesh checks failed.";
  $("quality-status").style.color = target ? "var(--accent)" : "var(--orange)";
  const max = Math.max(...q.histogram);
  $("hist").replaceChildren(
    ...q.histogram.map((n, i) => {
      const bar = document.createElement("span");
      bar.style.height = `${(100 * n) / max}%`;
      bar.style.background = "#" + colorFor(i * 5 + 2.5).getHexString();
      bar.title = `${i * 5}–${(i + 1) * 5}°: ${fmt(n)} triangles`;
      return bar;
    }),
  );
  const diam = (a) => Math.sqrt((4 * a) / Math.PI) * c.cellSize;
  $("size-stats").textContent =
    `Equivalent triangle diameter, p01 / median / p99: ${num(diam(q.area_p01), 2)} / ${num(diam(q.area_median), 2)} / ${num(diam(q.area_p99), 2)} mm.`;
  $("stl").disabled = false;
  $("json").disabled = false;
  $("vtu").disabled = !m.tetrahedra.length;
}
function updateView() {
  if (!mesh) return;
  material.vertexColors = quality;
  material.color.set(quality ? 0xffffff : 0xc3dba6);
  material.needsUpdate = true;
  material.polygonOffset = wire || tetWire;
  material.polygonOffsetFactor = 1;
  material.polygonOffsetUnits = 1;
  if (wire && !wireMesh) {
    wireMesh = new THREE.LineSegments(
      new THREE.WireframeGeometry(mesh.geometry),
      new THREE.LineBasicMaterial({
        color: 0x15251b,
        transparent: true,
        opacity: 0.48,
        clippingPlanes: [plane],
      }),
    );
    group.add(wireMesh);
  }
  if (wireMesh) wireMesh.visible = wire;
  if (tetWire && !tetMesh) {
    const { points, tetrahedra } = current.mesh;
    const edges = new Set(),
      flat = [];
    for (const t of tetrahedra)
      for (let i = 0; i < 4; i++)
        for (let j = i + 1; j < 4; j++) {
          const a = Math.min(t[i], t[j]),
            b = Math.max(t[i], t[j]),
            key = `${a}:${b}`;
          if (!edges.has(key)) {
            edges.add(key);
            flat.push(
              ...points[a].map((x) => x / current.config.repeat),
              ...points[b].map((x) => x / current.config.repeat),
            );
          }
        }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.Float32BufferAttribute(flat, 3));
    tetMesh = new THREE.LineSegments(
      geo,
      new THREE.LineBasicMaterial({
        color: 0x9ec995,
        transparent: true,
        opacity: 0.36,
        clippingPlanes: [plane],
      }),
    );
    group.add(tetMesh);
  }
  if (tetMesh) tetMesh.visible = tetWire;
  mesh.visible = !tetWire;
  $("legend").hidden = !quality || tetWire;
  $("quality").classList.toggle("on", quality);
  $("solid").classList.toggle("on", !quality && !tetWire);
  $("wire").classList.toggle("on", wire);
  $("volume-wire").classList.toggle("on", tetWire);
}
$("quality").onclick = () => {
  quality = !quality;
  tetWire = false;
  updateView();
};
$("wire").onclick = () => {
  wire = !wire;
  updateView();
};
$("solid").onclick = () => {
  quality = false;
  tetWire = false;
  updateView();
};
$("volume-wire").onclick = () => {
  tetWire = !tetWire;
  updateView();
};
$("slice").oninput = () => {
  plane.constant = +$("slice").value / 100 - 0.5 + 0.001;
};
$("reset").onclick = () => {
  camera.position.set(1.9, 1.45, 2.1);
  controls.target.set(0, 0, 0);
  controls.update();
  $("slice").value = 100;
  plane.constant = 0.51;
};
function renderHistory() {
  const tbody = $("history");
  tbody.replaceChildren();
  for (const r of history) {
    const c = r.config,
      m = r.mesh,
      q = m.metrics;
    const row = document.createElement("tr");
    const cells = [
      `${names[c.shape]} / ${c.repeat}³${c.grade ? " graded" : ""}`,
      `${presets[c.preset]} / ${c.resolution}`,
      `${num(r.seconds, 3)} s`,
      fmt(m.triangles.length),
      fmt(m.tetrahedra.length),
      `${num(q.angle_p01, 1)}°`,
      `${small(q.centroid_distance_p95 * c.cellSize)} mm`,
      q.volume ? num(q.volume.minimum_quality) : "—",
    ];
    for (const value of cells) {
      const td = document.createElement("td");
      td.textContent = value;
      row.append(td);
    }
    const td = document.createElement("td"),
      button = document.createElement("button");
    button.textContent = "View";
    button.disabled = busy;
    button.onclick = () => {
      if (!busy) {
        show(r);
        status("");
      }
    };
    td.append(button);
    row.append(td);
    tbody.append(row);
  }
}
$("clear").onclick = () => {
  history = [];
  renderHistory();
};
function save(blob, suffix) {
  const c = current.config,
    url = URL.createObjectURL(blob),
    a = document.createElement("a");
  a.href = url;
  a.download = `${names[c.shape].toLowerCase().replaceAll(" ", "-")}-${presets[c.preset].toLowerCase().replaceAll(" ", "-")}.${suffix}`;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
$("stl").onclick = () => {
  if (current)
    save(
      new Blob([stlBytes(current.mesh, current.config.cellSize)], {
        type: "application/octet-stream",
      }),
      "stl",
    );
};
$("vtu").onclick = () => {
  if (current)
    save(
      new Blob([vtuText(current.mesh, current.config.cellSize)], {
        type: "application/xml",
      }),
      "vtu",
    );
};
$("json").onclick = () =>
  save(
    new Blob(
      [
        JSON.stringify(
          {
            configuration: current.config,
            seconds: current.seconds,
            metrics: current.mesh.metrics,
            units:
              "mesh coordinates in unit cells; display and mesh exports scaled by cellSize in mm",
          },
          null,
          2,
        ),
      ],
      { type: "application/json" },
    ),
    "json",
  );
labels();
launch(config());
