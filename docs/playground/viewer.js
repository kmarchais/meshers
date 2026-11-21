import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { RoomEnvironment } from "three/addons/environments/RoomEnvironment.js";

const lab = document.querySelector(".mesh-lab");
const stage = lab.querySelector(".mesh-stage");
const status = lab.querySelector(".mesh-status");
const canvas = stage.querySelector("canvas");
const settings = lab.querySelector("fieldset");
const field = document.querySelector("#mesh-field");
const cells = document.querySelector("#mesh-cells");
const width = document.querySelector("#mesh-width");
const lines = document.querySelector("#mesh-lines");

function start() {
  const renderer = new THREE.WebGLRenderer({ canvas, alpha: true, antialias: true });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 1.22;
  const scene = new THREE.Scene();
  const studio = new RoomEnvironment();
  const generator = new THREE.PMREMGenerator(renderer);
  scene.environment = generator.fromScene(studio, 0.04).texture;
  studio.dispose();
  generator.dispose();
  const key = new THREE.DirectionalLight(0xffebcf, 3);
  key.position.set(2, 4, 3);
  scene.add(key);
  const camera = new THREE.PerspectiveCamera(33, 1, 0.05, 20);
  camera.position.set(1.85, 1.25, 2.05);
  const orbit = new OrbitControls(camera, canvas);
  orbit.enablePan = false;
  orbit.minDistance = 1.25;
  orbit.maxDistance = 5;
  orbit.update();
  orbit.saveState();
  const render = () => renderer.render(scene, camera);
  orbit.addEventListener("change", render);
  new ResizeObserver(() => {
    const { width, height } = stage.getBoundingClientRect();
    renderer.setSize(width, height, false);
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
    render();
  }).observe(stage);
  document.querySelector("#mesh-reset").onclick = () => orbit.reset();
  canvas.addEventListener("keydown", event => {
    const pose = new THREE.Spherical().setFromVector3(camera.position);
    switch (event.key) {
      case "ArrowLeft": pose.theta -= 0.15; break;
      case "ArrowRight": pose.theta += 0.15; break;
      case "ArrowUp": pose.phi -= 0.15; break;
      case "ArrowDown": pose.phi += 0.15; break;
      case "+": case "=": pose.radius /= 1.1; break;
      case "-": pose.radius *= 1.1; break;
      default: return;
    }
    event.preventDefault();
    pose.makeSafe();
    camera.position.setFromSpherical(pose);
    orbit.update();
  });
  const materials = [
    new THREE.MeshStandardMaterial({color: 0x389d90, roughness: 0.33, metalness: 0.48, side: THREE.DoubleSide}),
    new THREE.MeshStandardMaterial({color: 0xd38b4c, roughness: 0.41, metalness: 0.55, side: THREE.DoubleSide}),
  ];
  const ink = new THREE.LineBasicMaterial({color: 0x12292a, transparent: true, opacity: 0.23});
  const solid = new THREE.Mesh(new THREE.BufferGeometry(), materials);
  const edges = new THREE.LineSegments(new THREE.BufferGeometry(), ink);
  solid.position.setScalar(-0.5);
  edges.position.copy(solid.position);
  scene.add(solid, edges);
  lines.onchange = () => { edges.visible = lines.checked; render(); };

  const worker = new Worker(new URL("../mesh-worker.js", import.meta.url));
  let latest = 0;
  let timer;
  let busy = false;
  function request() {
    if (busy) return;
    busy = true;
    worker.postMessage({ id: latest, field: Number(field.value), cells: Number(cells.value), width: Number(width.value) });
  }
  function change() {
    latest++;
    lab.querySelector('output[for="mesh-cells"]').value = cells.value;
    lab.querySelector('output[for="mesh-width"]').value = Number(width.value).toFixed(2);
    status.textContent = "Generating mesh…";
    clearTimeout(timer);
    timer = setTimeout(request, 100);
  }
  worker.onmessage = ({ data }) => {
    busy = false;
    if (data.id !== latest) { clearTimeout(timer); request(); return; }
    if (data.error) { status.textContent = "Mesh rejected. " + data.error; return; }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(data.coords, 3));
    geometry.setIndex(new THREE.BufferAttribute(data.indices, 1));
    geometry.addGroup(0, data.wallCount, 0);
    geometry.addGroup(data.wallCount, data.indices.length - data.wallCount, 1);
    geometry.computeVertexNormals();
    solid.geometry.dispose();
    solid.geometry = geometry;
    edges.geometry.dispose();
    edges.geometry = new THREE.WireframeGeometry(geometry);
    canvas.hidden = false;
    stage.querySelector("img").hidden = true;
    status.textContent = data.tets.toLocaleString() + " tetrahedra · " + data.nodes.toLocaleString() + " nodes";
    document.querySelector("#mesh-code").textContent =
      'import meshers as mrs\n\nmesh = mrs.generate(\n    "' + (field.value === "0" ? "gyroid" : "schwarz_p") +
      '", cells=' + cells.value + ', band=(-' + width.value + ', ' + width.value + '),\n' +
      '    periodic=(True,) * 3, geometry_tolerance=0.05, optimize_passes=0,\n)';
    render();
  };
  worker.onerror = () => {
    status.textContent = "The meshing worker could not start. Reload to try again.";
    settings.disabled = true;
  };
  [field, cells, width].forEach(input => input.addEventListener("input", change));
  settings.disabled = false;
  change();
}

try {
  start();
} catch (error) {
  status.textContent = "The 3D viewer could not start. Check WebGL 2 support or try the Python example.";
  console.error(error);
}
