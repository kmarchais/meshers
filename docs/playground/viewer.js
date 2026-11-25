import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { RoomEnvironment } from "three/addons/environments/RoomEnvironment.js";

export function createViewer(stage) {
  const canvas = stage.querySelector("canvas");
  const renderer = new THREE.WebGLRenderer({
    canvas, alpha: true, antialias: true, powerPreference: "low-power",
  });
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 1.22;
  const scene = new THREE.Scene();
  const studio = new RoomEnvironment();
  const generator = new THREE.PMREMGenerator(renderer);
  const environment = generator.fromScene(studio, 0.04);
  scene.environment = environment.texture;
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
  // High-rate pointers can deliver several orbit updates between two paints.
  let frame;
  function render() {
    if (frame !== undefined) return;
    frame = requestAnimationFrame(() => {
      frame = undefined;
      renderer.render(scene, camera);
    });
  }
  orbit.addEventListener("change", render);
  const resize = new ResizeObserver(() => {
    const { width, height } = stage.getBoundingClientRect();
    if (!width || !height) return;
    const ratio = Math.min(devicePixelRatio, 1.5, Math.sqrt(1_200_000 / (width * height)));
    renderer.setDrawingBufferSize(width, height, ratio);
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
    render();
  });
  resize.observe(stage);
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
  function showMesh(data, showEdges) {
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
    edges.visible = showEdges;
    canvas.hidden = false;
    stage.querySelector("img").hidden = true;
    render();
  }
  return {
    showMesh,
    reset() { orbit.reset(); },
    setEdges(visible) { edges.visible = visible; render(); },
    dispose() {
      cancelAnimationFrame(frame);
      resize.disconnect();
      orbit.dispose();
      solid.geometry.dispose();
      edges.geometry.dispose();
      for (const material of materials) material.dispose();
      ink.dispose();
      environment.dispose();
      renderer.dispose();
      renderer.forceContextLoss();
      const replacement = canvas.cloneNode(false);
      replacement.hidden = true;
      canvas.replaceWith(replacement);
      stage.querySelector("img").hidden = false;
    },
  };
}
