import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { RoomEnvironment } from "three/addons/environments/RoomEnvironment.js";

export function createViewer(stage) {
  const canvas = stage.querySelector("canvas");
  const renderer = new THREE.WebGLRenderer({
    canvas, alpha: false, antialias: false, powerPreference: "default",
  });
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 1.22;
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x263b3b);
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
  let interacting = false;
  let reduced = false;
  let previous = 0;
  let late = 0;
  let settle;
  let width = 0;
  let height = 0;
  function sizeBuffer() {
    if (!width || !height) return;
    const ratio = reduced
      ? Math.min(devicePixelRatio, 0.72, Math.sqrt(200_000 / (width * height)))
      : Math.min(devicePixelRatio, 1.5, Math.sqrt(1_200_000 / (width * height)));
    if (canvas.width !== Math.floor(width * ratio) || canvas.height !== Math.floor(height * ratio)) {
      renderer.setDrawingBufferSize(width, height, ratio);
    }
  }
  function render() {
    if (frame !== undefined) return;
    frame = requestAnimationFrame(now => {
      frame = undefined;
      if (interacting && !reduced) {
        const gap = now - previous;
        // Consecutive late frames during motion justify a smaller drawing buffer.
        late = previous && gap > 22 && gap < 100 ? late + 1 : 0;
        previous = now;
        if (late >= 3) {
          reduced = true;
          sizeBuffer();
        }
      }
      renderer.render(scene, camera);
    });
  }
  orbit.addEventListener("change", render);
  orbit.addEventListener("start", () => {
    clearTimeout(settle);
    interacting = true;
  });
  orbit.addEventListener("end", () => {
    // Wheel events come in bursts, each with its own start and end.
    settle = setTimeout(() => {
      interacting = false;
      previous = late = 0;
      if (reduced) {
        reduced = false;
        sizeBuffer();
        render();
      }
    }, 180);
  });
  const resize = new ResizeObserver(([entry]) => {
    ({ width, height } = entry.contentRect);
    if (!width || !height) return;
    sizeBuffer();
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
    new THREE.MeshStandardMaterial({color: 0x389d90, roughness: 0.33, metalness: 0.48, polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1}),
    new THREE.MeshStandardMaterial({color: 0xd38b4c, roughness: 0.41, metalness: 0.55, polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1}),
  ];
  const ink = new THREE.LineBasicMaterial({color: 0x12292a, transparent: true, opacity: 0.23, depthWrite: false});
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
      clearTimeout(settle);
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
