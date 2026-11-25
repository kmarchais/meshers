// The projection is fixed. Each face records its back and front tetrahedron.
export async function mountShapeRange(figure) {
  const canvas = figure.querySelector("canvas");
  const context = canvas.getContext("2d");
  if (!context) return;
  const response = await fetch(new URL("./gyroid-shape.json", import.meta.url));
  if (!response.ok) throw new Error("Shape data could not load.");
  const { points, faces, scores } = await response.json();
  const low = figure.querySelector("#shape-min");
  const high = figure.querySelector("#shape-max");
  const track = figure.querySelector(".range-track");
  const count = figure.querySelector('[role="status"]');
  const colors = scores.map((q) => {
    const mix = (a, b) => Math.round(a + q * (b - a));
    return `rgb(${mix(199, 45)} ${mix(100, 159)} ${mix(60, 143)})`;
  });
  context.scale(2, 2);
  context.strokeStyle = "#31463f";
  context.lineWidth = 0.55;
  let frame;

  function draw() {
    frame = undefined;
    const min = low.valueAsNumber,
      max = high.valueAsNumber;
    const selected = scores.map((q) => q >= min && q <= max);
    context.clearRect(0, 0, canvas.width / 2, canvas.height / 2);
    for (const [a, b, c, back, front] of faces) {
      if (!selected[back] || selected[front]) continue;
      context.beginPath();
      context.moveTo(...points[a]);
      context.lineTo(...points[b]);
      context.lineTo(...points[c]);
      context.closePath();
      context.fillStyle = colors[back];
      context.fill();
      context.stroke();
    }
    const total = selected.filter(Boolean).length;
    count.textContent = `${total.toLocaleString()} of ${scores.length.toLocaleString()} cells in range.`;
    canvas.setAttribute(
      "aria-label",
      `Gyroid tetrahedra with shape quality from ${low.value} to ${high.value}. ${count.textContent}`,
    );
  }
  function update(handle) {
    if (low.valueAsNumber > high.valueAsNumber) handle.value = handle === low ? high.value : low.value;
    low.setAttribute("aria-valuemax", high.value);
    high.setAttribute("aria-valuemin", low.value);
    for (const input of [low, high]) {
      figure.querySelector(`output[for="${input.id}"]`).value = input.valueAsNumber.toFixed(2);
    }
    track.style.setProperty("--start", low.valueAsNumber * 100 + "%");
    track.style.setProperty("--end", high.valueAsNumber * 100 + "%");
    if (frame === undefined) frame = requestAnimationFrame(draw);
  }
  for (const input of [low, high]) input.addEventListener("input", () => update(input));
  // The thumbs are native controls. A track click moves the nearest endpoint,
  // including when both thumbs coincide at either end.
  track.addEventListener("pointerdown", (event) => {
    if (event.target !== track) return;
    const rect = track.getBoundingClientRect();
    const value = Math.max(0, Math.min(1, (event.clientX - rect.left) / rect.width));
    const handle = value < (low.valueAsNumber + high.valueAsNumber) / 2 ? low : high;
    handle.value = value.toFixed(2);
    handle.focus();
    update(handle);
  });
  figure.querySelector("img").hidden = true;
  canvas.hidden = false;
  figure.querySelector("fieldset").hidden = false;
  update(low);
}
