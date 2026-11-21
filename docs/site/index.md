---
hide:
  - toc
---

# Tetrahedral mesh explorer

Generate a periodic solid and inspect its tetrahedral mesh. Drag to orbit, scroll to zoom.

<div class="mesh-lab">
  <div class="mesh-stage">
    <img class="mesh-poster" src="assets/explorer-gyroid.webp" alt="Gyroid sheet with triangulated walls and box caps">
    <canvas tabindex="0" aria-label="3D mesh. Drag or use arrow keys to rotate. Use plus and minus to zoom." hidden></canvas>
    <div class="mesh-key"><span>Implicit wall</span><span>Box cap</span></div>
  </div>
  <form class="mesh-controls">
    <h2>Mesh settings</h2>
    <fieldset class="mesh-fields">
      <legend>Implicit field</legend>
      <div class="mesh-segments">
        <label><input type="radio" name="field" value="0" checked><span>Gyroid</span></label>
        <label><input type="radio" name="field" value="1"><span>Schwarz P</span></label>
      </div>
    </fieldset>
    <div class="mesh-slider">
      <label for="mesh-cells">Cells per axis</label><output for="mesh-cells">10</output>
      <input id="mesh-cells" type="range" min="6" max="16" value="10">
    </div>
    <div class="mesh-slider">
      <label for="mesh-width">Band half-width</label><output for="mesh-width">0.50</output>
      <input id="mesh-width" type="range" min="0.30" max="0.80" step="0.05" value="0.50">
    </div>
    <div class="mesh-display">
      <label class="mesh-switch"><input id="mesh-lines" type="checkbox" checked><span>Triangle edges</span></label>
      <button type="button" id="mesh-reset" disabled>Reset view</button>
    </div>
  </form>
  <p class="mesh-status" role="status">Preparing the 3D view…</p>
</div>
<noscript>Enable JavaScript to rotate and regenerate the mesh.</noscript>

<div class="language-python highlight"><pre><code id="mesh-code"><span class="kn">import</span> <span class="mesh-module">meshers</span> <span class="k">as</span> <span class="mesh-module">mrs</span>

mesh = <span class="mesh-module">mrs</span>.<span class="nf">generate</span>(
    <span class="s2" data-param="field">"gyroid"</span>, cells=<span class="mi" data-param="cells">10</span>, band=(-<span class="mf" data-param="width">0.5</span>, <span class="mf" data-param="width">0.5</span>),
    periodic=(<span class="kc">True</span>,) * <span class="mi">3</span>, geometry_tolerance=<span class="mf">0.05</span>, optimize_passes=<span class="mi">2</span>, threads=<span class="mi">0</span>,
)</code></pre></div>

## Explore the guides

<div class="mesh-gallery" markdown="1">

[![Tetrahedral sphere boundary](assets/sphere.svg) **Generate your first mesh**](getting-started.md)

[![Named planar cut through a sphere](assets/cut-sphere.svg) **Label an intersection**](intersections.md)

[![Displacement and energy convergence](assets/playground/displacement.svg) **Verify a finite element solution**](solver-verification.md)

[![Gyroid cutaway colored by tetrahedron quality](assets/gyroid-quality-section.webp) **Inspect boundary quality**](quality.md)

</div>

<script type="module" src="assets/explorer.js"></script>
