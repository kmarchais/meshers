/** Instantiate the experimental TPMS package from a URL or wasm bytes. */
export async function createMeshers(
  source = new URL("./meshers.wasm", import.meta.url),
) {
  const bytes =
    source instanceof Uint8Array || source instanceof ArrayBuffer
      ? source
      : await (await fetch(source)).arrayBuffer();
  const { instance } = await WebAssembly.instantiate(bytes, {});
  const e = instance.exports;
  return {
    generate({
      shape = 0,
      preset = 0,
      repeat = 1,
      resolution = 12,
      thickness = 0.6,
      grade = 0,
    } = {}) {
      if (![shape, preset, repeat, resolution].every(Number.isInteger))
        throw new Error(
          "Shape, preset, repeat and resolution must be integers.",
        );
      const length = e.meshers_generate(
        shape,
        preset,
        repeat,
        resolution,
        thickness,
        grade,
      );
      try {
        const output = JSON.parse(
          new TextDecoder().decode(
            new Uint8Array(e.memory.buffer, e.meshers_output_ptr(), length),
          ),
        );
        if (!output.ok) throw new Error(output.error);
        return output.mesh;
      } finally {
        e.meshers_clear();
      }
    },
  };
}
