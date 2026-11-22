export function fillSlider(input) {
  const fraction = (input.valueAsNumber - Number(input.min)) / (Number(input.max) - Number(input.min));
  input.style.setProperty("--fill", fraction * 100 + "%");
}
