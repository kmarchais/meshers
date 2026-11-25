// BashLexer leaves external command names and their options as plain text.
for (const block of document.querySelectorAll(
  ".language-sh code, .language-bash code, .language-console code",
)) {
  const walker = document.createTreeWalker(block, NodeFilter.SHOW_TEXT);
  const nodes = [];
  while (walker.nextNode()) nodes.push(walker.currentNode);
  for (const node of nodes) {
    if (node.parentElement.closest(".c,.c1,.s,.s1,.s2,a")) continue;
    const parts = node.textContent.split(
      /(\B--?[a-zA-Z][\w-]*\b|\b(?:cargo|uv|rustup|python|pytest|ruff|maturin|pip)\b)/g,
    );
    if (parts.length === 1) continue;
    const fragment = document.createDocumentFragment();
    for (const [i, part] of parts.entries()) {
      if (i % 2 === 0) fragment.append(part);
      else {
        const span = document.createElement("span");
        span.className = part.startsWith("-") ? "cli-option" : "cli-command";
        span.textContent = part;
        fragment.append(span);
      }
    }
    node.replaceWith(fragment);
  }
}
