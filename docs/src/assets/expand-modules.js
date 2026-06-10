// Keep module-level docstrings expanded even on pages that set
// `CollapsedDocStrings = true`. Documenter collapses every `details.docstring`
// on load; this reopens just the entries whose category is "Module", which run
// after that collapse because `load` fires after Documenter's `ready` handler.
window.addEventListener("load", function () {
  document.querySelectorAll("details.docstring").forEach(function (d) {
    var category = d.querySelector(".docstring-category");
    if (category && category.textContent.trim() === "Module") {
      d.setAttribute("open", "true");
    }
  });
});
