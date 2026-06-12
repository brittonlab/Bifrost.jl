// Keep module- and type-level docstrings expanded even on pages that set
// `CollapsedDocStrings = true`. Documenter collapses every `details.docstring`
// on load; this reopens the entries whose category is "Module" or "Type", and
// runs after that collapse because `load` fires after Documenter's `ready`
// handler.
window.addEventListener("load", function () {
  var EXPANDED_CATEGORIES = ["Module", "Type"];
  document.querySelectorAll("details.docstring").forEach(function (d) {
    var category = d.querySelector(".docstring-category");
    if (category && EXPANDED_CATEGORIES.indexOf(category.textContent.trim()) !== -1) {
      d.setAttribute("open", "true");
    }
  });
});
