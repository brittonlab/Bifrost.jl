```@meta
CollapsedDocStrings = true
```

# [PathGeometry](@id path-geometry-api)

Three-dimensional path authoring and differential geometry: straight, bend, catenary,
helix, and `JumpBy` / `JumpTo` connectors, plus material-twist resolution, sampling, and
global path diagnostics. This module also carries the per-segment metadata vocabulary
(`Nickname`, `MCMadd`, `MCMmul`). The Theory section covers the integral geometric
quantities (turning angle, torsion, material twist, frame rotation).

```@autodocs
Modules = [Bifrost.PathGeometry]
Filter = is_public_api
```
