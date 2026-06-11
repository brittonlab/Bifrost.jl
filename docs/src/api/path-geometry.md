# [PathGeometry](@id path-geometry-api)

Three-dimensional path authoring and differential geometry: straight, bend, catenary,
helix, and `JumpBy` / `JumpTo` connectors, plus material-twist resolution, sampling, and
global path diagnostics. The transverse frame exposed by `normal`/`binormal`/`frame` is
the parallel-transported (Bishop) pair. This module also carries the per-segment
metadata vocabulary (`Nickname`, `MCMadd`, `MCMmul`). The Theory section covers the
integral geometric quantities (turning angle, torsion, material twist) and the frame
convention.

```@autodocs
Modules = [Bifrost.PathGeometry]
```
