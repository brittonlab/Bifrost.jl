# Usage

## Building Blocks

These files are intentionally useful on their own:

| File | Standalone role |
| --- | --- |
| `material-properties.jl` | Together with the other files in `src/material/`, provides material constants and spectra; no path or fiber geometry. |
| `path-geometry.jl` | Three-dimensional path construction and geometric queries; no optics. |
| `path-integral.jl` | Adaptive propagation for callable `K(s)` and `Kω(s)` generators. |

The fiber-specific layers combine those pieces:

| File | How it extends the standalone pieces |
| --- | --- |
| `cross-section.jl` | Together with the other files in `src/fiber-cross-section/`, adds step-index fiber optics and birefringence responses. |
| `fiber-path.jl` | Binds path geometry to a cross section and assembles bend/twist `K` and `Kω`. |

## Current Model

High-level authoring is path based:

1. Build geometry with `PathSpecBuilder`.
2. Freeze and place it with `build(...)`, producing a `PathSpecCached`.
3. Bind that path to a `FiberCrossSection` *subtype* with `Fiber(path; cross_section,
   T_ref_K)`.
4. Propagate at a requested wavelength with `propagate_fiber(fiber; λ_m=...)`.

The current path primitives include:

- `StraightSegment`
- `BendSegment`
- `CatenarySegment`
- `HelixSegment`
- `JumpBy`
- `JumpTo`

`JumpBy` and `JumpTo` are authoring conveniences. At build time they are
resolved into a G2 quintic connector implemented in
`src/geometry/path-geometry-connector.jl`.

Material twist is attached as per-segment metadata using `Twist`. A twist run
starts on the segment carrying the annotation and continues until the next
twist annotation or the end of the path. Twist rates may be constant or
callable functions of run-local arc length.

The optical `Fiber` stores:

- the built `PathSpecCached`,
- the `FiberCrossSection`,
- a reference temperature `T_ref_K`,
- the fiber domain `[s_start, s_end]`.

The operating wavelength is not stored on `Fiber`. It is supplied per query to
`generator_K(fiber, λ_m)`, `generator_Kω(fiber, λ_m)`, `propagate_fiber`, and
`propagate_fiber_sensitivity`.

For runnable code, see [Examples](@ref). The [Theory](@ref "Deriving the Generators")
pages cover the generator formalism, propagation, and DGD sensitivity.
