```@meta
CurrentModule = Bifrost
```

# BIFROST

BIFROST (Birefringence In Fiber: Research and Optical Simulation Toolkit) is a Julia
package for simulating polarization mode dispersion in step-index, silica-based,
germania-doped optical fibers.

Rather than representing a fiber as a pre-sliced list of Jones matrices, BIFROST builds a
continuous centerline path, binds it to a transverse fiber cross section, and integrates a
local Jones generator

```math
\frac{dJ}{ds}=K(s,\omega)\,J,\qquad J(s_0)=I,
```

with an adaptive, phase-insensitive, exponential-midpoint Lie-group integrator that never
steps across path breakpoints.

## Architecture

The package is one umbrella module, `Bifrost`, composed of layered submodules:

| Module | Role |
| --- | --- |
| [`MaterialProperties`](@ref material-properties) | Intrinsic optical material properties and spectral responses. |
| [`PathGeometry`](@ref path-geometry-api) | 3D path authoring, differential geometry, and twist resolution. |
| [`FiberCS`](@ref fiber-cross-section) | Step-index cross-section optics and local birefringence responses. |
| [`FiberPath`](@ref fiber-path) | Binds a path to a cross section and assembles `K(s)` / `Kω(s)`. |
| [`PathIntegral`](@ref path-integral) | Adaptive Jones and DGD-sensitivity propagation. |
| [`Plots`](@ref plots) | Opt-in visualization helpers (`using Bifrost.Plots`). |

The **Theory** pages derive the generators and explain the integral geometric quantities.
The **API** pages list the documented, exported symbols of each module.

## Quick start

```julia
using Bifrost

xs = StepIndexCrossSection(
    GermaniaSilicaGlass(0.036),
    GermaniaSilicaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like",
)

spec = PathSpecBuilder()
straight!(spec; length = 0.5, meta = [Nickname("lead-in")])
bend!(spec; radius = 0.05, angle = pi / 2, meta = [Nickname("90 deg bend")])
straight!(spec; length = 0.5, meta = [Nickname("lead-out")])

path  = build(spec)
fiber = Fiber(path; cross_section = xs, T_ref_K = 297.15)

J, stats = propagate_fiber(fiber; λ_m = 1550e-9, rtol = 1e-9, verbose = false)
```

For differential group delay (DGD):

```julia
J, G, stats = propagate_fiber_sensitivity(fiber; λ_m = 1550e-9, rtol = 1e-9)
dgd = output_dgd(J, G)   # use output_dgd_2x2(J, G) on MonteCarloMeasurements paths
```

See [Examples](@ref) for the full Julia and Python walkthroughs, [Usage](@ref) for the
authoring model, and §1 of `test/human/bifrost-demos.ipynb` for the smallest runnable
example.
Installation instructions live in the project `README.md`.

## Regime of Operation

This library models step-index silica-based germania-doped optical fibers. It
includes chromatic dispersion effects. At this time, the library does not model
other possible dopants, specially engineered materials such as
dispersion-compensating fiber, or other index profiles such as graded-index
fibers.

At present, BIFROST models birefringence from four mechanisms:

- Core noncircularity
- Asymmetric thermal stress due to differing coefficients of thermal expansion
  between core and cladding when the core is noncircular
- Bending
- Twisting

It does not model birefringence due to:

- Cladding noncircularity
- Non-concentric cladding and core
- External asymmetric stress, such as pushing on the fiber in one direction
- Transverse electric fields
- Axial magnetic fields

Based on validation work, as well as the limits of the approximations made and
the validity range of the data used in BIFROST, we believe the codebase
correctly computes supported contributions to birefringence in the following
regime:

- Single-mode operation, $V<2.405$
- The weakly guiding regime, $n_{\text{co}}-n_{\text{cl}} \ll 1$, which
  implicitly requires weak germanium doping
- The nearly circular-core regime, $e^2 \ll 1$
- Bend radii much larger than the cladding radius, $R \gg r_{\text{cl}}$
- Temperatures 200 K $\lesssim T \lesssim$ 300 K, limited by the model for
  the thermo-optic coefficient $dn/dT$ of bulk germania glass
- Telecom wavelengths 1 $\mu$m $\lesssim \lambda \lesssim$ 2 $\mu$m

We do not model the temperature dependence of the coefficients of thermal
expansion or the photoelastic constants $p_{11}$ and $p_{12}$ in fused
silica and germania, as the variation is small within the above parameter
regime. Polarization-dependent loss and nonlinear scattering effects are also
outside the current model.

### Known limitations

- `src/path-integral.jl` assumes lossless Jones propagation. Do not introduce
  gain, loss, or polarization-dependent loss there without a separate design.
- `modify(fiber)` handles geometry and metadata perturbations, including
  thermal scaling, but twist remapping through modification remains a known
  caveat area.
- The Julia code has substantial internal tests, but it is not yet a validated
  replacement for the legacy Python model or for published fiber data.
