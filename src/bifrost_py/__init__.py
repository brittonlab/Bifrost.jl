"""bifrost: Polarization ray-tracing in optical fibers.

High-level Python API wrapping the Bifrost.jl Julia library.

Quick Start
===========

Install::

    pip install bifrost

Basic usage with deterministic fiber::

    import bifrost as bf
    import numpy as np
    
    # Create path: straight → bend → straight
    spec = bf.SubpathBuilder()
    bf.start_b(spec)
    bf.straight_b(spec, length_m=0.5)
    bf.bend_b(spec, radius_m=0.05, angle_rad=np.pi/2)
    bf.straight_b(spec, length_m=0.5)
    bf.seal_b(spec)
    path = bf.build(spec)
    
    # Define cross-section
    xs = bf.StepIndexCrossSection(
        core_material=bf.SilicaGermaniaGlass(0.036),
        cladding_material=bf.SilicaGermaniaGlass(0.0),
        core_diameter_m=8.2e-6,
        cladding_diameter_m=125e-6,
    )
    
    # Create fiber
    fiber = bf.Fiber(path, cross_section=xs, T_ref_K=297.15)
    
    # Propagate
    J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9)
    print(J.shape)  # (2, 2)

With Monte Carlo uncertainties::

    import bifrost as bf
    import numpy as np
    from scipy.stats import norm
    
    # Temperature uncertainty: 20°C ± 5°C
    T_K = bf.mcm.Particles(100, norm(293.15, 5))
    
    # Build path with meta (Nickname for documentation)
    spec = bf.SubpathBuilder(spin_rate=lambda s: np.sin(2*np.pi*s/100))
    bf.start_b(spec)
    bf.straight_b(spec, length_m=0.5, meta=[bf.Nickname("lead-in")])
    bf.bend_b(spec, radius_m=0.05, angle_rad=np.pi/2,
              meta=[bf.Nickname("90° bend")])
    bf.straight_b(spec, length_m=0.5, meta=[bf.Nickname("lead-out")])
    bf.seal_b(spec)
    path = bf.build(spec)
    
    xs = bf.StepIndexCrossSection(...)
    
    # Fiber with uncertain temperature (MCMadd meta)
    spec2 = bf.SubpathBuilder()
    bf.start_b(spec2)
    bf.helix_b(spec2, radius_m=0.025, pitch_m=0.05, turns=1000,
               meta=[bf.Nickname("temperature-sensitive helix"),
                     bf.MCMadd('T_K', norm(0, 5))])
    bf.straight_b(spec2, length_m=5.0)
    bf.seal_b(spec2)
    path2 = bf.build(spec2)
    
    fiber = bf.Fiber(path2, cross_section=xs, T_ref_K=T_K)
    
    # Propagate: Julia vectorizes across all 100 samples automatically
    J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9)
    
    # Results are ensemble statistics
    print(J.particles.shape)  # (2, 2, 100) - raw samples
    print(J.mean)             # (2, 2) - mean Jones matrix
    print(J.std)              # (2, 2) - std of Jones matrix elements

API Overview
============

**Path Building:**
  - SubpathBuilder() — create a mutable path specification
  - start_b(spec; point=(0,0,0), outgoing_tangent=(0,0,1), ...)
  - straight_b(spec; length_m, meta=None)
  - bend_b(spec; radius_m, angle_rad, axis_angle_rad=0, meta=None)
  - helix_b(spec; radius_m, pitch_m, turns, axis_angle_rad=0, meta=None)
  - catenary_b(spec; a_m, length_m, axis_angle_rad=0, meta=None)
  - seal_b(spec; extra_m=0, meta=None)
  - build(spec) → PathBuilt (Julia object)

**Fiber Creation:**
  - Fiber(path, cross_section, T_ref_K=None)

**Propagation:**
  - propagate_fiber(fiber, wavelength_m, ...)

**Monte Carlo (MCM):**
  - mcm.Particles(n, distribution, seed=None)
  - mcm.StaticParticles(n, distribution, seed=None)

**Metadata Annotations:**
  - Nickname(label) — human-readable segment label
  - MCMadd(symbol, distribution) — additive perturbation
  - MCMmul(symbol, distribution) — multiplicative perturbation

**Result Types:**
  - ParticlesMatrix — wraps 2×2 matrix of samples with .particles, .mean, .std

**Dynamic Forwarding:**
  All names exported from Bifrost.jl are available dynamically:
  - Fiber, StepIndexCrossSection, GradedIndexCrossSection
  - SilicaGermaniaGlass, SilicaFluorinatedGlass, etc.

Ensemble Propagation Strategy
==============================

The Python API retains Julia's native vectorization for Monte Carlo:

1. **Wrap uncertain parameters in Particles objects**::

    T = bf.mcm.Particles(100, norm(293.15, 5))  # 100 temperature samples

2. **Pass to Fiber constructor or use MCM meta**::

    fiber = bf.Fiber(path, cross_section=xs, T_ref_K=T)
    # or use meta on segments:
    bf.helix_b(spec, radius_m=0.025, turns=1000,
               meta=[bf.MCMadd('T_K', norm(0, 5))])

3. **Call propagate_fiber — Julia handles vectorization**::

    J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9)

4. **Access results as ensemble statistics**::

    J.particles  # (2, 2, 100) — all 100 samples
    J.mean       # (2, 2) — mean output matrix
    J.std        # (2, 2) — std of each element

**Why efficient:** Julia's Particles layout is SIMD-friendly; all 100 samples
of a scalar are contiguous, enabling CPU vectorization. Python looping would
destroy this locality.

Metadata System
===============

The geometry layer stores and applies metadata annotations per-segment.
Three types are built-in:

**Nickname(label)** — Document and label segments for plotting/diagnostics::

    meta = [bf.Nickname("90° bend"), bf.Nickname("lead-in")]
    bf.bend_b(spec, radius_m=0.05, angle_rad=np.pi/2, meta=meta)

**MCMadd(symbol, distribution)** — Additive perturbation::

    # Temperature uncertainty: 5 K centered at reference
    T_unc = norm(loc=0, scale=5)
    meta = bf.MCMadd('T_K', T_unc)
    bf.helix_b(spec, radius_m=0.025, turns=1000, meta=[meta])

    # Result: perturbed_value = baseline_value + T_unc_sample

**MCMmul(symbol, distribution)** — Multiplicative scale factor::

    # Radius uncertainty: log-normal ±10%
    radius_scale = lognorm(s=0.1, scale=1.0)
    meta = bf.MCMmul('radius', radius_scale)
    bf.bend_b(spec, radius_m=0.05, angle_rad=np.pi/2, meta=[meta])

    # Result: perturbed_radius = baseline_radius * radius_scale_sample

Per-segment MCM meta is applied **once at build time** by the Julia layer,
not interpreted by Python. This keeps all MCM logic in Julia where it belongs.

Examples
========

**Simple deterministic fiber**::

    import bifrost as bf
    import numpy as np
    
    spec = bf.SubpathBuilder()
    bf.start_b(spec)
    bf.straight_b(spec, length_m=1.0)
    bf.seal_b(spec)
    path = bf.build(spec)
    
    xs = bf.StepIndexCrossSection(
        core_material=bf.SilicaGermaniaGlass(0.036),
        cladding_material=bf.SilicaGermaniaGlass(0.0),
        core_diameter_m=8.2e-6,
        cladding_diameter_m=125e-6,
    )
    
    fiber = bf.Fiber(path, cross_section=xs, T_ref_K=297.15)
    J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9)

**Path with labeled segments**::

    import bifrost as bf
    import numpy as np
    
    spec = bf.SubpathBuilder()
    bf.start_b(spec)
    bf.straight_b(spec, length_m=0.5, 
                  meta=[bf.Nickname("lead-in")])
    bf.bend_b(spec, radius_m=0.05, angle_rad=np.pi/2, 
              meta=[bf.Nickname("90° bend")])
    bf.straight_b(spec, length_m=0.5, 
                  meta=[bf.Nickname("lead-out")])
    bf.seal_b(spec)
    path = bf.build(spec)

**Temperature-uncertain helix**::

    import bifrost as bf
    import numpy as np
    from scipy.stats import norm
    
    # Temperature perturbation: ±5 K centered at reference (20°C)
    T_perturbation = norm(loc=0, scale=5)
    
    spec = bf.SubpathBuilder()
    bf.start_b(spec)
    bf.helix_b(spec, radius_m=0.025, pitch_m=0.05, turns=1000,
               meta=[bf.Nickname("temperature-sensitive"),
                     bf.MCMadd('T_K', T_perturbation)])
    bf.straight_b(spec, length_m=5.0, 
                  meta=[bf.Nickname("reference")])
    bf.seal_b(spec)
    path = bf.build(spec)
    
    xs = bf.StepIndexCrossSection(...)
    
    # Build 100-sample ensemble with temperature variation
    T_ensemble = bf.mcm.Particles(100, norm(293.15, 5))
    fiber = bf.Fiber(path, cross_section=xs, T_ref_K=T_ensemble)
    
    # Single call processes all 100 samples
    J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9)
    print(f"Mean output power (dBm): {20*np.log10(np.abs(J.mean[0,0])):.2f}")
    print(f"Uncertainty (dBm): {20*np.log10(J.std[0,0]):.2f}")

**Multi-parameter uncertainty study**::

    import bifrost as bf
    import numpy as np
    from scipy.stats import norm, lognorm
    
    # Two independent uncertainty sources
    T_unc = norm(loc=0, scale=5)          # Temperature: ±5 K
    radius_unc = lognorm(s=0.1, scale=1)  # Radius: ±10% log-normal
    
    spec = bf.SubpathBuilder()
    bf.start_b(spec)
    
    # Temperature-sensitive helix
    bf.helix_b(spec, radius_m=0.025, pitch_m=0.05, turns=1000,
               meta=[bf.Nickname("T-sensitive"),
                     bf.MCMadd('T_K', T_unc)])
    
    # Bend with uncertain radius
    bf.bend_b(spec, radius_m=0.05, angle_rad=np.pi/2,
              meta=[bf.Nickname("variable bend"),
                    bf.MCMmul('radius', radius_unc)])
    
    bf.straight_b(spec, length_m=5.0)
    bf.seal_b(spec)
    path = bf.build(spec)
    
    xs = bf.StepIndexCrossSection(...)
    
    # Build ensemble with both T and radius varying (same 100 samples)
    T_ensemble = bf.mcm.Particles(100, norm(293.15, 5), seed=42)
    fiber = bf.Fiber(path, cross_section=xs, T_ref_K=T_ensemble)
    
    J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9)
    
    # Correlate output power with temperature
    J00_samples = J.particles[0, 0, :]  # 100 samples
    T_samples = T_ensemble.particles
    correlation = np.corrcoef(np.abs(J00_samples), T_samples)[0, 1]
    print(f"J[0,0] vs Temperature correlation: {correlation:.3f}")

Notes
=====

- All path building uses Python functions (start_b, straight_b, etc.)
- Metadata is attached via the `meta=` keyword argument to any builder function
- Multiple metadata items can be passed as a list
- Julia names (λ_m, π, etc.) are available as bf.lambda_m, bf.pi, etc.
- MCM perturbations defined on segments are applied by Julia at build time
- For reproducible ensembles, pass a `seed=` to mcm.Particles()

See Also
========

- Julia Bifrost.jl documentation: https://github.com/...
- Monte Carlo Measurements tutorial: https://github.com/baggepinnen/MonteCarloMeasurements.jl
- scipy.stats distributions: https://docs.scipy.org/doc/scipy/reference/stats.html

"""

__version__ = "0.1.0"

# Utilities
import numpy as np
from numpy import pi

# Core API
from .bifrost_py import start, info, load_plots, __getattr__, __dir__

# Path builders with meta support
from ._builder_wrappers import start_b, straight_b, bend_b, helix_b, catenary_b, seal_b

# Wrapped high-level functions
from ._propagate import propagate_fiber
from ._fiber_proxy import Fiber

# Result type for ensemble propagation
from ._particles_matrix import ParticlesMatrix

# MCM submodule
from . import _mcm as mcm

# Meta annotation system
from ._meta import AbstractMeta, Nickname, MCMadd, MCMmul

__all__ = [
    # Bootstrap
    "start",
    "info",
    "load_plots",
    # Wrapped entry point
    "propagate_fiber",
    "Fiber",
    # Result type for ensembles
    "ParticlesMatrix",
    # MCM support
    "mcm",
    # Metadata system
    "AbstractMeta",
    "Nickname",
    "MCMadd",
    "MCMmul",
    # Path builders
    "start_b",
    "straight_b", 
    "bend_b",
    "helix_b",
    "catenary_b",
    "seal_b",
    # Utilities
    "np",
    "pi",
    # Dynamic forwarding
    "__getattr__",
    "__dir__",
]