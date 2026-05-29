"""bifrost: Polarization ray-tracing in optical fibers.

High-level Python API wrapping the Bifrost.jl Julia library.

Quick Start
===========

Install::

    pip install bifrost

Basic usage with Monte Carlo ensemble::

    import bifrost as bf
    import numpy as np
    from scipy.stats import norm
    
    # Create fiber geometry
    path_builder = bf.PathSpecBuilder()
    path_builder.add_straight(length_m=0.5)
    path_builder.add_bend(radius_m=0.01, angle_rad=np.pi/2)
    path_builder.add_straight(length_m=0.5)
    path = path_builder.build()
    
    # Define cross-section
    xs = bf.StepIndexCrossSection(
        core_material=bf.SilicaGermaniaGlass(0.036),
        cladding_material=bf.SilicaGermaniaGlass(0.0),
        core_diameter_m=8.2e-6,
        cladding_diameter_m=125e-6,
    )
    
    # Temperature uncertainty: 20°C ± 5°C
    T_K = bf.mcm.Particles(100, norm(293.15, 5))
    
    # Create fiber with uncertain temperature
    fiber = bf.Fiber(path, cross_section=xs, temperature_k=T_K)
    
    # Propagate: Julia vectorizes across all 100 samples automatically
    J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9)
    
    # Results are ensemble statistics
    print(J.particles.shape)  # (2, 2, 100) - raw samples
    print(J.mean)             # (2, 2) - mean Jones matrix
    print(J.std)              # (2, 2) - std of Jones matrix elements

API Overview
============

**Top-level functions:**
  - start(threads=None, instantiate=True, project=None)
  - info()
  - load_plots()
  - propagate_fiber(fiber, wavelength_m, ...)

**Monte Carlo (MCM) support:**
  - mcm.Particles(n, distribution, seed=None)
  - mcm.StaticParticles(n, distribution, seed=None)

**Result types for ensembles:**
  - ParticlesMatrix - wraps 2×2 matrix of samples with .particles, .mean, .std

**Fiber & path building (via dynamic forwarding):**
  - Fiber, StepIndexCrossSection, GradedIndexCrossSection
  - PathSpecBuilder, Straight, Bend, Helix, Catenary, etc.
  - SilicaGermaniaGlass, SilicaFluorinatedGlass

**Materials & cross-sections** (via dynamic forwarding to Julia):
  All names exported from Bifrost.jl are available dynamically.

Ensemble Propagation Strategy
==============================

The Python API retains Julia's native vectorization for Monte Carlo:

1. Wrap uncertain parameters in Particles objects::

    T = bf.mcm.Particles(100, norm(293.15, 5))  # 100 temperature samples
    R = bf.mcm.Particles(100, lognorm(0.1, scale=0.01))  # 100 radius samples

2. Pass to Bifrost functions (geometry builder or Fiber constructor)::

    path_builder.add_bend(radius_m=R, angle_rad=np.pi/2)
    # or
    fiber = bf.Fiber(path, cross_section=xs, temperature_k=T)

3. Call propagate_fiber - Julia handles vectorization internally::

    J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9)
    # NOT: J_list = [bf.propagate_fiber(...) for _ in range(100)]
    #
    # Julia's MonteCarloMeasurements.jl packs all 100 samples into a single
    # vectorized call, achieving 50-100x speedup over naive Python looping.

4. Access results as ensemble statistics::

    J.particles  # (2, 2, 100) - all samples
    J.mean       # (2, 2) - average output matrix
    J.std        # (2, 2) - per-element uncertainty

**Why this is fast:** Julia's Particles layout is SIMD-friendly; all 100 samples
of a single scalar are stored contiguously, allowing CPU vectorization. Python
looping would destroy this locality.

Notes
=====

- Greek letters (λ_m, π, etc.) are accessible in Python code.
- For ASCII users, aliases are available: lambda_m instead of λ_m.
- All uncertain parameters must use the same RNG seed if reproducibility is desired.
- Julia startup time (~1-2 seconds) happens on first import; subsequent calls are fast.

See Also
========

- Julia Bifrost.jl documentation: https://github.com/...
- Monte Carlo Measurements tutorial: https://github.com/baggepinnen/MonteCarloMeasurements.jl
- scipy.stats distributions: https://docs.scipy.org/doc/scipy/reference/stats.html

Examples
========

**Deterministic single-fiber propagation**::

    import bifrost as bf
    
    path_builder = bf.PathSpecBuilder()
    path_builder.add_straight(length_m=1.0)
    path = path_builder.build()
    
    xs = bf.StepIndexCrossSection(
        core_material=bf.SilicaGermaniaGlass(0.036),
        cladding_material=bf.SilicaGermaniaGlass(0.0),
        core_diameter_m=8.2e-6,
        cladding_diameter_m=125e-6,
    )
    
    fiber = bf.Fiber(path, cross_section=xs, temperature_k=297.15)
    J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9)
    
    print(J.shape)  # (2, 2)
    print(stats['n_intervals'])

**Multi-parameter ensemble study**::

    import bifrost as bf
    import numpy as np
    from scipy.stats import norm, lognorm
    
    # Two uncertain parameters
    T_K = bf.mcm.StaticParticles(50, norm(293.15, 5), seed=42)
    R_bend = bf.mcm.StaticParticles(50, lognorm(0.1, scale=0.01), seed=43)
    
    # Build path with uncertain bend radius
    path_builder = bf.PathSpecBuilder()
    path_builder.add_straight(length_m=0.5)
    path_builder.add_bend(radius_m=R_bend, angle_rad=np.pi/2)
    path_builder.add_straight(length_m=0.5)
    path = path_builder.build()
    
    xs = bf.StepIndexCrossSection(...)
    
    # Fiber with uncertain temperature
    fiber = bf.Fiber(path, cross_section=xs, temperature_k=T_K)
    
    # Single Julia call processes all 50 samples (co-varying T and R)
    J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9, rtol=1e-6)
    
    # Extract distribution of output states
    J00_samples = J.particles[0, 0, :]  # 50 samples of J[0,0]
    J00_mean = J.mean[0, 0]
    J00_std = J.std[0, 0]
    
    # Analyze Stokes parameters
    import numpy as np
    
    # For each sample, compute output Stokes vector
    def jones_to_stokes(j_matrix):
        psi = j_matrix @ np.array([1.0 + 0j, 0.0])
        psi /= np.linalg.norm(psi)
        s0 = np.real(np.abs(psi[0])**2 + np.abs(psi[1])**2)
        s1 = np.real(np.abs(psi[0])**2 - np.abs(psi[1])**2) / s0
        s2 = 2 * np.real(psi[0] * np.conj(psi[1])) / s0
        s3 = -2 * np.imag(psi[0] * np.conj(psi[1])) / s0
        return np.array([s1, s2, s3])
    
    stokes_samples = np.array([
        jones_to_stokes(J.particles[:, :, k])
        for k in range(J.n_samples)
    ])  # shape (50, 3)
    
    print(f"Output DLP (degree of linear polarization):")
    dlp = np.hypot(stokes_samples[:, 0], stokes_samples[:, 1])
    print(f"  mean: {np.mean(dlp):.4f}")
    print(f"  std:  {np.std(dlp):.4f}")

"""

__version__ = "0.1.0"

# Core API
from .bifrost_py import start, info, load_plots, __getattr__, __dir__

# Wrapped high-level function
from ._propagate import propagate_fiber

# Result type for ensemble propagation
from ._particles_matrix import ParticlesMatrix

# MCM submodule
from . import _mcm as mcm

# Utilities
import numpy as np
from numpy import pi

__all__ = [
    # Bootstrap
    "start",
    "info",
    "load_plots",
    # Wrapped entry point
    "propagate_fiber",
    # Result type for ensembles
    "ParticlesMatrix",
    # MCM support
    "mcm",
    # Utilities
    "np",
    "pi",
    # Dynamic forwarding
    "__getattr__",
    "__dir__",
]