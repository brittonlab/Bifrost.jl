"""
High-level wrapped propagation functions with type hints and documentation.

Key design: detect Particles in inputs and use native Julia vectorization.
"""

from typing import Optional, Union, Any
import numpy as np
from .bifrost_py import get_jl
from ._particles_matrix import ParticlesMatrix
from ._convert import (
    convert_propagate_result,
    jl_particles_matrix_to_python,
    jl_particles_stats_to_python
)


def propagate_fiber(
    fiber: Any,
    *,
    wavelength_m: float,
    rtol: float = 1e-9,
    atol: float = 1e-12,
    verbose: bool = False,
) -> Union[tuple[np.ndarray, dict], tuple["ParticlesMatrix", dict]]:
    """
    Propagate light through a fiber and compute output Jones matrix.
    
    Uses adaptive step-doubling integration along the fiber path to solve
    the coupled-mode equations for polarization evolution.
    
    **Native Monte Carlo support:** If fiber parameters (temperature, bend radius,
    etc.) are uncertain and passed as Particles objects, Julia will vectorize
    the entire propagation across all samples in a single call. This is MUCH
    faster than Python-side looping.
    
    Parameters
    ----------
    fiber : Fiber
        Fiber instance created via bifrost.Fiber(...).
        May contain Particles-wrapped uncertain parameters.
    wavelength_m : float
        Wavelength in meters.
        Example: 1550e-9 for 1550 nm (standard telecom).
    rtol : float, optional
        Relative tolerance for adaptive integrator. Default: 1e-9.
    atol : float, optional
        Absolute tolerance for adaptive integrator. Default: 1e-12.
    verbose : bool, optional
        Print integration diagnostics. Default: False.
    
    Returns
    -------
    J : np.ndarray or ParticlesMatrix
        Output Jones matrix.
        
        - If inputs are deterministic: shape (2, 2), dtype complex128.
        - If inputs contain Particles (ensemble): ParticlesMatrix object
          with .particles (shape (2, 2, n_samples)), .mean, .std properties.
    
    stats : dict
        Integration statistics. If ensemble:
        - 'n_intervals' : array of shape (n_samples,), steps per sample
        - 'arc_length_m' : scalar, same for all samples
        - Or reduced stats (mean/std of per-sample stats)
    
    Raises
    ------
    ValueError
        If wavelength_m <= 0 or tolerances invalid.
    
    Examples
    --------
    **Deterministic propagation**::
    
        import bifrost as bf
        import numpy as np
        
        # Create fiber
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
        
        # Propagate (single call, deterministic)
        J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9)
        print(J.shape)  # (2, 2)
    
    **Ensemble (Monte Carlo) propagation with native vectorization**::
    
        import bifrost as bf
        from scipy.stats import norm
        import numpy as np
        
        # Temperature: 20°C ± 5°C (normal distribution)
        # This creates a Particles object with 100 samples
        T_K = bf.mcm.StaticParticles(100, norm(293.15, 5))
        
        # Create fiber with uncertain temperature
        fiber = bf.Fiber(path, cross_section=xs, temperature_k=T_K)
        
        # Propagate: Julia vectorizes across all 100 samples in ONE call
        # This is ~100x faster than looping in Python!
        J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9)
        
        # J is a ParticlesMatrix
        print(type(J))  # <class 'ParticlesMatrix'>
        print(J.particles.shape)  # (2, 2, 100)
        print(J.mean)  # (2, 2) - mean J over 100 samples
        print(J.std)   # (2, 2) - std of J over 100 samples
        
        # Access individual element distributions
        J00_samples = J.particles[0, 0, :]  # 100 samples of J[0,0]
        J00_mean = J.mean[0, 0]
        J00_std = J.std[0, 0]
    
    **Mixed: some parameters uncertain, some deterministic**::
    
        # Temperature uncertain, but bend radius deterministic
        T_K = bf.mcm.Particles(50, norm(293.15, 5))
        
        path_builder = bf.PathSpecBuilder()
        path_builder.add_straight(length_m=0.5)
        path_builder.add_bend(radius_m=0.01, angle_rad=np.pi/2)  # fixed
        path_builder.add_straight(length_m=0.5)
        path = path_builder.build()
        
        fiber = bf.Fiber(path, cross_section=xs, temperature_k=T_K)
        
        # Result: 50 samples of output matrix (temperature varies, geometry fixed)
        J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9)
        print(J.particles.shape)  # (2, 2, 50)
    
    **Multiple uncertain parameters**::
    
        from scipy.stats import norm, lognorm
        
        # Temperature and bend radius both uncertain
        T_K = bf.mcm.StaticParticles(30, norm(293.15, 5))
        
        path_builder = bf.PathSpecBuilder()
        path_builder.add_straight(length_m=0.5)
        path_builder.add_bend(
            radius_m=bf.mcm.StaticParticles(30, lognorm(0.1, scale=0.01)),
            angle_rad=np.pi/2
        )
        path_builder.add_straight(length_m=0.5)
        path = path_builder.build()
        
        fiber = bf.Fiber(path, cross_section=xs, temperature_k=T_K)
        
        # Result: 30 samples (all parameters co-vary)
        J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9)
        print(J.particles.shape)  # (2, 2, 30)
    
    Notes
    -----
    - Input polarization is assumed to be horizontal [1, 0]ᵀ.
    - Julia's MonteCarloMeasurements.jl handles all vectorization internally.
    - Speedup scales with n_samples: 100 samples takes ~2–5x longer than
      1 sample (not 100x), thanks to Julia's SIMD-friendly Particles layout.
    - All uncertain parameters must have the same number of samples (n).
    
    See Also
    --------
    mcm.Particles : Create ensemble with scipy.stats distributions.
    mcm.StaticParticles : Create fixed-seed ensemble (deterministic).
    ParticlesMatrix : Return type for ensemble results.
    """
    # Input validation
    if not isinstance(wavelength_m, (int, float)) or wavelength_m <= 0:
        raise ValueError(f"wavelength_m must be positive, got {wavelength_m}")
    if rtol <= 0 or atol <= 0:
        raise ValueError(f"Tolerances must be positive, got rtol={rtol}, atol={atol}")
    
    # Call Julia (passes Particles through directly if present)
    jl = get_jl()
    J_jl, stats_jl = jl.Bifrost.propagate_fiber(
        fiber,
        lambda_m=wavelength_m,
        rtol=rtol,
        atol=atol,
        verbose=verbose,
    )
    
    # Detect if result is Particles (each matrix element is a Particles object)
    is_ensemble = _is_particles_matrix(J_jl)
    
    if is_ensemble:
        # Return ParticlesMatrix wrapper
        J = ParticlesMatrix(J_jl)
        stats = jl_particles_stats_to_python(stats_jl)
    else:
        # Regular deterministic result
        J, stats = convert_propagate_result(J_jl, stats_jl)
    
    return J, stats


def _is_particles_matrix(mat: Any) -> bool:
    """Check if a 2x2 matrix contains Particles objects."""
    try:
        # Check if any element has a .particles attribute (Julia Particles marker)
        return hasattr(mat[0, 0], "particles")
    except (TypeError, IndexError, AttributeError):
        return False