"""
Lightweight wrapper for Fiber constructor.

Converts Python Particles objects to Julia equivalents before passing to Julia.
"""

from typing import Any, Optional
from .bifrost_py import get_jl
from ._meta import AbstractMeta, validate_meta, _meta_list_to_julia


def Fiber(path, cross_section, T_ref_K=None):
    """
    Create a Fiber with automatic Particles and meta conversion.
    
    Parameters
    ----------
    path : PathBuilt
        Built path from SubpathBuilder.build()
    cross_section : StepIndexCrossSection or GradedIndexCrossSection
        Fiber cross-section
    T_ref_K : float or Particles, optional
        Reference temperature in Kelvin. Can be a scalar or Particles ensemble.
    
    Returns
    -------
    Fiber (Julia object via juliacall)
    
    Examples
    --------
    >>> spec = bf.SubpathBuilder()
    >>> bf.start_b(spec)
    >>> bf.straight_b(spec, length_m=1.0)
    >>> bf.seal_b(spec)
    >>> path = bf.build(spec)
    
    >>> xs = bf.StepIndexCrossSection(...)
    
    >>> # Deterministic fiber:
    >>> fiber = bf.Fiber(path, cross_section=xs, T_ref_K=297.15)
    
    >>> # Fiber with uncertain temperature (100 samples):
    >>> T = bf.mcm.Particles(100, norm(297.15, 5))
    >>> fiber_mcm = bf.Fiber(path, cross_section=xs, T_ref_K=T)
    
    >>> # Propagate
    >>> J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9)
    """
    jl = get_jl()
    
    # Convert Particles to Julia type if present
    if T_ref_K is not None and hasattr(T_ref_K, '_julia_type'):
        jl_temperature = T_ref_K._julia_type
    else:
        jl_temperature = T_ref_K
    
    # Call Julia Fiber constructor
    return jl.Bifrost.Fiber(
        path,
        cross_section=cross_section,
        T_ref_K=jl_temperature,
    )