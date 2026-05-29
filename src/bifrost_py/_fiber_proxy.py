"""
Lightweight wrapper for Fiber constructor.

Converts Python Particles objects to Julia equivalents before passing to Julia.
"""

from typing import Any, Optional
from .bifrost_py import get_jl


def Fiber(path, cross_section, T_ref_K=None):
    """
    Create a Fiber with automatic Particles conversion.
    
    Parameters
    ----------
    path : Path
        PathSpecCached from builder.build()
    cross_section : StepIndexCrossSection or GradedIndexCrossSection
        Fiber cross-section
    temperature_k : float or Particles, optional
        Reference temperature in Kelvin. Can be a scalar or Particles ensemble.
    
    Returns
    -------
    Fiber (Julia object via juliacall)
    
    Examples
    --------
    >>> fiber = bf.Fiber(path, cross_section=xs, temperature_k=297.15)
    >>> # With uncertain temperature:
    >>> T = bf.mcm.Particles(50, norm(297.15, 5))
    >>> fiber_mcm = bf.Fiber(path, cross_section=xs, temperature_k=T)
    """
    jl = get_jl()
    
    # Convert Particles to Julia type if present
    if T_ref_K is not None and hasattr(T_ref_K, '_julia_type'):
        print("Yeah, it's got it")
        jl_temperature = T_ref_K._julia_type
    else:
        jl_temperature = T_ref_K
    
    # Call Julia Fiber constructor
    return jl.Bifrost.Fiber(
        path,
        cross_section=cross_section,
        T_ref_K=jl_temperature,
    )