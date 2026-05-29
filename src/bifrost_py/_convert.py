"""
Return type converters: Julia → Python/numpy.

When Julia functions return complex objects (matrices, named tuples, etc.),
convert them to standard Python types for better usability.
"""

from typing import Any, Union
import numpy as np
from ._particles_matrix import ParticlesMatrix


def jl_matrix_to_numpy(jl_matrix: Any) -> np.ndarray:
    """
    Convert Julia matrix to numpy array.
    
    Parameters
    ----------
    jl_matrix : Any
        Julia matrix object (from juliacall).
    
    Returns
    -------
    np.ndarray
        Numpy array with appropriate dtype (complex128 if complex input).
    
    Examples
    --------
    >>> J_jl = jl.Bifrost.some_function()  # Returns Julia matrix
    >>> J = jl_matrix_to_numpy(J_jl)
    >>> print(type(J), J.dtype)
    <class 'numpy.ndarray'> complex128
    """
    # juliacall matrices implement numpy's array interface
    arr = np.asarray(jl_matrix)
    
    # Ensure proper dtype
    if np.iscomplexobj(arr):
        arr = arr.astype(np.complex128)
    else:
        arr = arr.astype(np.float64)
    
    return arr


def jl_namedtuple_to_dict(jl_nt: Any) -> dict:
    """
    Convert Julia NamedTuple to Python dict.
    
    Parameters
    ----------
    jl_nt : Any
        Julia NamedTuple object.
    
    Returns
    -------
    dict
        Dictionary with same keys/values (recursively converts nested structures).
    
    Examples
    --------
    >>> stats_jl = jl.Bifrost.get_stats()  # Returns NamedTuple
    >>> stats = jl_namedtuple_to_dict(stats_jl)
    >>> print(stats["n_intervals"])
    42
    """
    result = {}
    
    # Julia NamedTuple exposes keys via indexing
    try:
        keys = list(jl_nt._fields)  # Most reliable way
    except AttributeError:
        try:
            keys = list(jl_nt.keys())
        except (AttributeError, TypeError):
            # Fallback: try to iterate
            keys = list(jl_nt)
    
    for key in keys:
        val = getattr(jl_nt, key, None) or jl_nt[key]
        
        # Recursively convert nested structures
        if hasattr(val, "_fields"):  # Another NamedTuple
            val = jl_namedtuple_to_dict(val)
        elif hasattr(val, "__array__"):  # Matrix-like
            val = jl_matrix_to_numpy(val)
        elif isinstance(val, (list, tuple)):
            val = [
                jl_matrix_to_numpy(v) if hasattr(v, "__array__") else v
                for v in val
            ]
        
        result[str(key)] = val
    
    return result


def convert_propagate_result(J_jl: Any, stats_jl: Any) -> tuple[np.ndarray, dict]:
    """
    Convert propagate_fiber return values to Python types.
    
    Parameters
    ----------
    J_jl : Any
        Julia Jones matrix (2×2 complex).
    stats_jl : Any
        Julia NamedTuple with propagation statistics.
    
    Returns
    -------
    J : np.ndarray
        Shape (2, 2), dtype complex128.
    stats : dict
        Dictionary with keys like 'n_intervals', 'arc_length_m', etc.
    """
    J = jl_matrix_to_numpy(J_jl)
    stats = jl_namedtuple_to_dict(stats_jl)
    return J, stats


def jl_particles_to_numpy(jl_p: Any) -> np.ndarray:
    """
    Extract samples from a Julia Particles object.
    
    Parameters
    ----------
    jl_p : Any
        Julia Particles object (from MonteCarloMeasurements.jl).
    
    Returns
    -------
    np.ndarray
        1D array of samples.
    """
    # Julia Particles exposes .particles field
    return np.asarray(jl_p.particles, dtype=np.complex128)


def jl_particles_matrix_to_python(jl_mat: Any) -> "ParticlesMatrix":
    """
    Convert a 2x2 matrix of Julia Particles objects to a Python ParticlesMatrix.
    
    Parameters
    ----------
    jl_mat : Any
        2x2 matrix where each element is a Particles object.
    
    Returns
    -------
    ParticlesMatrix
        Wrapper providing .particles, .mean, .std properties.
    """
    return ParticlesMatrix(jl_mat)


def jl_particles_stats_to_python(stats_jl: Any) -> dict:
    """
    Convert Julia stats (which may contain Particles) to Python dict.
    
    For each field that is a Particles object, extract samples and compute
    mean/std. For scalar fields, keep as-is.
    
    Parameters
    ----------
    stats_jl : Any
        Julia NamedTuple with stats (some fields may be Particles).
    
    Returns
    -------
    dict
        Reduced statistics: {field: mean or value}.
    """
    result = {}
    
    try:
        keys = list(stats_jl._fields)
    except AttributeError:
        keys = list(stats_jl.keys()) if hasattr(stats_jl, 'keys') else []
    
    for key in keys:
        val = getattr(stats_jl, key, None) or stats_jl[key]
        
        # If Particles, reduce to mean
        if hasattr(val, "particles"):
            particles_array = np.asarray(val.particles)
            result[str(key)] = np.mean(particles_array)
            result[str(key) + "_std"] = np.std(particles_array)
        elif hasattr(val, "__array__"):
            result[str(key)] = np.asarray(val)
        else:
            result[str(key)] = val
    
    return result


def convert_propagate_result(
    J_jl: Any,
    stats_jl: Any
) -> tuple[np.ndarray, dict]:
    """
    Convert deterministic propagate_fiber result to numpy.
    
    (For ensemble results, use ParticlesMatrix directly.)
    """
    J = jl_matrix_to_numpy(J_jl)
    stats = jl_namedtuple_to_dict(stats_jl)
    return J, stats