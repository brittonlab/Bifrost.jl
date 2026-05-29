"""
Return type converters: Julia → Python/numpy.

When Julia functions return complex objects (matrices, named tuples, etc.),
convert them to standard Python types for better usability.
"""

from typing import Any, Union, TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
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
    """
    arr = np.asarray(jl_matrix)
    
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
    """
    result = {}
    
    # Extract keys and values from Julia NamedTuple.
    # Julia's _fields returns Symbols, so we index by position to be safe.
    try:
        fields = jl_nt._fields
        n_fields = len(fields)
        
        for i in range(n_fields):
            # Get key as string
            key = str(fields[i])
            # Index by position (more reliable than getattr with converted strings)
            val = jl_nt[i]
            
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
            
            result[key] = val
            
    except (AttributeError, TypeError, IndexError):
        # Fallback: try dict-like access
        try:
            for key in jl_nt.keys():
                key_str = str(key)
                val = jl_nt[key]
                result[key_str] = val
        except (AttributeError, TypeError):
            pass
    
    return result


def jl_vector_of_structs_to_list(jl_vec: Any) -> list:
    """
    Convert a Julia Vector of structs to a Python list of dicts.
    
    Parameters
    ----------
    jl_vec : Any
        Julia Vector where each element is a struct (has ._fields).
    
    Returns
    -------
    list
        List of dicts, one per struct.
    
    Examples
    --------
    Julia:  [PropagatorStats(1, 0), PropagatorStats(2, 1)]
    Python: [{'accepted_steps': 1, 'rejected_steps': 0}, 
             {'accepted_steps': 2, 'rejected_steps': 1}]
    """
    result = []
    try:
        for item in jl_vec:
            if hasattr(item, "_fields"):
                # It's a struct-like object
                result.append(jl_namedtuple_to_dict(item))
            elif hasattr(item, "__array__"):
                result.append(jl_matrix_to_numpy(item))
            else:
                result.append(item)
    except (TypeError, AttributeError):
        pass
    
    return result


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
    from ._particles_matrix import ParticlesMatrix
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
        keys = [str(k) for k in stats_jl._fields]
    except AttributeError:
        keys = [str(k) for k in stats_jl.keys()] if hasattr(stats_jl, 'keys') else []
    
    for key in keys:
        try:
            val = getattr(stats_jl, key)
        except (AttributeError, TypeError):
            try:
                val = stats_jl[key]
            except (KeyError, TypeError, IndexError):
                continue
        
        # If Particles, reduce to mean
        if hasattr(val, "particles"):
            particles_array = np.asarray(val.particles)
            result[key] = np.mean(particles_array)
            result[key + "_std"] = np.std(particles_array)
        elif hasattr(val, "__array__"):
            result[key] = np.asarray(val)
        else:
            result[key] = val
    
    return result


def convert_propagate_result(J_jl: Any, stats_jl: Any) -> tuple[np.ndarray, Union[dict, list]]:
    """
    Convert propagate_fiber return values to Python types.
    
    Parameters
    ----------
    J_jl : Any
        Julia Jones matrix (2×2 complex).
    stats_jl : Any
        Julia stats. Can be:
        - A NamedTuple (single struct)
        - A Vector of structs (list of PropagatorStats)
    
    Returns
    -------
    J : np.ndarray
        Shape (2, 2), dtype complex128.
    stats : dict or list
        If NamedTuple: dict with keys like 'n_intervals', 'arc_length_m', etc.
        If Vector of structs: list of dicts, one per interval.
    """
    J = jl_matrix_to_numpy(J_jl)
    
    # Detect if stats_jl is a vector or a NamedTuple
    if hasattr(stats_jl, "__len__") and not hasattr(stats_jl, "_fields"):
        # It's a vector/list (has __len__ but no _fields)
        stats = jl_vector_of_structs_to_list(stats_jl)
    else:
        # It's a NamedTuple or dict-like
        stats = jl_namedtuple_to_dict(stats_jl)
    
    return J, stats