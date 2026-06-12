"""
Meta annotation system for fiber path geometry.

Provides Python equivalents of Julia's AbstractMeta types:
  - Nickname: human-readable labels for segments
  - MCMadd: additive Monte Carlo perturbations (baseline + perturbation)
  - MCMmul: multiplicative Monte Carlo perturbations (baseline * scale)

These are attached to segments via the `meta=` parameter in builder functions
and serialized to Julia at build time.

Usage
=====

Import the types::

    import bifrost as bf
    from scipy.stats import norm
    import numpy as np

Attach to segments using the `meta=` parameter::

    spec = bf.SubpathBuilder()
    bf.start_b(spec)
    
    # Label a segment
    bf.straight_b(spec, length_m=0.5, 
                  meta=[bf.Nickname("lead-in")])
    
    # Add temperature uncertainty (additive)
    T_unc = norm(loc=0, scale=5)  # ±5K centered at reference
    bf.helix_b(spec, radius_m=0.025, pitch_m=0.05, turns=1000,
               meta=[bf.MCMadd('T_K', T_unc),
                     bf.Nickname("temperature-sensitive")])
    
    # Add multiplicative scaling (e.g., bend radius ±10%)
    radius_scale = lognorm(s=0.1, scale=1.0)
    bf.bend_b(spec, radius_m=0.05, angle_rad=np.pi/2,
              meta=[bf.MCMmul('radius', radius_scale)])
    
    bf.seal_b(spec)
    path = bf.build(spec)

When passed to `bf.Fiber()`, the meta is carried through to Julia, where
MCMadd/MCMmul perturbations are applied during build at the geometry layer.
Foreign meta (labels the geometry layer doesn't interpret) are stored and
available to consuming layers (e.g., the fiber analyzer).

Perturbation Convention
=======================

Consumers combine all MCMadd and MCMmul entries as::

    perturbed = baseline * product(all_MCMmul) + sum(all_MCMadd)

Examples:

- **Additive:** MCMadd('T_K', norm(0, 5)) means the field gets ±5 K uncertainty
  centered at the baseline value (the reference temperature).

- **Multiplicative:** MCMmul('length', 0.95) means the length is scaled by 0.95
  (i.e., 5% shorter). MCMmul('length', lognorm(0.1)) means ±10% log-normal
  variation around the baseline length.

"""

from abc import ABC, abstractmethod
from typing import Any, Union, Callable, List
import numpy as np


class AbstractMeta(ABC):
    """
    Base class for per-segment annotations.
    
    Subclasses carry metadata that the geometry layer stores verbatim and
    applies according to their type. Foreign meta (annotations this layer
    doesn't interpret) are carried through inertly for consuming layers.
    """
    
    @abstractmethod
    def to_julia(self):
        """Convert to a Julia object for serialization across the bridge."""
        pass
    
    def __repr__(self):
        return f"{self.__class__.__name__}({self._repr_args()})"
    
    @abstractmethod
    def _repr_args(self) -> str:
        """Return string representation of constructor arguments."""
        pass


class Nickname(AbstractMeta):
    """
    Attach a human-readable label to a path segment.
    
    Plotting and diagnostic code use this metadata for display labels. The 
    geometry layer stores it in the meta bag and consuming layers use it for
    visualization and reporting.
    
    Parameters
    ----------
    label : str
        The human-readable label for the segment.
        
    Example
    -------
    >>> meta = bf.Nickname("temperature-sensitive helix")
    >>> bf.helix_b(spec, radius_m=0.025, pitch_m=0.05, turns=1000, 
    ...            meta=[meta])
    """
    
    def __init__(self, label: str):
        if not isinstance(label, str):
            raise TypeError(f"Nickname label must be str, got {type(label)}")
        self.label = label
    
    def to_julia(self):
        """Return Julia Nickname object for serialization."""
        return ("Nickname", self.label)
    
    def _repr_args(self) -> str:
        return repr(self.label)


class MCMadd(AbstractMeta):
    """
    Attach an additive Monte Carlo perturbation to a segment field.
    
    Consumers combine matching entries as `baseline * product(mul) + sum(add)`;
    MCMadd contributes to the additive sum. The distribution may be any object
    the consumer knows how to sample, including a scalar, Particles ensemble,
    or scipy.stats distribution.
    
    The perturbation is applied **once at Julia build time**, not in Python.
    Julia's geometry layer handles all MCM logic for optimal vectorization.
    
    Parameters
    ----------
    symbol : str or Symbol
        The field name to perturb (e.g., 'T_K', ':T_K', or just T_K).
        This should match a field name on a Julia segment type.
    distribution : float, Particles, or scipy.stats distribution
        The additive perturbation. For a scalar, the field is offset by that
        amount. For a distribution (scipy.stats or Particles), each sample
        provides a different offset.
        
    Example
    -------
    >>> # Temperature uncertainty: ±5 K centered at reference
    >>> from scipy.stats import norm
    >>> T_unc = norm(loc=0, scale=5)
    >>> meta = bf.MCMadd('T_K', T_unc)
    >>> bf.helix_b(spec, radius_m=0.025, pitch_m=0.05, turns=1000, 
    ...            meta=[meta])
    
    >>> # With Particles ensemble (100 samples):
    >>> T_ensemble = bf.mcm.Particles(100, norm(293.15, 5))
    >>> # Note: MCMadd on a segment is for offsets; use T_ensemble for Fiber():
    >>> fiber = bf.Fiber(path, cross_section=xs, T_ref_K=T_ensemble)
    """
    
    def __init__(self, symbol: Union[str, Any], distribution: Any):
        self.symbol = _normalize_symbol(symbol)
        self.distribution = distribution
    
    def to_julia(self):
        """Return Julia MCMadd object for serialization."""
        dist_julia = _distribution_to_julia(self.distribution)
        return ("MCMadd", self.symbol, dist_julia)
    
    def _repr_args(self) -> str:
        dist_repr = _distribution_repr(self.distribution)
        return f":{self.symbol}, {dist_repr}"


class MCMmul(AbstractMeta):
    """
    Attach a multiplicative Monte Carlo perturbation to a segment field.
    
    Consumers combine matching entries as `baseline * product(mul) + sum(add)`;
    MCMmul contributes a direct scale factor to the multiplicative product.
    
    The perturbation is applied **once at Julia build time**, not in Python.
    Julia's geometry layer handles all MCM logic for optimal vectorization.
    
    Parameters
    ----------
    symbol : str or Symbol
        The field name to scale (e.g., 'radius', ':radius').
        This should match a field name on a Julia segment type.
    distribution : float, Particles, or scipy.stats distribution
        The multiplicative scale factor. A value of 1.0 means no change;
        0.95 means 5% reduction; 1.1 means 10% increase.
        Note: MCMmul(:length, -0.4) flips the sign and shortens the segment.
        
    Example
    -------
    >>> # Bend radius uncertainty: ±10% log-normal
    >>> from scipy.stats import lognorm
    >>> radius_scale = lognorm(s=0.1, scale=1.0)
    >>> meta = bf.MCMmul('radius', radius_scale)
    >>> bf.bend_b(spec, radius_m=0.05, angle_rad=np.pi/2, meta=[meta])
    
    >>> # Deterministic scale-down
    >>> meta = bf.MCMmul('pitch', 0.5)  # Half the pitch
    >>> bf.helix_b(spec, radius_m=0.025, pitch_m=0.05, turns=1000, 
    ...            meta=[meta])
    """
    
    def __init__(self, symbol: Union[str, Any], distribution: Any):
        self.symbol = _normalize_symbol(symbol)
        self.distribution = distribution
    
    def to_julia(self):
        """Return Julia MCMmul object for serialization."""
        dist_julia = _distribution_to_julia(self.distribution)
        return ("MCMmul", self.symbol, dist_julia)
    
    def _repr_args(self) -> str:
        dist_repr = _distribution_repr(self.distribution)
        return f":{self.symbol}, {dist_repr}"


# ============================================================================
# Helpers for normalization and serialization
# ============================================================================

def _normalize_symbol(symbol: Union[str, Any]) -> str:
    """
    Normalize a symbol to a string without leading colon.
    
    Accepts:
      - 'T_K' (str without colon)
      - ':T_K' (str with leading colon)
      - Symbol objects (converted to string)
      - Any object with __name__ (converted to string)
    
    Returns the symbol as a string without the leading colon.
    """
    if isinstance(symbol, str):
        return symbol.lstrip(':')
    elif hasattr(symbol, 'name'):
        return str(symbol.name).lstrip(':')
    else:
        s = str(symbol).lstrip(':')
        return s


def _distribution_to_julia(distribution: Any) -> Any:
    """
    Convert a Python distribution to a Julia-compatible form.
    
    Handles:
      - float/int scalars → pass through
      - Particles objects → extract Julia representation
      - scipy.stats distributions → pass through (Julia bridge handles)
      - numpy arrays → convert to list
      - functions → pass through
    """
    try:
        from . import _mcm
        if isinstance(distribution, (_mcm.Particles, _mcm.StaticParticles)):
            return distribution._julia_type
    except (ImportError, AttributeError):
        pass
    
    if isinstance(distribution, (int, float, complex)):
        return float(distribution)
    
    if isinstance(distribution, np.ndarray):
        return distribution.tolist()
    
    if callable(distribution):
        return distribution
    
    return distribution


def _distribution_repr(distribution: Any) -> str:
    """Return a readable representation of a distribution."""
    if isinstance(distribution, (int, float)):
        return repr(distribution)
    elif isinstance(distribution, np.ndarray):
        return f"array{distribution.shape}"
    elif callable(distribution):
        if hasattr(distribution, 'dist') and hasattr(distribution, 'args'):
            return f"{distribution.dist.name}(...)"
        elif hasattr(distribution, '__name__'):
            return distribution.__name__
        else:
            return "<function>"
    else:
        class_name = distribution.__class__.__name__
        if hasattr(distribution, '__repr__'):
            full_repr = repr(distribution)
            if len(full_repr) > 60:
                return f"{class_name}(...)"
            return full_repr
        return class_name


def _meta_list_to_julia(meta_list: Union[List[AbstractMeta], None]) -> Any:
    """
    Convert a Python meta list to Julia representation.
    
    Parameters
    ----------
    meta_list : list of AbstractMeta or None
        The metadata list from Python.
        
    Returns
    -------
    list
        A list of tuples ready for Julia bridge serialization.
    """
    if meta_list is None:
        return []
    
    if not isinstance(meta_list, (list, tuple)):
        meta_list = [meta_list]
    
    return [m.to_julia() for m in meta_list]


# ============================================================================
# Utilities for working with meta in the path builder
# ============================================================================

def validate_meta(meta: Any) -> List[AbstractMeta]:
    """
    Validate and normalize a meta argument.
    
    Accepts:
      - None → []
      - Single AbstractMeta → [meta]
      - List of AbstractMeta → meta (validated)
    
    Raises TypeError if meta contains non-AbstractMeta objects.
    
    Parameters
    ----------
    meta : None, AbstractMeta, or list of AbstractMeta
        The input metadata.
        
    Returns
    -------
    list of AbstractMeta
        Validated metadata list (empty if input was None).
        
    Raises
    ------
    TypeError
        If meta contains items that are not AbstractMeta subclasses.
    """
    if meta is None:
        return []
    
    if isinstance(meta, AbstractMeta):
        return [meta]
    
    if isinstance(meta, (list, tuple)):
        for m in meta:
            if not isinstance(m, AbstractMeta):
                raise TypeError(
                    f"meta items must be AbstractMeta subclasses; got {type(m)}")
        return list(meta)
    
    raise TypeError(
        f"meta must be None, AbstractMeta, or list of AbstractMeta; got {type(meta)}")

def _meta_to_julia_vector(meta_list: Union[List[AbstractMeta], None]) -> Any:
    """
    Convert a Python meta list to a Julia vector of AbstractMeta objects.
    
    Uses Julia's vect() function to construct the vector from individual objects.
    """
    from .bifrost_py import get_jl
    
    jl = get_jl()
    
    # Empty list case
    if meta_list is None or len(meta_list) == 0:
        return jl.eval("AbstractMeta[]")
    
    # Build each Julia meta object
    julia_meta_objs = []
    for m in meta_list:
        if isinstance(m, Nickname):
            jl_meta = jl.Bifrost.Nickname(m.label)
        elif isinstance(m, MCMadd):
            jl_dist = _distribution_to_julia(m.distribution)
            jl_meta = jl.Bifrost.MCMadd(jl.Symbol(m.symbol), jl_dist)
        elif isinstance(m, MCMmul):
            jl_dist = _distribution_to_julia(m.distribution)
            jl_meta = jl.Bifrost.MCMmul(jl.Symbol(m.symbol), jl_dist)
        else:
            raise TypeError(f"Unknown meta type: {type(m)}")
        
        julia_meta_objs.append(jl_meta)
    
    # Use Julia's vect() function to construct the vector from variadic args
    # vect unpacks the list into individual arguments
    return jl.Base.vect(*julia_meta_objs)


# ============================================================================
# Export
# ============================================================================

__all__ = [
    'AbstractMeta',
    'Nickname',
    'MCMadd',
    'MCMmul',
]