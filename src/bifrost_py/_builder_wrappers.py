"""
Wrappers around Julia builder functions that handle Python↔Julia meta conversion.

The builders (start_b, straight_b, etc.) are dynamically forwarded from Julia,
but they need to convert Python meta lists to Julia vectors before calling.
"""

from typing import Any, Optional, List
from .bifrost_py import get_jl
from ._meta import _meta_to_julia_vector, validate_meta


def _wrap_builder_call(jl_func, *args, meta=None, **kwargs):
    """
    Generic wrapper for Julia builder functions that handles meta conversion.
    
    Converts the Python meta list to a Julia vector before passing to Julia.
    """
    # Validate and convert meta
    if meta is not None:
        validated_meta = validate_meta(meta)
        jl_meta_vector = _meta_to_julia_vector(validated_meta)
        kwargs['meta'] = jl_meta_vector
    
    # Call the Julia function
    return jl_func(*args, **kwargs)


def start_b(spec, point=None, outgoing_tangent=None, outgoing_curvature=None, 
            spin_rate=None, meta=None):
    """
    Wrapper for Julia start! with meta handling.
    
    Parameters
    ----------
    spec : SubpathBuilder
        The builder to initialize.
    point : tuple, optional
        Start point (x, y, z). Default (0, 0, 0).
    outgoing_tangent : tuple, optional
        Start tangent direction. Default (0, 0, 1).
    outgoing_curvature : tuple, optional
        Start curvature vector. Default (0, 0, 0).
    spin_rate : float, callable, ':inherit', or None
        Material spin rate over the whole Subpath.
    meta : AbstractMeta or list of AbstractMeta, optional
        Annotations (Nickname, MCMadd, MCMmul).
    """
    jl = get_jl()
    kwargs = {}
    if point is not None:
        kwargs['point'] = point
    if outgoing_tangent is not None:
        kwargs['outgoing_tangent'] = outgoing_tangent
    if outgoing_curvature is not None:
        kwargs['outgoing_curvature'] = outgoing_curvature
    if spin_rate is not None:
        kwargs['spin_rate'] = spin_rate
    
    return _wrap_builder_call(jl.start_b, spec, meta=meta, **kwargs)


def straight_b(spec, length=None, twist_rate=None, meta=None):
    """
    Wrapper for Julia straight! with meta handling.
    
    Parameters
    ----------
    spec : SubpathBuilder
        The builder to append to.
    length : float
        Segment length (m). Negative length walks backward.
    twist_rate : float, callable, or None
        Mechanical twist rate (rad/m).
    meta : AbstractMeta or list of AbstractMeta, optional
        Annotations.
    """
    jl = get_jl()
    kwargs = {}
    if length is not None:
        kwargs['length'] = length
    if twist_rate is not None:
        kwargs['twist'] = twist_rate
    
    return _wrap_builder_call(jl.straight_b, spec, meta=meta, **kwargs)


def bend_b(spec, radius=None, angle=None, axis_angle=0.0, 
           twist_rate=None, meta=None):
    """
    Wrapper for Julia bend! with meta handling.
    
    Parameters
    ----------
    spec : SubpathBuilder
        The builder to append to.
    radius : float
        Bend radius (m).
    angle : float
        Total angle swept (rad).
    axis_angle : float
        Orientation in transverse plane (rad). Default 0.
    twist_rate : float, callable, or None
        Mechanical twist rate (rad/m).
    meta : AbstractMeta or list of AbstractMeta, optional
        Annotations.
    """
    jl = get_jl()
    kwargs = {}
    if radius is not None:
        kwargs['radius'] = radius
    if angle is not None:
        kwargs['angle'] = angle
    if axis_angle != 0.0:
        kwargs['axis_angle'] = axis_angle
    if twist_rate is not None:
        kwargs['twist'] = twist_rate
    
    return _wrap_builder_call(jl.bend_b, spec, meta=meta, **kwargs)


def helix_b(spec, radius=None, pitch=None, turns=None, axis_angle=0.0,
            twist_rate=None, meta=None):
    """
    Wrapper for Julia helix! with meta handling.
    
    Parameters
    ----------
    spec : SubpathBuilder
        The builder to append to.
    radius : float
        Helix radius (m).
    pitch : float
        Pitch per turn (m).
    turns : float
        Number of turns.
    axis_angle : float
        Orientation (rad). Default 0.
    twist_rate : float, callable, or None
        Mechanical twist rate (rad/m).
    meta : AbstractMeta or list of AbstractMeta, optional
        Annotations.
    """
    jl = get_jl()
    kwargs = {}
    if radius is not None:
        kwargs['radius'] = radius
    if pitch is not None:
        kwargs['pitch'] = pitch
    if turns is not None:
        kwargs['turns'] = turns
    if axis_angle != 0.0:
        kwargs['axis_angle'] = axis_angle
    if twist_rate is not None:
        kwargs['twist'] = twist_rate
    
    return _wrap_builder_call(jl.helix_b, spec, meta=meta, **kwargs)


def catenary_b(spec, a=None, length=None, axis_angle=0.0, 
               twist_rate=None, meta=None):
    """
    Wrapper for Julia catenary! with meta handling.
    
    Parameters
    ----------
    spec : SubpathBuilder
        The builder to append to.
    a : float
        Catenary parameter (m).
    length : float
        Arc length (m).
    axis_angle : float
        Orientation (rad). Default 0.
    twist_rate : float, callable, or None
        Mechanical twist rate (rad/m).
    meta : AbstractMeta or list of AbstractMeta, optional
        Annotations.
    """
    jl = get_jl()
    kwargs = {}
    if a is not None:
        kwargs['a'] = a
    if length is not None:
        kwargs['length'] = length
    if axis_angle != 0.0:
        kwargs['axis_angle'] = axis_angle
    if twist_rate is not None:
        kwargs['twist'] = twist_rate
    
    return _wrap_builder_call(jl.catenary_b, spec, meta=meta, **kwargs)


def seal_b(spec, extra=0.0, twist_rate=None, meta=None):
    """
    Wrapper for Julia seal! with meta handling.
    
    Seals the Subpath at the natural exit (no global target bend).
    
    Parameters
    ----------
    spec : SubpathBuilder
        The builder to seal.
    extra : float
        Optional straight lead-out length (m). Default 0.
    twist_rate : float, callable, or None
        Mechanical twist rate on the terminal connector (rad/m).
    meta : AbstractMeta or list of AbstractMeta, optional
        Annotations on the terminal connector.
    """
    jl = get_jl()
    kwargs = {}
    if extra != 0.0:
        kwargs['extra'] = extra
    if twist_rate is not None:
        kwargs['twist'] = twist_rate
    
    return _wrap_builder_call(jl.seal_b, spec, meta=meta, **kwargs)


__all__ = [
    'start_b',
    'straight_b',
    'bend_b',
    'helix_b',
    'catenary_b',
    'seal_b',
]