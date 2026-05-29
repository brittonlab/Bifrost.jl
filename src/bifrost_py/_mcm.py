"""
Monte Carlo Measurements (MCM) integration.

Provides Python-friendly Particles and StaticParticles wrappers that integrate
with Julia's MonteCarloMeasurements.jl for uncertainty quantification.
"""

from typing import Optional, Callable, Any, Union
import numpy as np
from scipy import stats as scipy_stats
from .bifrost_py import get_jl


class Particles:
    """
    Ensemble of Monte Carlo samples (alias for MonteCarloMeasurements.Particles).
    
    Represents an uncertain quantity as a distribution of values. Can be used
    directly in Bifrost functions; each function call produces a distribution
    of outputs.
    
    Parameters
    ----------
    n : int
        Number of samples/particles.
    distribution : scipy.stats distribution
        Univariate continuous distribution (e.g., norm, uniform, lognorm).
        Must have .rvs(size=n) method.
    seed : int, optional
        Random seed for reproducibility.
    
    Attributes
    ----------
    particles : np.ndarray
        Array of n samples (readonly after creation).
    mean : float
        Mean of the samples.
    std : float
        Standard deviation of the samples.
    
    Examples
    --------
    Temperature uncertainty (normal distribution)::
    
        from scipy.stats import norm
        import bifrost as bf
        
        T_K = bf.mcm.Particles(100, norm(293.15, 5))  # 20°C ± 5°C
        
        # Use in fiber setup
        fiber = bf.Fiber(path, cross_section=xs, temperature_k=T_K)
        
        # Propagate and get output distribution
        J, stats = bf.propagate_fiber(fiber, wavelength_m=1550e-9)
        
        # J[0,0] now contains a distribution
        print(J[0, 0].mean)
        print(J[0, 0].std)
    
    Bend radius uncertainty (lognormal)::
    
        from scipy.stats import lognorm
        
        R = bf.mcm.Particles(50, lognorm(0.1, scale=0.01))  # Log-normal radius
        
        path_builder = bf.PathSpecBuilder()
        path_builder.add_bend(radius_m=R, angle_rad=np.pi/2)
        # ...
    """
    
    def __init__(
        self,
        n: int,
        distribution: Any,
        seed: Optional[int] = None,
    ):
        if seed is not None:
            np.random.seed(seed)
        
        self.n = n
        self.particles = distribution.rvs(size=n)
        self.distribution = distribution
    
    @property
    def mean(self) -> float:
        """Mean of the particles."""
        return float(np.mean(self.particles))
    
    @property
    def std(self) -> float:
        """Standard deviation of the particles."""
        return float(np.std(self.particles))
    
    @property
    def _julia_type(self) -> Any:
        """Convert to Julia Particles object."""
        jl = get_jl()
        # Use juliacall to pass to Julia's MCM.Particles
        return jl.MonteCarloMeasurements.Particles(self.particles)
    
    def __repr__(self) -> str:
        return f"Particles(n={self.n}, mean={self.mean:.4f}, std={self.std:.4f})"


class StaticParticles:
    """
    Static ensemble: creates particles once and reuses them (no randomness per call).
    
    Useful when you want fixed, reproducible samples that don't change if
    a computation is repeated.
    
    Parameters
    ----------
    n : int
        Number of particles.
    distribution : scipy.stats distribution
        Univariate continuous distribution.
    seed : int, optional
        Random seed.
    
    Examples
    --------
    >>> T_K = bf.mcm.StaticParticles(20, norm(293.15, 5))
    >>> fiber1 = bf.Fiber(path, temperature_k=T_K)
    >>> fiber2 = bf.Fiber(path, temperature_k=T_K)
    >>> # fiber1 and fiber2 have identical temperature ensembles
    """
    
    def __init__(
        self,
        n: int,
        distribution: Any,
        seed: Optional[int] = None,
    ):
        if seed is not None:
            np.random.seed(seed)
        
        self.n = n
        self.particles = distribution.rvs(size=n)
        self.distribution = distribution
    
    @property
    def mean(self) -> float:
        return float(np.mean(self.particles))
    
    @property
    def std(self) -> float:
        return float(np.std(self.particles))
    
    @property
    def _julia_type(self) -> Any:
        """Convert to Julia StaticParticles object."""
        jl = get_jl()
        return jl.MonteCarloMeasurements.StaticParticles(self.particles)
    
    def __repr__(self) -> str:
        return f"StaticParticles(n={self.n}, mean={self.mean:.4f}, std={self.std:.4f})"