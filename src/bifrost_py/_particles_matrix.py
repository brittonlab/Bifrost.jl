"""
ParticlesMatrix: wrapper for 2x2 matrices of MCM Particles.

Returned by propagate_fiber when inputs contain uncertain parameters.
"""

from typing import Any
import numpy as np


class ParticlesMatrix:
    """
    2×2 matrix of samples from Monte Carlo ensemble propagation.
    
    Provides convenient access to:
    - Raw samples: .particles (shape (2, 2, n_samples))
    - Statistical summaries: .mean, .std (shape (2, 2))
    - Individual element samples: pm[i, j] → (n_samples,)
    
    Examples
    --------
    >>> J, stats = bf.propagate_fiber(fiber_ensemble, ...)
    >>> print(J.particles.shape)  # (2, 2, 100)
    >>> print(J.mean)  # (2, 2) mean matrix
    >>> j00 = J.particles[0, 0, :]  # 100 samples of J[0,0]
    """
    
    def __init__(self, jl_particles_matrix: Any):
        """
        Initialize from Julia 2x2 matrix of Particles.
        
        Parameters
        ----------
        jl_particles_matrix : Any
            Julia matrix where each element is a Particles object.
        """
        self._jl_mat = jl_particles_matrix
        self.shape = (2, 2)
        
        # Extract samples to (2, 2, n) numpy array
        # Get n from first element
        first_element = jl_particles_matrix[0, 0]
        n_samples = len(first_element.particles)
        
        self._particles_array = np.zeros(
            (2, 2, n_samples),
            dtype=np.complex128
        )
        
        for i in range(2):
            for j in range(2):
                element_particles = jl_particles_matrix[i, j]
                self._particles_array[i, j, :] = np.asarray(
                    element_particles.particles,
                    dtype=np.complex128
                )
        
        # Compute statistics
        self._mean_array = np.mean(self._particles_array, axis=2)
        self._std_array = np.std(self._particles_array, axis=2)
    
    @property
    def particles(self) -> np.ndarray:
        """
        Raw samples: shape (2, 2, n_samples), dtype complex128.
        
        Access via `J.particles[i, j, :]` for individual element samples.
        """
        return self._particles_array
    
    @property
    def mean(self) -> np.ndarray:
        """
        Mean over samples: shape (2, 2).
        
        Equivalent to `J.particles.mean(axis=2)`.
        """
        return self._mean_array.copy()
    
    @property
    def std(self) -> np.ndarray:
        """
        Standard deviation over samples: shape (2, 2).
        
        Equivalent to `J.particles.std(axis=2)`.
        """
        return self._std_array.copy()
    
    @property
    def n_samples(self) -> int:
        """Number of Monte Carlo samples."""
        return self._particles_array.shape[2]
    
    def __repr__(self) -> str:
        return (
            f"ParticlesMatrix(shape=(2,2), n_samples={self.n_samples}, "
            f"mean=\n{self.mean})"
        )
    
    def __getitem__(self, key):
        """
        Index into the matrix or particles.
        
        Examples
        --------
        >>> J[0, 0]  # First element samples: shape (n_samples,)
        >>> J[0, 0, 5]  # 5th sample of J[0,0]: scalar
        """
        return self._particles_array[key]