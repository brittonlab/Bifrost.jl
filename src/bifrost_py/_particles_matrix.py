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
    """
    
    def __init__(self, jl_particles_matrix: Any):
        """
        Initialize from Julia 2x2 matrix of Complex{Particles}.
        
        Parameters
        ----------
        jl_particles_matrix : Any
            Julia matrix where each element is Complex{Particles}.
        """
        self._jl_mat = jl_particles_matrix
        self.shape = (2, 2)
        
        # Extract samples to (2, 2, n) numpy array
        # For Complex{Particles}, extract from real part to get n_samples
        first_element_real = jl_particles_matrix[0, 0].real
        n_samples = len(first_element_real.particles)
        
        self._particles_array = np.zeros(
            (2, 2, n_samples),
            dtype=np.complex128
        )
        
        for i in range(2):
            for j in range(2):
                element = jl_particles_matrix[i, j]  # This is Complex{Particles}
                # Extract real and imaginary parts (each is Particles)
                real_part = element.real.particles
                imag_part = element.imag.particles
                self._particles_array[i, j, :] = (
                    np.asarray(real_part, dtype=np.float64) +
                    1j * np.asarray(imag_part, dtype=np.float64)
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