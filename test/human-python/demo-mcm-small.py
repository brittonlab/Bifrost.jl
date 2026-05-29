# demo_mcm_minimal.py
"""Minimal MCM demo: temperature-dependent Jones matrix rotation angle."""

import bifrost_py as bf
import numpy as np
from scipy.stats import norm


def jones_to_rotation_angle(J):
    """
    Extract rotation angle (rad) from a 2×2 Jones matrix.
    
    For a unitary matrix J = exp(i·θ·n̂·σ), the rotation angle θ
    satisfies: tr(J) = 2·cos(θ)
    
    Parameters
    ----------
    J : np.ndarray
        2×2 complex Jones matrix.
    
    Returns
    -------
    float
        Rotation angle in radians.
    """
    tr_J = np.trace(J)
    cos_theta = np.real(tr_J) / 2
    # Clamp to avoid numerical issues with arccos
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return np.arccos(cos_theta)


# =============================
# Setup: fiber with uncertain temperature
# =============================

print("Setting up fiber...")

xs = bf.StepIndexCrossSection(
    bf.SilicaGermaniaGlass(0.036),
    bf.SilicaGermaniaGlass(0.0),
    8.2e-6,
    125.0e-6,
)

# Build a simple path: straight 1 meter
spec = bf.PathSpecBuilder()
bf.straight_b(spec, length = 0.5)
bf.bend_b(spec, radius = 0.01, angle = bf.pi / 2)
bf.straight_b(spec, length = 0.5)
path = bf.build(spec)

# Reference temperature
T_ref_K = 297.15  # 24°C

# =============================
# MCM: Temperature uncertainty
# =============================

print("Creating Monte Carlo ensemble...")

# Temperature: 24°C ± 5°C (normal distribution, 50 samples)
T_dist = norm(loc=297.15, scale=5.0)  # mean, std
T_K_particles = bf.mcm.StaticParticles(n=50, distribution=T_dist, seed=42)

print(f"  Temperature: {T_K_particles.mean:.2f} ± {T_K_particles.std:.2f} K")
print(f"  Samples: {T_K_particles.n}")
print(T_K_particles._julia_type)

# Fiber with uncertain temperature
fiber_mcm = bf.Fiber(path, cross_section=xs, T_ref_K=T_K_particles)

# =============================
# Propagate
# =============================

print("Propagating fiber with uncertain temperature...")
J_particles, stats = bf.propagate_fiber(
    fiber_mcm,
    lambda_m=1550e-9,
    rtol=1e-9,
    verbose=False,
)

print(f"Done. J is a ParticlesMatrix with shape {J_particles.particles.shape}")
print(f"  (2, 2, {J_particles.n_samples}) samples\n")

# =============================
# Extract per-sample rotation angles
# =============================

print("Extracting rotation angles per sample...")

angles = []
for k in range(J_particles.n_samples):
    # Extract k-th sample: shape (2, 2)
    J_k = J_particles.particles[:, :, k]
    theta_k = jones_to_rotation_angle(J_k)
    angles.append(theta_k)

angles = np.array(angles)

# =============================
# Results
# =============================

print("\nResults:")
print(f"  Temperature range: {T_K_particles.particles.min():.2f}–{T_K_particles.particles.max():.2f} K")
print(f"  Rotation angle (rad):")
print(f"    Mean:  {np.mean(angles):.6f}")
print(f"    Std:   {np.std(angles):.6f}")
print(f"    Min:   {np.min(angles):.6f}")
print(f"    Max:   {np.max(angles):.6f}")

print("\nFirst 5 samples:")
for k in range(min(5, J_particles.n_samples)):
    T_k = T_K_particles.particles[k]
    theta_k = angles[k]
    print(f"  Sample {k+1}: T={T_k-273.15:.2f}°C, θ={theta_k:.6f} rad")