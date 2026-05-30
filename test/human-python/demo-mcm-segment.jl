# demo_mcm_segment_minimal.jl
"""
Minimal MCM demo with segment-level temperature perturbation.

Instead of applying temperature uncertainty to the entire fiber,
we apply MCMadd(:T_K, ΔT) to a single bend segment via meta.
"""

using Bifrost
using Distributions
using MonteCarloMeasurements
using Printf

const N_SAMPLES = 100

# =============================
# Setup: fiber with multi-segment path
# =============================

xs = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),
    8.2e-6,
    125.0e-6,
    manufacturer = "Corning",
    model_number = "SMF-like",
)

# Reference temperature: 24°C (fixed for the whole fiber)
T_ref_K = 297.15

# =============================
# MCM: Segment-level temperature perturbation
# =============================

# Temperature offset: ΔT ~ N(0, 5 K)  [±5°C variation around reference]
# This is passed to the bend via MCMadd(:T_K, ΔT_K)
ΔT_dist = Normal(0.0, 5.0)  # mean offset, std
ΔT_K_particles = Particles(N_SAMPLES, ΔT_dist)

println("Temperature perturbation ensemble:")
@printf "  Mean ΔT:  %.2f K\n" mean(ΔT_K_particles.particles)
@printf "  Std ΔT:   %.2f K\n" std(ΔT_K_particles.particles)
@printf "  Min ΔT:   %.2f K\n" minimum(ΔT_K_particles.particles)
@printf "  Max ΔT:   %.2f K\n" maximum(ΔT_K_particles.particles)

# Apply segment-level MCM perturbation via meta
# The bend segment carries MCMadd(:T_K, ΔT_K_particles) in its meta,
# which modify() will interpret as "add ΔT_K to the local temperature"
spec_with_meta = PathSpecBuilder()
straight!(spec_with_meta; length = 0.5, meta = [Nickname("lead-in")])
bend!(spec_with_meta; radius = 0.01, angle = π/2,
      meta = [Nickname("sensitive bend"), MCMadd(:T_K, ΔT_K_particles)])
straight!(spec_with_meta; length = 0.5, meta = [Nickname("lead-out")])

path_with_meta = build(spec_with_meta)
fiber_with_meta = Fiber(path_with_meta; cross_section = xs, T_ref_K = T_ref_K)

# =============================
# Modify: apply MCM perturbations to the path
# =============================

println("\nApplying MCM perturbations via modify()...")
modified_path = modify(fiber_with_meta)
fiber_modified = Fiber(modified_path; cross_section = xs, T_ref_K = T_ref_K)

# =============================
# Propagate
# =============================

println("Propagating with segment-level temperature variation...")
J_particles, stats = propagate_fiber(
    fiber_modified;
    λ_m = 1550e-9,
    rtol = 1e-9,
    verbose = false,
)

println("Done.\n")

# =============================
# Convert Jones matrix to rotation angle
# =============================

function jones_to_rotation_angle(J)
    tr_J = J[1,1] + J[2,2]
    cos_theta = real(tr_J) / 2
    cos_theta = clamp(cos_theta, -1.0, 1.0)
    return acos(cos_theta)
end

# =============================
# Extract per-sample angles
# =============================

println("Extracting rotation angles per sample...")
n_samples = length(real(J_particles[1, 1]).particles)
angles = Float64[]
T_eff = Float64[]  # Effective temperature for each sample

for k in 1:n_samples
    # Extract k-th sample
    J_k = ComplexF64[
        real(J_particles[1,1]).particles[k] + im*imag(J_particles[1,1]).particles[k]   real(J_particles[1,2]).particles[k] + im*imag(J_particles[1,2]).particles[k];
        real(J_particles[2,1]).particles[k] + im*imag(J_particles[2,1]).particles[k]   real(J_particles[2,2]).particles[k] + im*imag(J_particles[2,2]).particles[k]
    ]
    θ_k = jones_to_rotation_angle(J_k)
    push!(angles, θ_k)
    push!(T_eff, T_ref_K + ΔT_K_particles.particles[k])
end

# =============================
# Results
# =============================

@printf "\nResults (segment-level MCM perturbation):\n"
@printf "  Rotation angle (°):\n"
@printf "    Mean:  %.6f\n" mean(angles)*180/π
@printf "    Std:   %.6f\n" std(angles)*180/π
@printf "    Min:   %.6f\n" minimum(angles)*180/π
@printf "    Max:   %.6f\n" maximum(angles)*180/π

@printf "\nFirst 5 samples:\n"
for k in 1:min(5, n_samples)
    @printf "  Sample %d: T_eff=%.2f K (ΔT=%+.2f K), θ=%.6f°\n" k T_eff[k] ΔT_K_particles.particles[k] angles[k]*180/π
end