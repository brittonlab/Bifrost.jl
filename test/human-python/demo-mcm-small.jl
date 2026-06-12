"""
This is a VERY simple demo of the MCM capabilities of the code. 
We don't even use the MCMAdd metas yet; right now the temperature is
varied over the *entire* fiber, which is a straight + bend + straight
combo. The MCMAdd meta is required for the task of varying the temperature
of a single element; we'll do that next.
"""

using Bifrost
using Distributions
using Printf
using MonteCarloMeasurements

# =============================
# Setup: fiber with uncertain temperature
# =============================

println("Starting...")

xs = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),
    8.2e-6,
    125.0e-6,
)

# Build a simple path
sb = SubpathBuilder(); start!(sb)
straight!(sb; length = 0.5, meta = [Nickname("lead-in")])
bend!(sb;     radius = 0.05, angle = π / 2, meta = [Nickname("90 deg bend")])
straight!(sb; length = 0.5, meta = [Nickname("lead-out")])
seal!(sb)
path = build(sb)

# =============================
# MCM: Temperature uncertainty
# =============================

# Temperature: 24°C ± 5°C (normal distribution, 50 samples)
const N_SAMPLES = 100

T_dist = Normal(297.15, 5.0)  # mean, std
T_K_samples = Particles(N_SAMPLES, T_dist)

# Fiber with uncertain temperature
fiber_mcm = Fiber(path; cross_section = xs, T_ref_K = T_K_samples)

# =============================
# Propagate
# =============================

println("Propagating fiber with uncertain temperature ($(N_SAMPLES) samples)...")
J_particles, stats = propagate_fiber(
    fiber_mcm;
    λ_m = 1550e-9,
    verbose = false,
)

println("Done. J shape: $(size(J_particles))")

# =============================
# Convert Jones matrix to rotation angle
# =============================

function jones_to_rotation_angle(J::AbstractMatrix)
    """
    Extract rotation angle (rad) from a 2×2 Jones matrix.
    
    For a unitary matrix J = exp(i·θ·n̂·σ), the rotation angle θ
    satisfies: tr(J) = 2·cos(θ)
    """
    tr_J = J[1,1] + J[2,2]
    # Clamp to avoid numerical issues with arccos
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

for k in 1:n_samples
    # Extract k-th sample from each matrix element
    J_k = ComplexF64[
        real(J_particles[1,1]).particles[k] + im*imag(J_particles[1,1]).particles[k]   real(J_particles[1,2]).particles[k] + im*imag(J_particles[1,2]).particles[k];
        real(J_particles[2,1]).particles[k] + im*imag(J_particles[2,1]).particles[k]   real(J_particles[2,2]).particles[k] + im*imag(J_particles[2,2]).particles[k]
    ]
    θ_k = jones_to_rotation_angle(J_k)
    push!(angles, θ_k)
    #k <= 3 ? display(J_k) : nothing
end

# =============================
# Results
# =============================

@printf "\nResults:\n"
@printf "Temperature range: %.3f K - %.3f K \n" minimum(T_K_samples.particles) maximum(T_K_samples.particles)
@printf "Rotation angles (°): \n"
@printf "  Mean ± Std:  %.3f ± %.4f \n" mean(angles)*180/π std(angles)*180/π
@printf "  Min - Max:   %.3f - %.3f\n" minimum(angles)*180/π maximum(angles)*180/π
# println("\nFirst 5 samples:")
# for k in 1:min(5, n_samples)
#     T_k = T_K_samples.particles[k]
#     θ_k = angles[k]
#     @printf "  Sample %d: T=%.2f°C, θ=%.3f rad\n" k (T_k-273.15) θ_k*180/π
# end