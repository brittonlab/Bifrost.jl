"""
During testing, I noticed an unusual behavior that this script is designed
to replicate. Consider a simple straight - bend - straight path for a fiber.
I noticed that if the temperature of the entire fiber was varied, the variation
in total polarization rotation angle change was much larger than if only the
temperature of the bend was changed. But this shouldn't be the case, because
only the bend has birefringence; the straight segments should contribute nothing
at all to the Jones matrix of the fiber.

This script builds both and tests both, converts the Jones matrices to rotation
angles, and then prints out the mean ± std of both cases. 
The straight segments should have NO birefringence, so it should not matter how
we parametrize the Jones matrices. 

To help reduce variation between tests, we even use the same sample of temperatures
for both tests.
"""

using Bifrost
using Distributions
using Printf
using MonteCarloMeasurements

const N_SAMPLES = 100
const λ_0 = 1550e-9

# Reference temperature: 24°C (fixed for the whole fiber)
T_ref_K = 297.15
# Temperature: 24°C ± 5°C (normal distribution, 50 samples)
T_dist = Normal(T_ref_K, 5.0)  # mean, std
T_K_samples0 = Particles(N_SAMPLES, T_dist)
T_K_samples1 = T_K_samples0 - T_ref_K
println("Checking temperature sample values:")
@printf "For T_K_samples0:  %.2f ± %.2f K\n" mean(T_K_samples0.particles) std(T_K_samples0.particles)
@printf "For T_K_samples1:  %.2f ± %.2f K\n" mean(T_K_samples1.particles) std(T_K_samples1.particles)

# Needed utility
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


##### Cross section
xs = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),
    8.2e-6,
    125.0e-6,
    manufacturer = "Corning",
    model_number = "SMF-like",
)


println("Now beginning whole-fiber temperature variation...")
##### Set up and build the fiber whose temperature will be wholly varied
spec0 = PathSpecBuilder()
straight!(spec0; length = 0.5, meta = [Nickname("lead-in")])
bend!(spec0; radius = 0.01, angle = pi / 2, meta = [Nickname("90 deg bend")])
straight!(spec0; length = 0.5, meta = [Nickname("lead-out")])
path0 = build(spec0)

# Fiber with uncertain temperature
fiber_mcm0 = Fiber(path0; cross_section = xs, T_ref_K = T_K_samples0)

println("    Propagating fiber with uncertain temperature ($(N_SAMPLES) samples)...")
J_particles0, stats0 = propagate_fiber(
    fiber_mcm0;
    λ_m = λ_0,
    rtol = 1e-9,
    verbose = false,
)

println("    Extracting rotation angles per sample...")
n_samples = length(real(J_particles0[1, 1]).particles)
angles0 = Float64[]
on_diags0 = ComplexF64[]

for k in 1:n_samples
    # Extract k-th sample from each matrix element
    J_k = ComplexF64[
        real(J_particles0[1,1]).particles[k] + im*imag(J_particles0[1,1]).particles[k]   real(J_particles0[1,2]).particles[k] + im*imag(J_particles0[1,2]).particles[k];
        real(J_particles0[2,1]).particles[k] + im*imag(J_particles0[2,1]).particles[k]   real(J_particles0[2,2]).particles[k] + im*imag(J_particles0[2,2]).particles[k]
    ]
    θ_k = jones_to_rotation_angle(J_k)
    push!(angles0, θ_k)
    push!(on_diags0, J_k[1,1])
end
println("    Done.")


##### Now do the bend-only fiber with MCMAdd meta
println("Now beginning bend-only temperature variation...")
# Apply segment-level MCM perturbation via meta
# The bend segment carries MCMadd(:T_K, ΔT_K_particles) in its meta,
# which modify() will interpret as "add ΔT_K to the local temperature"
spec_with_meta = PathSpecBuilder()
straight!(spec_with_meta; length = 0.5, meta = [Nickname("lead-in")])
bend!(spec_with_meta; radius = 0.01, angle = π/2,
      meta = [Nickname("sensitive bend"), MCMadd(:T_K, T_K_samples1)])
straight!(spec_with_meta; length = 0.5, meta = [Nickname("lead-out")])

path_with_meta = build(spec_with_meta)
fiber_with_meta = Fiber(path_with_meta; cross_section = xs, T_ref_K = T_ref_K)

# Modify: Apply MCM perturbations to path
println("    Applying MCM perturbations via modify()...")
modified_path = modify(fiber_with_meta)
fiber_modified = Fiber(modified_path; cross_section = xs, T_ref_K = T_ref_K)

# Propagate
println("    Propagating with segment-level temperature variation...")
J_particles1, stats1 = propagate_fiber(
    fiber_modified;
    λ_m = λ_0,
    rtol = 1e-9,
    verbose = false,
)

# Extract angles
println("    Extracting rotation angles per sample...")
n_samples = length(real(J_particles1[1, 1]).particles)
angles1 = Float64[]
on_diags1 = ComplexF64[]

for k in 1:n_samples
    # Extract k-th sample
    J_k = ComplexF64[
        real(J_particles1[1,1]).particles[k] + im*imag(J_particles1[1,1]).particles[k]   real(J_particles1[1,2]).particles[k] + im*imag(J_particles1[1,2]).particles[k];
        real(J_particles1[2,1]).particles[k] + im*imag(J_particles1[2,1]).particles[k]   real(J_particles1[2,2]).particles[k] + im*imag(J_particles1[2,2]).particles[k]
    ]
    θ_k = jones_to_rotation_angle(J_k)
    push!(angles1, θ_k)
    push!(on_diags1, J_k[1,1])
end
println("    Done.")


##### Printing of final results
println("Rotation angles (°) for the two cases:")
println("Whole-fiber temperature varied:")
@printf "  Mean ± Std:  %.5f ± %.6f \n" mean(angles0)*180/π std(angles0)*180/π
println("Bend-only temperature varied:")
@printf "  Mean ± Std:  %.5f ± %.6f \n" mean(angles1)*180/π std(angles1)*180/π

# This currently returns:
#
# Rotation angles (°) for the two cases:
# Whole-fiber temperature varied:
#   Mean ± Std:  9.50209 ± 0.008657 
# Bend-only temperature varied:
#   Mean ± Std:  9.50256 ± 0.000051
#
# The fact that the SDs are very different (and the means are different at all)
# is an indication of an underlying physics failure.

println(" ")
println(" ")

# Two possible failure points: MCMAdd meta interpretation and generator construction.
#
# The generators should be easy to check. So let's start there.
fiber_base = Fiber(path0; cross_section = xs, T_ref_K = T_ref_K)

# The propagate_fiber() call that should work:
# J_base, stats_base = propagate_fiber(
#     fiber_base;
#     λ_m = 1550e-9,
#     rtol = 1e-9,
#     verbose = false,
# )
# But we won't do that...
K_base = generator_K(fiber_base, fiber_base.cross_section, λ_0)
println("Base generator checks:")
println(K_base(0.1))
println(K_base(0.5005))
println(K_base(0.9))
result = K_base(0.1)
println("Type of result[1,1]: $(typeof(result[1,1]))")
println(result[1,1])
result = K_base(0.5005)
println("Type of result[1,1]: $(typeof(result[1,1]))")
println(result[1,1])
#
# Running this:
# println(K_base(0.1))
# println(K_base(0.5005))
# println(K_base(0.7))
# This returns:
# ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]
# ComplexF64[0.0 - 10.558396096235342im -0.0 - 0.0im; 0.0 - 0.0im 0.0 + 10.558396096235342im]
# ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]
#
# All of this is reasonable. It might seem weird that K = [0], but the generator is 
# K ∝ ϵ - n_avg^2. So that makes sense.
# Let's try it with the MCM-modified fiber:
K_MCM = generator_K(fiber_modified, fiber_modified.cross_section, λ_0)
println("MCM'd generator checks:")
println(K_MCM(0.1))
println(K_MCM(0.5005))
println(K_MCM(0.9))
result = K_MCM(0.1)
println("Type of result[1,1]: $(typeof(result[1,1]))")
println(result[1,1])
result = K_MCM(0.5005)
println("Type of result[1,1]: $(typeof(result[1,1]))")
println(result[1,1])

# This returns the fact that, for the bend-only gend generator, we get a
# Complex{Particles{Float64, 100}} return type.
# Are the off-diagonals varying the same? Let's see.
@printf "Whole fiber: %.6f + %.6f i\n" std(real.(on_diags0)) std(imag.(on_diags0))
@printf "Bend only:   %.6f + %.6f i\n" std(real.(on_diags1)) std(imag.(on_diags1))