using Test
using LinearAlgebra
using Bifrost

# End-to-end convergence study for the adaptive step-doubling solver (issue #26
# task 2: "how small is small enough" for rtol/atol, and pin a recommended
# default). The study has two halves that together answer the question:
#
#   1. Realistic fibers are tolerance-insensitive. A weakly-birefringent
#      step-index fiber built from straights and constant-radius bends has a
#      piecewise-constant generator K(s); the exponential-midpoint integrator is
#      *exact* for constant K at any step size, so the Jones matrix is identical
#      (to round-off) whether rtol is 1e-4 or 1e-9. The default is therefore
#      already far smaller than such fibers need.
#
#   2. A demanding generator reveals the controller. To actually exercise the
#      tolerance machinery we drive the same propagator with a strong, smooth,
#      non-commuting K(s) (O(1) accumulated rotation, the structure used in
#      test_dgd.jl). There the error tracks rtol as designed, and we confirm the
#      package default (rtol = 1e-9) reaches a documented accuracy target.

const _SX = ComplexF64[0 1; 1 0]
const _SY = ComplexF64[0 -im; im 0]
const _SZ = ComplexF64[1 0; 0 -1]

# Demanding, smooth, non-commuting generator and its angular-frequency
# derivative. Coefficients are O(1) so the accumulated rotation over the domain
# is large and the integrator's local error genuinely depends on the tolerance.
_conv_K(s) = im * (2.0 * sin(3.0 * s)) .* _SX +
             im * (1.5 * cos(2.0 * s)) .* _SY +
             im * 0.8 .* _SZ
_conv_Kω(s) = im * (0.3 * sin(2.0 * s)) .* _SX +
              im * (0.2 * cos(1.5 * s)) .* _SZ
const _CONV_BREAKS = [0.0, 1.0, 2.5, 4.0]   # multi-interval: exercises piecewise

# Recommended-default accuracy targets on the demanding generator. At the package
# default (rtol = 1e-9) the solver reaches ~3e-7 on the Jones matrix and ~1e-7 on
# the DGD; the targets below sit just above those with margin. Update them
# intentionally if the solver or step controller changes.
const _CONV_J_TARGET = 1e-6
const _CONV_DGD_TARGET = 1e-6

# Realistic SMF-like fiber, mirroring test/human/demo-smallest.jl.
const _CONV_XS = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036), SilicaGermaniaGlass(0.0),
    8.2e-6, 125e-6,
)
const _CONV_T_REF = 297.15
const _CONV_λ = 1550e-9

function _conv_fiber()
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.5, meta = [Nickname("lead-in")])
    bend!(sb; radius = 0.05, angle = π / 2, meta = [Nickname("90 deg bend")])
    straight!(sb; length = 0.5, meta = [Nickname("lead-out")])
    seal!(sb)
    return Fiber(build(sb); cross_section = _CONV_XS, T_ref_K = _CONV_T_REF)
end

@testset "solver tolerance convergence" begin
    @testset "realistic fiber is tolerance-insensitive" begin
        # T-SIM-REGRESSION: a piecewise-constant fiber generator integrates
        # exactly, so loosening rtol by five orders of magnitude must not move
        # the Jones matrix. This documents that the default is more than small
        # enough for realistic fibers.
        fiber = _conv_fiber()
        J_loose, _ = propagate_fiber(
            fiber; λ_m = _CONV_λ, verbose = false,
            params = SolverParams(rtol = 1e-4, atol = 1e-7),
        )
        J_tight, _ = propagate_fiber(
            fiber; λ_m = _CONV_λ, verbose = false,
            params = SolverParams(rtol = 1e-9, atol = 1e-12),
        )
        @test phase_insensitive_error(J_tight, J_loose) < 1e-10
    end

    # Ground-truth reference for the demanding generator, tighter than the sweep.
    ref_params = SolverParams(rtol = 1e-13, atol = 1e-15)
    J_ref, _ = propagate_piecewise(
        _conv_K, _CONV_BREAKS; verbose = false, params = ref_params,
    )
    Js_ref, G_ref, _ = propagate_piecewise_sensitivity(
        _conv_K, _conv_Kω, _CONV_BREAKS; verbose = false, params = ref_params,
    )
    dgd_ref = output_dgd(Js_ref, G_ref)

    # Sweep from loose to tight. atol tracks rtol so the relative term dominates.
    rtols = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
    j_errs = Float64[]
    dgd_errs = Float64[]
    for rt in rtols
        p = SolverParams(rtol = rt, atol = rt * 1e-3)
        J, _ = propagate_piecewise(_conv_K, _CONV_BREAKS; verbose = false, params = p)
        Js, G, _ = propagate_piecewise_sensitivity(
            _conv_K, _conv_Kω, _CONV_BREAKS; verbose = false, params = p,
        )
        push!(j_errs, phase_insensitive_error(J_ref, J))
        push!(dgd_errs, abs(output_dgd(Js, G) - dgd_ref))
    end

    @testset "errors shrink as tolerance tightens" begin
        # T-SIM-REGRESSION: tightening rtol must improve accuracy. We require a
        # large end-to-end gain and forbid any single tightening step from making
        # the Jones error meaningfully worse. The additive 1e-13 slack absorbs
        # floor noise once the error nears the reference's own truncation level.
        @test j_errs[end] < j_errs[1] / 10
        @test dgd_errs[end] < dgd_errs[1] / 10
        for i in 1:(length(j_errs) - 1)
            @test j_errs[i + 1] <= 2 * j_errs[i] + 1e-13
        end
    end

    @testset "recommended default is small enough" begin
        # T-SIM-REGRESSION: the package default SolverParams() uses rtol = 1e-9.
        # Confirm that default already meets the documented accuracy targets on
        # the demanding generator, so users need not hand-tune tolerances.
        @test SolverParams().rtol == 1e-9
        idx_default = findfirst(==(1e-9), rtols)
        @test idx_default !== nothing
        @test j_errs[idx_default] < _CONV_J_TARGET
        @test dgd_errs[idx_default] < _CONV_DGD_TARGET
    end
end
