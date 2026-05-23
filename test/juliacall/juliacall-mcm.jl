module MCMDemo

using MonteCarloMeasurements
using Random
using Printf

if !isdefined(Main, :refractive_index)
    include(joinpath(@__DIR__, "..", "..", "src", "material-properties.jl"))
end

# Pure-silica refractive index at fixed wavelength, evaluated on M independent
# T_i ~ Normal(T_nom, T_sigma) Particles ensembles. Returns reduced scalars
# (pmean, pstd) so the result is trivially serializable and byte-comparable
# across runs.
function run(; N::Integer = 2000,
               λ::Real = 1550e-9,
               T_nom::Real = 293.0,
               T_sigma::Real = 10.0,
               M::Integer = 8,
               seed::Integer = 0xB1F205,
               threaded::Bool = false)
    Random.seed!(seed)
    # Build T inputs deterministically and serially, so RNG draws don't
    # depend on thread scheduling. Use the explicit N-particle constructor
    # so callers can vary ensemble size.
    Ts = [Particles(N, MonteCarloMeasurements.Normal(T_nom, T_sigma))
          for _ in 1:M]

    T_means = Vector{Float64}(undef, M)
    T_stds  = Vector{Float64}(undef, M)
    n_means = Vector{Float64}(undef, M)
    n_stds  = Vector{Float64}(undef, M)

    prev = MonteCarloMeasurements.unsafe_comparisons()
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        if threaded
            Threads.@threads for i in 1:M
                _eval!(i, Ts[i], λ, T_means, T_stds, n_means, n_stds)
            end
        else
            for i in 1:M
                _eval!(i, Ts[i], λ, T_means, T_stds, n_means, n_stds)
            end
        end
    finally
        MonteCarloMeasurements.unsafe_comparisons(prev)
    end

    return [(i = i,
             T_mean = T_means[i], T_std = T_stds[i],
             n_mean = n_means[i], n_std = n_stds[i]) for i in 1:M]
end

@inline function _eval!(i, T, λ, T_means, T_stds, n_means, n_stds)
    n = refractive_index(PURE_SILICA, λ, T)
    T_means[i] = pmean(T)
    T_stds[i]  = pstd(T)
    n_means[i] = pmean(n)
    n_stds[i]  = pstd(n)
    return nothing
end

# Stable CSV serialization. Header is fixed; values are %.17g so that any
# Float64 round-trips exactly. No trailing whitespace.
function write_table(path::AbstractString, rows)
    open(path, "w") do io
        println(io, "i,T_mean,T_std,n_mean,n_std")
        for r in rows
            @printf(io, "%d,%.17g,%.17g,%.17g,%.17g\n",
                    r.i, r.T_mean, r.T_std, r.n_mean, r.n_std)
        end
    end
    return path
end

nthreads() = Threads.nthreads()

end # module MCMDemo
