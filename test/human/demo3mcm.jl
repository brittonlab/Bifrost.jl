# =====================================================================
# demo3mcm.jl — MCM temperature-dependent PTF demos
# =====================================================================
#
# Two demonstrations of temperature sensitivity on a multi-segment fiber:
#
#   * demo_mcm_temperature_ptf        — StaticParticles(50), single propagation,
#                                       Stokes parameters and polarisation angle vs T
#   * demo_mcm_temperature_ptf_scatter — Particles(2000), single propagation,
#                                       Poincaré equatorial scatter
#
# `demo3mcm_all()` runs both and writes `output/demo-index.html`.
#
# This file expects to be `include`d after demo2.jl is in scope (it
# reuses `_sample_segment_xyz` and the path-builder API).

if !isdefined(Main, :_sample_segment_xyz)
    include(joinpath(@__DIR__, "demo2.jl"))
end

# =====================================================================
# MCM demo — temperature-dependent PTF on a multi-segment fiber
# =====================================================================
#
# Illustrates end-to-end use of MCMadd(:T_K, …) on a helix segment.
# Segment structure:
#   straight 5 m
#   helix D=0.05 m (R=0.025 m), 10000 turns, pitch=0.05 m  ← MCMadd(:T_K, ΔT) here
#   straight 5 m
#   helix D=0.05 m (R=0.025 m), 10000 turns, pitch=0.05 m
#   straight 5 m
#
# Why is temperature sensitivity visible here?
# ─────────────────────────────────────────────
# Bending birefringence: Δβ ∝ 1/R².  At R=0.025 m, Δβ ≈ 3.38 rad/m —
# roughly 1500× larger than at R=5 m.  After 10000 turns the accumulated
# retardation |Δβ|·L ≈ 1775·2π rad.  A small change in Δβ (via temperature)
# shifts this by many radians, which rotates the output polarisation state.
#
# Polarisation crossing / mid-fringe operating point
# ───────────────────────────────────────────────────
# The Jones matrix of a linear retarder with retardation Γ(T) = Δβ(T)·L
# maps the Poincaré sphere in a way that maximises PTF sensitivity when the
# operating point is at a "mid-fringe" condition: mod(Γ, 2π) = π.  There
# the output state is maximally far from both the "crossing" (Γ mod 2π = 0)
# and its complement (2π), and the output Stokes vector moves at its fastest
# rate per degree of temperature.  By contrast at the crossing itself
# (Γ mod 2π = 0) the PTF is locally flat in temperature: dPTF/dT = 0.
#
# The turns count below (10001.892…) is chosen so that at T_nom = 30°C the
# first helix sits precisely at mid-fringe.  Over ΔT = ±5°C the accumulated
# retardation swings by ≈ ±0.75π rad, producing a large, observable rotation
# of the output state on the Poincaré sphere.
#
# Two plots are written:
#   1. mcm-temperature-ptf.html  — Stokes (S1,S2,S3) and DLP vs temperature
#      from a StaticParticles(50) run (T_nom=30°C ± 5°C).
#   2. mcm-temperature-ptf-scatter.html — DLP vs temperature from a single
#      MCM Particles run (T_nom=30°C ± 5°C), showing the ensemble scatter.

const _MCM_DEMO_XS = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),
    8.2e-6,
    125e-6,
)
const _MCM_DEMO_T_REF_K = 303.15   # 30°C reference temperature (= design point)
const _MCM_DEMO_λ_M     = 1550e-9
const _MCM_DEMO_T_NOM_C = 30.0     # nominal operating temperature (°C)
const _MCM_DEMO_T_SIG_C = 5.0      # ±5°C uncertainty (1σ)

# turns chosen so Δβ(30°C)·L = odd multiple of π (mid-fringe):
# k=887 → (2·887+1)·π / (|Δβ(30°C)| · arc_per_turn) ≈ 10001.892
const _MCM_DEMO_HELIX1_TURNS = 10001.892069208387

# Build the 5-segment fiber spec.  `ΔT_K` is the temperature offset applied
# to the first helix via MCMadd(:T_K, ΔT_K); pass 0.0 for the baseline.
function _mcm_demo_fiber(ΔT_K)
    spec = SubpathBuilder()
    start!(spec)
    # Sinusoidal spinning: amplitude 1 rad/m, period 10 m, starts here and runs
    # to end of fiber (no subsequent Spinning annotation resets it).
    straight!(spec; length = 5.0,
              meta = AbstractMeta[
                  Nickname("lead-in"),
                  Spinning(; rate = s -> sin(2π * s / 100.0)),
              ])
    helix!(spec; radius = 0.025, pitch = 0.05, turns = _MCM_DEMO_HELIX1_TURNS,
           axis_angle = 0.0,
           meta = AbstractMeta[Nickname("temperature-sensitive helix"),
                               MCMadd(:T_K, ΔT_K)])
    straight!(spec; length = 5.0, meta = [Nickname("spacer")])
    helix!(spec; radius = 0.025, pitch = 0.05, turns = 10000.0,
           axis_angle = 0.0, meta = [Nickname("reference helix")])
    straight!(spec; length = 5.0, meta = [Nickname("lead-out")])
    seal!(spec)
    # Fiber(builder) applies the :T_K thermal scaling (via the cladding CTE at
    # T_ref) during the single build; no separate modify step.
    return Fiber(spec; cross_section = _MCM_DEMO_XS, T_ref_K = _MCM_DEMO_T_REF_K)
end

# Apply a 2×2 Jones matrix to a normalised input state [1,0] and return
# normalised (s1,s2,s3).  Works for scalar ComplexF64 entries only.
function _jones_to_stokes(J::AbstractMatrix{ComplexF64})
    ψ = J * ComplexF64[1.0, 0.0]
    ψ ./= sqrt(abs2(ψ[1]) + abs2(ψ[2]))
    s0 = real(abs2(ψ[1]) + abs2(ψ[2]))
    s1 = real(abs2(ψ[1]) - abs2(ψ[2])) / s0
    s2 = 2 * real(ψ[1] * conj(ψ[2])) / s0
    s3 = -2 * imag(ψ[1] * conj(ψ[2])) / s0
    dlp = hypot(s1, s2)
    return (s1 = s1, s2 = s2, s3 = s3, dlp = dlp)
end

"""
    demo_mcm_temperature_ptf(; output_dir = …)

StaticParticles(50) temperature run: S1, S2, S3, and degree of linear
polarization (DLP) vs temperature for the MCM demo fiber.
Single propagation with T ~ N(30°C, 5°C), 50 particles.
"""
function demo_mcm_temperature_ptf(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    desc = "MCM temperature demo: DLP and Stokes parameters vs temperature " *
           "(StaticParticles(50), T = 30°C ± 5°C (1σ), single propagation) for a 5-segment fiber " *
           "with MCMadd(:T_K, ΔT) on the first large helix."

    MonteCarloMeasurements.unsafe_comparisons(true)
    T_C_particles  = StaticParticles(50, Normal(Float64(_MCM_DEMO_T_NOM_C), Float64(_MCM_DEMO_T_SIG_C)))
    T_K_particles  = T_C_particles + 273.15
    ΔT_K_particles = T_K_particles - _MCM_DEMO_T_REF_K

    fiber_mcm     = _mcm_demo_fiber(ΔT_K_particles)
    modified_path = fiber_mcm.path  # thermal already applied at construction
    fiber_mod     = Fiber(modified_path; cross_section = _MCM_DEMO_XS, T_ref_K = T_K_particles)
    J_p, _        = propagate_fiber(
        fiber_mod;
        λ_m = _MCM_DEMO_λ_M,
    )
    MonteCarloMeasurements.unsafe_comparisons(false)

    N = 50
    T_C_samples = T_C_particles.particles
    s1_p = zeros(N); s2_p = zeros(N); s3_p = zeros(N)
    dlp_p = zeros(N); angle_p = zeros(N)

    for k in 1:N
        J_k = ComplexF64[
            real(J_p[1,1]).particles[k] + im*imag(J_p[1,1]).particles[k]   real(J_p[1,2]).particles[k] + im*imag(J_p[1,2]).particles[k];
            real(J_p[2,1]).particles[k] + im*imag(J_p[2,1]).particles[k]   real(J_p[2,2]).particles[k] + im*imag(J_p[2,2]).particles[k]
        ]
        st = _jones_to_stokes(J_k)
        s1_p[k]    = st.s1
        s2_p[k]    = st.s2
        s3_p[k]    = st.s3
        dlp_p[k]   = st.dlp
        angle_p[k] = mod(0.5 * atan(st.s2, st.s1) * 180.0 / π, 180.0)
    end

    js_arr(v) = "[" * join(string.(v), ",") * "]"

    html = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>MCM temperature PTF — StaticParticles(50)</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
html,body{margin:0;padding:20px;background:#111;color:#eee;font-family:sans-serif;}
h2{color:#aaa;font-size:1.1em;margin:0 0 8px 0;}
.row{display:flex;gap:10px;flex-wrap:wrap;}
.cell{flex:1 1 420px;}
</style>
</head>
<body>
<h2>MCM demo — temperature-dependent PTF (StaticParticles(50), T = $(_MCM_DEMO_T_NOM_C)°C ± $(_MCM_DEMO_T_SIG_C)°C, λ = 1550 nm)</h2>
<p style="color:#888;font-size:0.9em;margin:0 0 12px 0;">
  Each point is one Monte-Carlo particle. Colour encodes temperature. N = $(N) particles.
  Fiber: straight 5 m → helix (D=0.05 m, ~10002 turns, p=0.05 m, T-sensitive) → straight 5 m →
  helix (D=0.05 m, 10000 turns, p=0.05 m) → straight 5 m.
  T_nom = 30°C, mid-fringe design point. Input state: horizontal (H) polarisation.
  The output remains nearly linearly polarised (DLP ≈ 1) because the helix birefringence
  rotates the linear polarisation angle without introducing significant ellipticity (S3 ≈ 0).
  The sensitive observable is therefore the <b>polarisation angle</b>, not DLP.
</p>
<div class="row">
  <div class="cell" id="angle_plot"></div>
  <div class="cell" id="stokes_plot"></div>
</div>
<script>
const T = $(js_arr(T_C_samples));
const s1 = $(js_arr(s1_p));
const s2 = $(js_arr(s2_p));
const s3 = $(js_arr(s3_p));
const dlp = $(js_arr(dlp_p));
const angle = $(js_arr(angle_p));

const layout_dark = {
  paper_bgcolor: '#1a1a1a',
  plot_bgcolor:  '#1a1a1a',
  font: {color: '#ccc'},
  xaxis: {gridcolor:'#333', color:'#aaa', title:'Temperature (°C)'},
  margin: {l:60, r:20, t:50, b:50},
  legend: {font:{color:'#ccc'}}
};

const marker_opts = {size: 6, color: T, colorscale: 'Viridis',
                     colorbar: {title: 'T (°C)', tickfont:{color:'#ccc'}, titlefont:{color:'#ccc'}},
                     opacity: 0.85};

Plotly.newPlot('angle_plot',
  [{x: T, y: angle, mode: 'markers', marker: marker_opts,
    name: 'pol. angle', hovertemplate: 'T=%{x:.2f}°C<br>angle=%{y:.2f}°<extra></extra>'}],
  Object.assign({}, layout_dark, {
    title: {text: 'Linear Polarisation Angle vs Temperature', font:{color:'#ddd'}},
    yaxis: {gridcolor:'#333', color:'#aaa', title:'Angle (deg)', range:[0,180]}
  })
);

Plotly.newPlot('stokes_plot',
  [
    {x: T, y: s1, mode: 'markers', name: 'S1',
     marker: Object.assign({}, marker_opts, {colorscale:'Reds'}),
     hovertemplate: 'T=%{x:.2f}°C<br>S1=%{y:.4f}<extra></extra>'},
    {x: T, y: s2, mode: 'markers', name: 'S2',
     marker: Object.assign({}, marker_opts, {colorscale:'Greens'}),
     hovertemplate: 'T=%{x:.2f}°C<br>S2=%{y:.4f}<extra></extra>'},
    {x: T, y: s3, mode: 'markers', name: 'S3',
     marker: Object.assign({}, marker_opts, {colorscale:'Blues'}),
     hovertemplate: 'T=%{x:.2f}°C<br>S3=%{y:.4f}<extra></extra>'},
    {x: T, y: dlp, mode: 'markers', name: 'DLP',
     marker: Object.assign({}, marker_opts, {colorscale:'Oranges', symbol:'diamond'}),
     hovertemplate: 'T=%{x:.2f}°C<br>DLP=%{y:.4f}<extra></extra>'}
  ],
  Object.assign({}, layout_dark, {
    title: {text: 'Stokes Parameters vs Temperature', font:{color:'#ddd'}},
    yaxis: {gridcolor:'#333', color:'#aaa', title:'Stokes parameter / DLP', range:[-1,1]}
  })
);
</script>
</body></html>
"""

    out = joinpath(output_dir, "mcm-temperature-ptf.html")
    open(out, "w") do io
        write(io, html)
    end
    println("Wrote MCM temperature PTF demo to: ", out)
    return (path = out, desc = desc)
end

"""
    demo_mcm_temperature_ptf_scatter(; output_dir = …)

MCM Particles run: build one fiber with MCMadd(:T_K, ΔT) where
ΔT ~ N(0, 5 K) (30 °C ± 5 °C), propagate once, and show the
per-particle DLP and Stokes scatter on a Poincaré equatorial projection.
"""
function demo_mcm_temperature_ptf_scatter(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    desc = "MCM temperature demo: per-particle DLP scatter from a single " *
           "Particles run with T = 30°C ± 5°C (1σ) on the first large helix."

    MonteCarloMeasurements.unsafe_comparisons(true)
    T_C_particles  = _MCM_DEMO_T_NOM_C ± _MCM_DEMO_T_SIG_C  # Particles{Float64,2000}
    T_K_particles  = T_C_particles + 273.15
    ΔT_K_particles = T_K_particles - _MCM_DEMO_T_REF_K

    fiber_mcm = _mcm_demo_fiber(ΔT_K_particles)
    modified_path = fiber_mcm.path  # thermal already applied at construction
    # T_ref_K = Particles so material birefringence is also uncertain.
    fiber_mod = Fiber(modified_path;
                      cross_section = _MCM_DEMO_XS,
                      T_ref_K       = T_K_particles)

    J_p, _ = propagate_fiber(
        fiber_mod;
        λ_m = _MCM_DEMO_λ_M,
    )
    MonteCarloMeasurements.unsafe_comparisons(false)

    # Extract per-particle output state: apply J to [1,0] for each particle.
    N = nparticles(real(J_p[1, 1]))

    T_C_samples = T_C_particles.particles
    s1_p = zeros(N); s2_p = zeros(N); s3_p = zeros(N)
    dlp_p = zeros(N); angle_p = zeros(N)

    for k in 1:N
        J_k = ComplexF64[
            real(J_p[1,1]).particles[k] + im*imag(J_p[1,1]).particles[k]   real(J_p[1,2]).particles[k] + im*imag(J_p[1,2]).particles[k];
            real(J_p[2,1]).particles[k] + im*imag(J_p[2,1]).particles[k]   real(J_p[2,2]).particles[k] + im*imag(J_p[2,2]).particles[k]
        ]
        st = _jones_to_stokes(J_k)
        s1_p[k]    = st.s1
        s2_p[k]    = st.s2
        s3_p[k]    = st.s3
        dlp_p[k]   = st.dlp
        angle_p[k] = mod(0.5 * atan(st.s2, st.s1) * 180.0 / π, 180.0)
    end

    js_arr(v)  = "[" * join(string.(v), ",") * "]"

    html = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>MCM temperature PTF — Particles scatter</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
html,body{margin:0;padding:20px;background:#111;color:#eee;font-family:sans-serif;}
h2{color:#aaa;font-size:1.1em;margin:0 0 8px 0;}
.row{display:flex;gap:10px;flex-wrap:wrap;}
.cell{flex:1 1 400px;}
</style>
</head>
<body>
<h2>MCM demo — temperature PTF scatter (T = $(_MCM_DEMO_T_NOM_C)°C ± $(_MCM_DEMO_T_SIG_C)°C, λ = 1550 nm)</h2>
<p style="color:#888;font-size:0.9em;margin:0 0 12px 0;">
  Each point is one Monte-Carlo particle.  Colour encodes temperature.  N = $(N) particles.
  The output is nearly linearly polarised throughout (DLP ≈ 1); the sensitive observable
  is the <b>linear polarisation angle</b>, shown left.  The right panel shows the
  Poincaré equatorial projection (S1–S2): particles trace an arc as temperature varies.
</p>
<div class="row">
  <div class="cell" id="angle_scatter"></div>
  <div class="cell" id="poincare_scatter"></div>
</div>
<script>
const T_samp = $(js_arr(T_C_samples));
const s1_p   = $(js_arr(s1_p));
const s2_p   = $(js_arr(s2_p));
const angle_p = $(js_arr(angle_p));

const layout_dark = {
  paper_bgcolor: '#1a1a1a',
  plot_bgcolor:  '#1a1a1a',
  font: {color: '#ccc'},
  margin: {l:60, r:20, t:50, b:50}
};

Plotly.newPlot('angle_scatter',
  [{
    x: T_samp, y: angle_p,
    mode: 'markers',
    marker: {size: 4, color: T_samp, colorscale: 'Viridis',
             colorbar: {title: 'T (°C)', tickfont:{color:'#ccc'}, titlefont:{color:'#ccc'}},
             opacity: 0.7},
    hovertemplate: 'T=%{x:.2f}°C<br>angle=%{y:.2f}°<extra></extra>'
  }],
  Object.assign({}, layout_dark, {
    title: {text: 'Linear Polarisation Angle vs Temperature (MCM scatter)', font:{color:'#ddd'}},
    xaxis: {gridcolor:'#333', color:'#aaa', title:'Temperature (°C)'},
    yaxis: {gridcolor:'#333', color:'#aaa', title:'Angle (deg)', range:[0,180]}
  })
);

Plotly.newPlot('poincare_scatter',
  [{
    x: s1_p, y: s2_p,
    mode: 'markers',
    marker: {size: 4, color: T_samp, colorscale: 'Viridis', opacity: 0.7,
             colorbar: {title: 'T (°C)', tickfont:{color:'#ccc'}, titlefont:{color:'#ccc'}}},
    hovertemplate: 'S1=%{x:.3f}<br>S2=%{y:.3f}<br>T=%{customdata:.2f}°C<extra></extra>',
    customdata: T_samp
  }],
  Object.assign({}, layout_dark, {
    title: {text: 'Poincaré equatorial projection (S1–S2) — MCM scatter', font:{color:'#ddd'}},
    xaxis: {gridcolor:'#333', color:'#aaa', title:'S1', range:[-1,1]},
    yaxis: {gridcolor:'#333', color:'#aaa', title:'S2', range:[-1,1],
            scaleanchor:'x', scaleratio:1}
  })
);
</script>
</body></html>
"""

    out = joinpath(output_dir, "mcm-temperature-ptf-scatter.html")
    open(out, "w") do io
        write(io, html)
    end
    println("Wrote MCM temperature PTF scatter demo to: ", out)
    return (path = out, desc = desc)
end

# =====================================================================
# Monolithic index entries
# =====================================================================

const DEMO3MCM_INDEX = [
    (group = "MCM", fn = demo_mcm_temperature_ptf,         kwargs = (;)),
    (group = "MCM", fn = demo_mcm_temperature_ptf_scatter, kwargs = (;)),
]

"""
    demo3mcm_all(; index_output)

Run every demo in `DEMO3MCM_INDEX` and write `demo-index.html`.
"""
function demo3mcm_entries()
    entries = Tuple{String, String, String, String}[]

    for d in DEMO3MCM_INDEX
        println("[ demo3mcm ] $(d.fn)")
        result = d.fn(; d.kwargs...)
        desc = _demo_result_desc(result, d)
        for path in _demo_html_paths(result)
            push!(entries, (d.group, basename(path), path, desc))
        end
    end
    return entries
end

function demo3mcm_all(; index_output::AbstractString = DEMO_MONOLITHIC_INDEX_OUTPUT)
    return _write_demo_index(
        [(title = "MCM temperature PTF demos",
          entries = demo3mcm_entries(),
          group_titles = Dict("MCM" => "MCM / temperature-dependent PTF"))];
        index_output,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo3mcm_all()
end
