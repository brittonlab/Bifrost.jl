# =====================================================================
# demo-fiber-padles.jl — fiber polarization-paddle controller demo
# =====================================================================
#
# Builds a three-paddle Lefèvre ("bat-ear") fiber polarization controller
# configured as λ/4, λ/2, λ/4 at 1550 nm and shows that sweeping the two
# inter-paddle orientations steers the output polarization across a wide
# region of the Poincaré sphere — the hallmark of a "universal" controller.
#
#   * demo_fiber_paddles — one propagation per (τ1, τ2) grid point; writes
#                          output/fiber-paddles.html (3D Poincaré scatter +
#                          Stokes-vs-sweep traces) and prints the achieved
#                          per-paddle retardance.
#
# `demo_fiber_paddles_all()` runs it and writes output/demo-fiber-padles.html.
#
# Physics
# ───────
# A paddle is a planar coil of fiber. Its bend-induced (photoelastic)
# birefringence Δβ ∝ 1/R² makes the coil act as a linear retarder (waveplate).
# For one full planar loop the retardance is
#
#     δ = |Δβ(R)| · (N · 2πR),         Δβ(R) = -A / R²,   A ≡ |Δβ|·R²,
#
# so δ = 2πN·A / R. With the SMF-28-like cross section below, A ≈ 2.1117e-3,
# and choosing R = 4A makes a single loop an exact λ/4 plate (δ = π/2) and two
# loops an exact λ/2 plate (δ = π). A 1-loop / 2-loop / 1-loop stack is then
# the λ/4, λ/2, λ/4 controller built in real 3-paddle hardware.
#
# A planar full loop is built with `bend!(radius = R, angle = N·2π)`, NOT
# `helix!`: a helix's curvature vector rotates through the loop and the
# waveplate averages out. The relative orientation of one paddle's fast axis
# to the next is set by a `Spinning` run on the connecting straight: spinning
# rotates the local frame, so the next bend's curvature azimuth (its c2φ/s2φ in
# the bend generator) is rotated — exactly how rotating a physical paddle
# rotates its waveplate axis.

using Bifrost
using LinearAlgebra

# demo_fiber_paddles_all uses the shared monolithic-index writer.
if !isdefined(Main, :_write_demo_index)
    include(joinpath(@__DIR__, "demo-index-helpers.jl"))
end

# =====================================================================
# Design constants (verified against the live model)
# =====================================================================

# Standard SMF-28-like step-index cross section.
const _PADDLE_XS = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-28-like",
)
const _PADDLE_λ_M = 1550e-9
const _PADDLE_T_K = 297.15

# Solved design radius: with this cross section A ≡ |Δβ|·R² ≈ 2.1117e-3, so
# R = 4A makes one planar loop an exact λ/4 plate at 1550 nm (verified:
# J = diag(e^{-iπ/4}, e^{+iπ/4})).
const _PADDLE_R_M = 8.446716876988274e-3

# Loop counts per paddle → λ/4, λ/2, λ/4.
const _PADDLE_LOOPS = (1, 2, 1)

# Length of each inter-paddle straight that carries the orientation `Spinning`.
const _PADDLE_TWIST_LEN_M = 0.05

# Straight lead-in / lead-out length (m).
const _PADDLE_LEAD_M = 0.02

# Inter-paddle twist-rate sweep (rad/m). A frame rotation of π covers every
# distinct waveplate azimuth (the bend generator depends on 2φ), reached at
# τ = π / _PADDLE_TWIST_LEN_M ≈ 63 rad/m; the upper bound here spans ≳ 2 full
# basis rotations so the swept states cover a wide region of the sphere.
const _PADDLE_TWIST_MAX = 130.0
const _PADDLE_SWEEP_N = 11

# =====================================================================
# Builder helper (local to this demo, not a shared library helper)
# =====================================================================

"""
    _paddle_controller_fiber(twist1, twist2)

Build the three-paddle λ/4, λ/2, λ/4 controller as a single `Fiber`:

    straight (lead-in)
    bend! angle = 1·2π          (paddle 1: λ/4)
    straight + Spinning(τ=twist1)
    bend! angle = 2·2π          (paddle 2: λ/2)
    straight + Spinning(τ=twist2)
    bend! angle = 1·2π          (paddle 3: λ/4)
    straight (lead-out)

`twist1` / `twist2` are constant spinning rates (rad/m) applied on the two
connecting straights; they set the relative fast-axis orientation of the
following paddle.
"""
function _paddle_controller_fiber(twist1::Real, twist2::Real)
    sb = SubpathBuilder()
    start!(sb)
    straight!(sb; length = _PADDLE_LEAD_M, meta = [Nickname("lead-in")])
    bend!(sb; radius = _PADDLE_R_M, angle = _PADDLE_LOOPS[1] * 2π,
          meta = [Nickname("paddle 1: λ/4")])
    straight!(sb; length = _PADDLE_TWIST_LEN_M,
              meta = AbstractMeta[Nickname("twist 1"),
                                  Spinning(; rate = s -> Float64(twist1))])
    bend!(sb; radius = _PADDLE_R_M, angle = _PADDLE_LOOPS[2] * 2π,
          meta = [Nickname("paddle 2: λ/2")])
    straight!(sb; length = _PADDLE_TWIST_LEN_M,
              meta = AbstractMeta[Nickname("twist 2"),
                                  Spinning(; rate = s -> Float64(twist2))])
    bend!(sb; radius = _PADDLE_R_M, angle = _PADDLE_LOOPS[3] * 2π,
          meta = [Nickname("paddle 3: λ/4")])
    # A final Spinning(0) closes the second run cleanly at the lead-out.
    straight!(sb; length = _PADDLE_LEAD_M,
              meta = AbstractMeta[Nickname("lead-out"),
                                  Spinning(; rate = s -> 0.0)])
    seal!(sb)
    return Fiber(build(sb); cross_section = _PADDLE_XS, T_ref_K = _PADDLE_T_K)
end

# Per-paddle retardance δ = |Δβ(R)| · (N · 2πR), computed from first principles.
function _paddle_retardance(loops::Integer)
    Δβ = bending_birefringence(_PADDLE_XS, _PADDLE_λ_M, _PADDLE_T_K;
                               bend_radius_m = _PADDLE_R_M)
    return abs(Δβ) * (loops * 2π * _PADDLE_R_M)
end

# Apply a 2×2 Jones matrix to the input state [1, 0] (horizontal) and return
# normalised Stokes parameters. Scalar ComplexF64 entries only.
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

# =====================================================================
# Demo — one function, one output file
# =====================================================================

"""
    demo_fiber_paddles(; output_dir = …)

Sweep the two inter-paddle twist rates over a grid, propagate the λ/4–λ/2–λ/4
controller once per grid point, and plot the swept output polarization states
on a 3D Poincaré sphere (plus S1/S2/S3 vs sweep index). The page and console
also report each paddle's achieved retardance, confirming the λ/4, λ/2, λ/4
targets from first principles.
"""
function demo_fiber_paddles(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    desc = "Fiber polarization paddles: a λ/4–λ/2–λ/4 Lefèvre controller " *
           "(bend-induced waveplates, R = $(round(_PADDLE_R_M * 1e3, digits = 3)) mm, " *
           "1/2/1 planar loops) swept over its two inter-paddle orientations, " *
           "showing the output state covering a wide region of the Poincaré sphere."

    # Retardances from first principles (independent of the propagation).
    δ1 = _paddle_retardance(_PADDLE_LOOPS[1])
    δ2 = _paddle_retardance(_PADDLE_LOOPS[2])
    δ3 = _paddle_retardance(_PADDLE_LOOPS[3])
    println("Paddle retardances at λ = 1550 nm, R = ",
            round(_PADDLE_R_M * 1e3, digits = 3), " mm:")
    println("  paddle 1 (", _PADDLE_LOOPS[1], " loop):  δ = ", round(δ1, digits = 6),
            " rad = ", round(δ1 / (π / 2), digits = 4), "·(π/2)   [target λ/4 = π/2]")
    println("  paddle 2 (", _PADDLE_LOOPS[2], " loops): δ = ", round(δ2, digits = 6),
            " rad = ", round(δ2 / π, digits = 4), "·π        [target λ/2 = π]")
    println("  paddle 3 (", _PADDLE_LOOPS[3], " loop):  δ = ", round(δ3, digits = 6),
            " rad = ", round(δ3 / (π / 2), digits = 4), "·(π/2)   [target λ/4 = π/2]")

    # Sweep the two inter-paddle twist rates.
    τs = collect(range(0.0, _PADDLE_TWIST_MAX, length = _PADDLE_SWEEP_N))
    npts = _PADDLE_SWEEP_N^2
    s1 = zeros(npts); s2 = zeros(npts); s3 = zeros(npts)
    idx = zeros(Int, npts)
    offdiag = zeros(npts)

    k = 0
    for (i, τ1) in enumerate(τs), (j, τ2) in enumerate(τs)
        k += 1
        fiber = _paddle_controller_fiber(τ1, τ2)
        J, _ = propagate_fiber(fiber; λ_m = _PADDLE_λ_M, verbose = false)
        st = _jones_to_stokes(J)
        s1[k] = st.s1; s2[k] = st.s2; s3[k] = st.s3
        idx[k] = k
        offdiag[k] = abs(J[1, 2])
    end

    # Sanity report: aligned (τ=0) stack is a pure retarder (J off-diagonal ≈ 0).
    println("Off-diagonal |J₁₂| at (τ1,τ2)=(0,0): ", round(offdiag[1], digits = 4),
            "  (pure retarder ⇒ ≈ 0)")
    println("Max off-diagonal |J₁₂| over sweep:    ", round(maximum(offdiag), digits = 4))

    js_arr(v) = "[" * join(string.(v), ",") * "]"

    # Translucent reference sphere (low-res mesh) for the Poincaré display.
    nu = 24; nv = 24
    us = range(0, 2π, length = nu)
    vs = range(0, π, length = nv)
    sx = [cos(u) * sin(v) for v in vs, u in us]
    sy = [sin(u) * sin(v) for v in vs, u in us]
    sz = [cos(v) for v in vs, u in us]
    js_mat(M) = "[" * join(["[" * join(string.(M[r, :]), ",") * "]"
                            for r in 1:size(M, 1)], ",") * "]"

    δ1s = round(δ1 / (π / 2), digits = 3)
    δ2s = round(δ2 / π, digits = 3)
    δ3s = round(δ3 / (π / 2), digits = 3)

    html = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>Fiber polarization paddles — λ/4, λ/2, λ/4 controller</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
html,body{margin:0;padding:20px;background:#111;color:#eee;font-family:sans-serif;}
h2{color:#aaa;font-size:1.1em;margin:0 0 8px 0;}
.row{display:flex;gap:10px;flex-wrap:wrap;}
.cell{flex:1 1 460px;}
.targets{color:#9ad;font-size:0.95em;margin:0 0 12px 0;}
code{color:#bbb;}
</style>
</head>
<body>
<h2>Fiber polarization paddles — λ/4, λ/2, λ/4 Lefèvre controller (λ = 1550 nm)</h2>
<p style="color:#888;font-size:0.9em;margin:0 0 6px 0;">
  Three planar coils of SMF-28-like fiber at R = $(round(_PADDLE_R_M*1e3, digits=3)) mm act as
  bend-induced waveplates: $(_PADDLE_LOOPS[1]) loop, $(_PADDLE_LOOPS[2]) loops, $(_PADDLE_LOOPS[3])
  loop ⇒ λ/4, λ/2, λ/4. Each point is one (τ1, τ2) inter-paddle orientation; the two connecting
  straights carry a constant <code>Spinning</code> rate swept over
  [0, $(_PADDLE_TWIST_MAX)] rad/m ($(_PADDLE_SWEEP_N)×$(_PADDLE_SWEEP_N) grid). Input state:
  horizontal (H). The swept output states cover a wide region of the Poincaré sphere — a
  universal polarization controller.
</p>
<p class="targets">
  Achieved retardance (from |Δβ(R)|·N·2πR):
  paddle 1 = $(δ1s)·(π/2),&nbsp; paddle 2 = $(δ2s)·π,&nbsp; paddle 3 = $(δ3s)·(π/2)
  &nbsp;⇒ λ/4, λ/2, λ/4.
</p>
<div class="row">
  <div class="cell" id="poincare"></div>
  <div class="cell" id="stokes"></div>
</div>
<script>
const s1 = $(js_arr(s1));
const s2 = $(js_arr(s2));
const s3 = $(js_arr(s3));
const idx = $(js_arr(idx));
const sx = $(js_mat(sx));
const sy = $(js_mat(sy));
const sz = $(js_mat(sz));

const sphere = {
  type: 'surface', x: sx, y: sy, z: sz,
  opacity: 0.18, showscale: false,
  colorscale: [[0,'#335'],[1,'#335']], hoverinfo: 'skip'
};

const pts = {
  type: 'scatter3d', mode: 'markers',
  x: s1, y: s2, z: s3,
  marker: {size: 4, color: idx, colorscale: 'Viridis',
           colorbar: {title: 'sweep index', tickfont:{color:'#ccc'},
                      titlefont:{color:'#ccc'}}, opacity: 0.9},
  hovertemplate: 'S1=%{x:.3f}<br>S2=%{y:.3f}<br>S3=%{z:.3f}<extra></extra>'
};

Plotly.newPlot('poincare', [sphere, pts], {
  paper_bgcolor: '#1a1a1a',
  font: {color: '#ccc'},
  title: {text: 'Output states on the Poincaré sphere', font:{color:'#ddd'}},
  margin: {l:0, r:0, t:40, b:0},
  scene: {
    xaxis: {title:'S1', range:[-1,1], color:'#aaa', gridcolor:'#333',
            backgroundcolor:'#1a1a1a', showbackground:true},
    yaxis: {title:'S2', range:[-1,1], color:'#aaa', gridcolor:'#333',
            backgroundcolor:'#1a1a1a', showbackground:true},
    zaxis: {title:'S3', range:[-1,1], color:'#aaa', gridcolor:'#333',
            backgroundcolor:'#1a1a1a', showbackground:true},
    aspectmode: 'cube'
  }
});

const layout_dark = {
  paper_bgcolor: '#1a1a1a',
  plot_bgcolor:  '#1a1a1a',
  font: {color: '#ccc'},
  margin: {l:60, r:20, t:50, b:50},
  legend: {font:{color:'#ccc'}},
  title: {text: 'Stokes parameters vs sweep index', font:{color:'#ddd'}},
  xaxis: {gridcolor:'#333', color:'#aaa', title:'sweep index (τ1 outer, τ2 inner)'},
  yaxis: {gridcolor:'#333', color:'#aaa', title:'Stokes parameter', range:[-1,1]}
};

Plotly.newPlot('stokes', [
  {x: idx, y: s1, mode:'markers', name:'S1',
   marker:{size:4, color:'#e66'}, hovertemplate:'idx=%{x}<br>S1=%{y:.3f}<extra></extra>'},
  {x: idx, y: s2, mode:'markers', name:'S2',
   marker:{size:4, color:'#6e6'}, hovertemplate:'idx=%{x}<br>S2=%{y:.3f}<extra></extra>'},
  {x: idx, y: s3, mode:'markers', name:'S3',
   marker:{size:4, color:'#69e'}, hovertemplate:'idx=%{x}<br>S3=%{y:.3f}<extra></extra>'}
], layout_dark);
</script>
</body></html>
"""

    mkpath(output_dir)
    out = joinpath(output_dir, "fiber-paddles.html")
    open(out, "w") do io
        write(io, html)
    end
    println("Wrote fiber-paddles demo to: ", out)
    return (path = out, desc = desc)
end

# =====================================================================
# Monolithic index entries
# =====================================================================

const DEMO_FIBER_PADDLES_INDEX = [
    (group = "paddles", fn = demo_fiber_paddles, kwargs = (;)),
]

"""
    demo_fiber_paddles_entries()

Run every demo in `DEMO_FIBER_PADDLES_INDEX` and return `(group, link_title,
path, desc)` tuples for the monolithic demo index.
"""
function demo_fiber_paddles_entries()
    entries = Tuple{String, String, String, String}[]
    for d in DEMO_FIBER_PADDLES_INDEX
        println("[ demo-fiber-padles ] $(d.fn)")
        result = d.fn(; d.kwargs...)
        desc = _demo_result_desc(result, d)
        for path in _demo_html_paths(result)
            push!(entries, (d.group, basename(path), path, desc))
        end
    end
    return entries
end

"""
    demo_fiber_paddles_all(; index_output)

Run the paddle demo and write a small index page.
"""
function demo_fiber_paddles_all(;
    index_output::AbstractString = joinpath(@__DIR__, "..", "..", "output",
                                            "demo-fiber-padles.html"),
)
    return _write_demo_index(
        [(title = "Fiber polarization-paddle demos",
          source_file = "demo-fiber-padles.jl",
          entries = demo_fiber_paddles_entries(),
          group_titles = Dict("paddles" => "Fiber polarization paddles"))];
        index_output,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo_fiber_paddles_all()
end
