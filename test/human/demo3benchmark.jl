# =====================================================================
# demo3benchmark.jl — MCM propagation speed benchmarks
# =====================================================================
#
# Extends the MCM temperature-PTF scenarios from demo3mcm.jl with
# timing measurements across four particle types:
#
#   * Float64          — deterministic baseline
#   * Particles(2000)  — heap-allocated, arbitrary N, standard MCM
#   * StaticParticles(50/100/200) — SVector-backed, SIMD-friendly
#
# For each variant two times are measured:
#   * first-call  : wall time of the very first call, i.e. JIT
#                   compilation overhead + one execution.
#   * steady-state: minimum wall time over many repetitions after JIT
#                   (using BenchmarkTools @belapsed).
#
# M1 SIMD note: Apple M1 NEON has 128-bit lanes (2×f64/lane).
# StaticParticles shines when N is small enough for Julia to emit
# SIMD-vectorised loops over the SVector.  On M1, the sweet spot is
# roughly N = 50–200; above ~300 the stack/register pressure grows
# and Particles(N) can match or beat StaticParticles(N).
#
# Two HTML files are written:
#   benchmark-mcm-propagate.html  — bar chart + table for propagate_fiber
#   benchmark-mcm-modify.html     — bar chart + table for modify + propagate_fiber
#
# `demo3benchmark_all()` runs both and writes `output/demo-index.html`.

if !isdefined(Main, :_mcm_demo_fiber)
    include(joinpath(@__DIR__, "demo3mcm.jl"))
end

using BenchmarkTools
using Distributions

# =====================================================================
# Helpers
# =====================================================================

# Construct an MCM-annotated fiber ready for propagate_fiber.
# `make_T` is a zero-argument function that returns a temperature value
# (Float64, Particles, or StaticParticles) at 30°C.
function _bench_build_fiber(make_T)
    T_val  = make_T()
    T_K    = T_val + 273.15
    ΔT_K   = T_K - _MCM_DEMO_T_REF_K
    fiber  = _mcm_demo_fiber(ΔT_K)            # thermal applied at construction
    return Fiber(fiber.path; cross_section = _MCM_DEMO_XS, T_ref_K = T_K)
end

# Time the first call to `f()` (JIT + one run) in milliseconds.
function _first_call_ms(f)
    t0 = time_ns()
    f()
    return (time_ns() - t0) / 1e6
end

# Steady-state minimum wall time in milliseconds using BenchmarkTools.
function _steady_ms(f)
    return @belapsed($f()) * 1e3
end

# =====================================================================
# Benchmark cases
# =====================================================================

# Each case: (label, N_particles_or_0, make_T_fn)
# make_T returns a fresh value each time it is called.
function _bench_cases()
    return [
        ("Float64",             0,    () -> Float64(_MCM_DEMO_T_NOM_C)),
        ("Particles(2000)",     2000, () -> Particles(2000, Normal(Float64(_MCM_DEMO_T_NOM_C), Float64(_MCM_DEMO_T_SIG_C)))),
        ("StaticParticles(50)", 50,   () -> StaticParticles(50,  Normal(Float64(_MCM_DEMO_T_NOM_C), Float64(_MCM_DEMO_T_SIG_C)))),
        ("StaticParticles(75)", 75,  () -> StaticParticles(75, Normal(Float64(_MCM_DEMO_T_NOM_C), Float64(_MCM_DEMO_T_SIG_C)))),
    ]
end

# Run all benchmarks.  Returns a Vector of NamedTuples with fields:
#   label, n_particles, first_ms, steady_ms
function _run_benchmarks(; scenario::Symbol = :propagate)
    cases  = _bench_cases()
    results = NamedTuple[]

    MonteCarloMeasurements.unsafe_comparisons(true)
    for (label, n_p, make_T) in cases
        println("  Benchmarking: $label ...")

        if scenario === :propagate
            # Benchmark: propagate_fiber only (fiber pre-built)
            fiber  = _bench_build_fiber(make_T)
            fn     = () -> propagate_fiber(fiber; λ_m = _MCM_DEMO_λ_M,
                              params = SolverParams(rtol = 1e-5, atol = 1e-9, h_min = 1e-12))
        else
            # Benchmark: fiber build (incl. :T_K thermal scaling) + propagate
            T_val  = make_T()
            T_K    = T_val + 273.15
            ΔT_K   = T_K - _MCM_DEMO_T_REF_K
            fn     = () -> begin
                # The build (with thermal scaling) happens here so the benchmark
                # covers fiber construction and propagation together.
                f2 = Fiber(_mcm_demo_fiber(ΔT_K).path;
                           cross_section = _MCM_DEMO_XS, T_ref_K = T_K)
                propagate_fiber(f2; λ_m = _MCM_DEMO_λ_M,
                                params = SolverParams(rtol = 1e-5, atol = 1e-9, h_min = 1e-12))
            end
        end

        fc  = _first_call_ms(fn)   # includes JIT
        fn()                        # second call warms any remaining caches
        ss  = _steady_ms(fn)

        push!(results, (; label, n_particles = n_p, first_ms = fc, steady_ms = ss))
    end
    MonteCarloMeasurements.unsafe_comparisons(false)

    return results
end

# =====================================================================
# HTML renderer
# =====================================================================

function _bench_html(results, title, subtitle, output_path)
    labels     = [r.label      for r in results]
    first_ms   = [r.first_ms   for r in results]
    steady_ms  = [r.steady_ms  for r in results]
    n_p        = [r.n_particles for r in results]

    js_strarr(v) = "[" * join(("\"$(x)\"" for x in v), ",") * "]"
    js_numarr(v) = "[" * join((string(round(x; digits=3)) for x in v), ",") * "]"

    # Build table rows
    table_rows = join([
        """<tr>
          <td>$(labels[i])</td>
          <td>$(n_p[i] == 0 ? "—" : string(n_p[i]))</td>
          <td>$(round(first_ms[i]; digits=1))</td>
          <td>$(round(steady_ms[i]; digits=3))</td>
          <td>$(round(first_ms[i]/first_ms[1]; digits=1))×</td>
          <td>$(round(steady_ms[i]/steady_ms[1]; digits=1))×</td>
        </tr>"""
        for i in eachindex(results)
    ], "\n")

    html = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>$(title)</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
html,body{margin:0;padding:20px;background:#111;color:#eee;font-family:sans-serif;}
h2{color:#aaa;font-size:1.1em;margin:0 0 6px 0;}
p.sub{color:#777;font-size:0.88em;margin:0 0 14px 0;}
#chart{height:420px;margin-bottom:24px;}
table{border-collapse:collapse;font-size:0.9em;width:100%;}
th,td{padding:6px 12px;border:1px solid #333;text-align:right;}
th{background:#222;color:#ccc;text-align:center;}
td:first-child{text-align:left;}
tr:hover td{background:#1a2a1a;}
</style>
</head>
<body>
<h2>$(title)</h2>
<p class="sub">$(subtitle)</p>
<div id="chart"></div>
<table>
  <thead><tr>
    <th>Variant</th><th>N particles</th>
    <th>First-call (ms)<br><small>JIT + run</small></th>
    <th>Steady-state (ms)<br><small>post-JIT min</small></th>
    <th>First-call ratio</th>
    <th>Steady-state ratio</th>
  </tr></thead>
  <tbody>$(table_rows)</tbody>
</table>
<script>
const labels = $(js_strarr(labels));
const first  = $(js_numarr(first_ms));
const steady = $(js_numarr(steady_ms));

const colors_first  = '#5af';
const colors_steady = '#f80';

const layout = {
  paper_bgcolor: '#1a1a1a',
  plot_bgcolor:  '#1a1a1a',
  font:   {color: '#ccc'},
  title:  {text: '$(title)', font: {color: '#ddd'}},
  barmode: 'group',
  xaxis:  {gridcolor: '#333', color: '#aaa'},
  yaxis:  {gridcolor: '#333', color: '#aaa', title: 'Wall time (ms)', type: 'log',
           tickformat: '.3g'},
  margin: {l:70, r:20, t:50, b:120},
  legend: {font: {color: '#ccc'}},
  bargap: 0.15, bargroupgap: 0.08
};

Plotly.newPlot('chart', [
  {
    name: 'First-call (JIT + run)',
    x: labels, y: first,
    type: 'bar',
    marker: {color: '#5af', opacity: 0.85},
    hovertemplate: '%{x}<br>First-call: %{y:.1f} ms<extra></extra>'
  },
  {
    name: 'Steady-state (post-JIT)',
    x: labels, y: steady,
    type: 'bar',
    marker: {color: '#f80', opacity: 0.85},
    hovertemplate: '%{x}<br>Steady-state: %{y:.3f} ms<extra></extra>'
  }
], layout);
</script>
</body></html>
"""
    open(output_path, "w") do io
        write(io, html)
    end
    return output_path
end

# =====================================================================
# Demo functions (one per HTML file)
# =====================================================================

"""
    demo_benchmark_mcm_propagate(; output_dir = …)

Benchmark `propagate_fiber` alone across Float64, Particles(2000),
and StaticParticles(50/100/200).  The fiber is pre-built; only the
ODE propagation is timed.
"""
function demo_benchmark_mcm_propagate(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    desc = "MCM benchmark: propagate_fiber wall time across Float64, " *
           "Particles(2000), StaticParticles(50/100/200).  " *
           "The reference fiber includes spinning sensitivity in addition to " *
           "temperature dependence.  First-call includes JIT; steady-state is " *
           "post-JIT minimum."

    println("Running propagate_fiber benchmarks …")
    results = _run_benchmarks(; scenario = :propagate)

    out = joinpath(output_dir, "benchmark-mcm-propagate.html")
    _bench_html(
        results,
        "MCM benchmark — propagate_fiber (log scale)",
        "Fiber: straight 5 m → helix (D=0.05 m, ~10002 turns, spinning + " *
        "temperature-sensitive) → straight 5 m → helix → straight 5 m. " *
        "λ = 1550 nm.  M1: NEON 128-bit SIMD (2×f64/lane); StaticParticles " *
        "sweet spot ≈ N = 50–100.",
        out,
    )
    println("Wrote propagate_fiber benchmark to: ", out)
    return (path = out, desc = desc)
end

"""
    demo_benchmark_mcm_modify_propagate(; output_dir = …)

Benchmark `modify + Fiber + propagate_fiber` together — the full
per-sample pipeline — across Float64, Particles(2000), and
StaticParticles(50/100/200).
"""
function demo_benchmark_mcm_modify_propagate(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    desc = "MCM benchmark: modify + Fiber + propagate_fiber wall time across " *
           "Float64, Particles(2000), StaticParticles(50/100/200).  " *
           "The reference fiber includes spinning sensitivity in addition to " *
           "temperature dependence, and timing covers the full modify + " *
           "propagate pipeline."

    println("Running modify + propagate_fiber benchmarks …")
    results = _run_benchmarks(; scenario = :modify_propagate)

    out = joinpath(output_dir, "benchmark-mcm-modify-propagate.html")
    _bench_html(
        results,
        "MCM benchmark — modify + propagate_fiber (log scale)",
        "Same spinning + temperature-sensitive fiber as above.  Timing includes " *
        "modify() geometry rebuild, Fiber() construction, and propagate_fiber().",
        out,
    )
    println("Wrote modify + propagate_fiber benchmark to: ", out)
    return (path = out, desc = desc)
end

# =====================================================================
# Monolithic index entries
# =====================================================================

const DEMO3BENCHMARK_INDEX = [
    (group = "benchmarks", fn = demo_benchmark_mcm_propagate,         kwargs = (;)),
    (group = "benchmarks", fn = demo_benchmark_mcm_modify_propagate,  kwargs = (;)),
]

"""
    demo3benchmark_all(; index_output)

Run every demo in `DEMO3BENCHMARK_INDEX` and write `demo-index.html`.
"""
function demo3benchmark_entries()
    entries = Tuple{String, String, String, String}[]

    for d in DEMO3BENCHMARK_INDEX
        println("[ demo3benchmark ] $(d.fn)")
        result = d.fn(; d.kwargs...)
        desc = _demo_result_desc(result, d)
        for path in _demo_html_paths(result)
            push!(entries, (d.group, basename(path), path, desc))
        end
    end
    return entries
end

function demo3benchmark_all(; index_output::AbstractString = DEMO_MONOLITHIC_INDEX_OUTPUT)
    return _write_demo_index(
        [(title = "MCM propagation benchmarks",
          entries = demo3benchmark_entries(),
          group_titles = Dict("benchmarks" => "Speed benchmarks"))];
        index_output,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo3benchmark_all()
end
