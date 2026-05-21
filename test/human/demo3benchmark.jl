# =====================================================================
# demo3benchmark.jl — MCM propagation speed benchmarks
# =====================================================================
#
# Extends the MCM temperature-PTF scenarios from demo4mcm.jl with
# timing measurements across four particle types:
#
#   * Float64          — deterministic baseline
#   * Particles(40)    — heap-allocated, arbitrary N, standard MCM
#   * StaticParticles(20/40) — SVector-backed, SIMD-friendly
#
# For each variant two times are measured:
#   * first-call  : wall time of the very first call, i.e. JIT
#                   compilation overhead + one execution.
#   * second-call : wall time of one repeated call after JIT.
#
# M1 SIMD note: Apple M1 NEON has 128-bit lanes (2×f64/lane).
# StaticParticles shines when N is small enough for Julia to emit
# SIMD-vectorised loops over the SVector. These human-demo benchmarks use
# small ensembles so the demo remains runnable during ordinary verification.
#
# Two HTML files are written:
#   benchmark-mcm-propagate.html  — bar chart + table for propagate_fiber
#   benchmark-mcm-modify.html     — bar chart + table for modify + propagate_fiber
#
# `demo3benchmark_all()` runs both and writes `output/demo3benchmark.html`.

if !isdefined(Main, :_mcm_demo_fiber)
    include(joinpath(@__DIR__, "demo4mcm.jl"))
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
    fiber  = _mcm_demo_fiber(ΔT_K)
    mpath  = modify(fiber)
    return Fiber(mpath; cross_section = _MCM_DEMO_XS, T_ref_K = T_K)
end

# Time the first call to `f()` (JIT + one run) in milliseconds.
function _first_call_ms(f)
    t0 = time_ns()
    f()
    return (time_ns() - t0) / 1e6
end

# Second-call wall time in milliseconds. This keeps the human demo quick.
function _steady_ms(f)
    return @elapsed(f()) * 1e3
end

# =====================================================================
# Benchmark cases
# =====================================================================

# Each case: (label, N_particles_or_0, make_T_fn)
# make_T returns a fresh value each time it is called.
function _bench_cases()
    return [
        ("Float64", 0, () -> Float64(_MCM_DEMO_T_NOM_C)),
        ("Particles(40)", _MCM_DEMO_PARTICLES_N, () -> Particles(
            _MCM_DEMO_PARTICLES_N,
            Normal(Float64(_MCM_DEMO_T_NOM_C), Float64(_MCM_DEMO_T_SIG_C)),
        )),
        ("StaticParticles(20)", _MCM_DEMO_STATIC_N, () -> StaticParticles(
            _MCM_DEMO_STATIC_N,
            Normal(Float64(_MCM_DEMO_T_NOM_C), Float64(_MCM_DEMO_T_SIG_C)),
        )),
        ("StaticParticles(40)", 40, () -> StaticParticles(
            40,
            Normal(Float64(_MCM_DEMO_T_NOM_C), Float64(_MCM_DEMO_T_SIG_C)),
        )),
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
            fn = () -> propagate_fiber(
                fiber;
                λ_m = _MCM_DEMO_λ_M,
                rtol = _MCM_DEMO_RTOL,
                verbose = false,
            )
        else
            # Benchmark: modify + propagate_fiber together
            T_val  = make_T()
            T_K    = T_val + 273.15
            ΔT_K   = T_K - _MCM_DEMO_T_REF_K
            base_f = _mcm_demo_fiber(ΔT_K)
            fn     = () -> begin
                mp = modify(base_f)
                f2 = Fiber(mp; cross_section = _MCM_DEMO_XS, T_ref_K = T_K)
                propagate_fiber(
                    f2;
                    λ_m = _MCM_DEMO_λ_M,
                    rtol = _MCM_DEMO_RTOL,
                    verbose = false,
                )
            end
        end

        fc  = _first_call_ms(fn)   # includes JIT
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

Benchmark `propagate_fiber` alone across Float64, Particles(40),
and StaticParticles(20/40). The fiber is pre-built; only the
ODE propagation is timed.
"""
function demo_benchmark_mcm_propagate(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    desc = "MCM benchmark: propagate_fiber wall time across Float64, " *
           "Particles(40), StaticParticles(20/40).  " *
           "The reference fiber now includes twist sensitivity in addition to " *
           "temperature dependence. First-call includes JIT; second-call is " *
           "one repeated call."

    println("Running propagate_fiber benchmarks …")
    results = _run_benchmarks(; scenario = :propagate)

    out = joinpath(output_dir, "benchmark-mcm-propagate.html")
    _bench_html(
        results,
        "MCM benchmark — propagate_fiber (log scale)",
        "Fiber: straight 5 m → helix (D=0.05 m, ~10002 turns, twist + " *
        "temperature-sensitive) → straight 5 m → helix → straight 5 m. " *
        "λ = 1550 nm. StaticParticles cases use small N so this human demo " *
        "runs quickly.",
        out,
    )
    println("Wrote propagate_fiber benchmark to: ", out)
    return (path = out, desc = desc)
end

"""
    demo_benchmark_mcm_modify_propagate(; output_dir = …)

Benchmark `modify + Fiber + propagate_fiber` together — the full
per-sample pipeline — across Float64, Particles(40), and StaticParticles(20/40).
"""
function demo_benchmark_mcm_modify_propagate(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    desc = "MCM benchmark: modify + Fiber + propagate_fiber wall time across " *
           "Float64, Particles(40), StaticParticles(20/40).  " *
           "The reference fiber now includes twist sensitivity in addition to " *
           "temperature dependence, and timing covers the full modify + " *
           "propagate pipeline."

    println("Running modify + propagate_fiber benchmarks …")
    results = _run_benchmarks(; scenario = :modify_propagate)

    out = joinpath(output_dir, "benchmark-mcm-modify-propagate.html")
    _bench_html(
        results,
        "MCM benchmark — modify + propagate_fiber (log scale)",
        "Same twist + temperature-sensitive fiber as above.  Timing includes " *
        "modify() geometry rebuild, Fiber() construction, and propagate_fiber().",
        out,
    )
    println("Wrote modify + propagate_fiber benchmark to: ", out)
    return (path = out, desc = desc)
end

# =====================================================================
# Index page (demo3benchmark.html)
# =====================================================================

const DEMO3BENCHMARK_INDEX = [
    (group = "benchmarks", fn = demo_benchmark_mcm_propagate,         kwargs = (;)),
    (group = "benchmarks", fn = demo_benchmark_mcm_modify_propagate,  kwargs = (;)),
]

"""
    demo3benchmark_all(; index_output)

Run every demo in `DEMO3BENCHMARK_INDEX` and write `demo3benchmark.html`.
"""
function demo3benchmark_all(;
    index_output::AbstractString = joinpath(@__DIR__, "..", "..", "output", "demo3benchmark.html"),
)
    entries = Tuple{String, String, String, String}[]

    for d in DEMO3BENCHMARK_INDEX
        println("[ demo3benchmark ] $(d.fn)")
        result = d.fn(; d.kwargs...)
        desc_inline = (result isa NamedTuple && haskey(result, :desc)) ?
                      String(result.desc) : ""
        desc_entry  = hasproperty(d, :desc) ? d.desc : ""
        desc        = isempty(desc_inline) ? desc_entry : desc_inline

        paths = result isa NamedTuple ? values(result) : (result,)
        for v in paths
            if v isa AbstractString && endswith(v, ".html")
                push!(entries, (d.group, basename(v), v, desc))
            end
        end
    end

    open(index_output, "w") do io
        println(io, """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>BIFROST MCM benchmarks</title>
  <style>
    body { font-family: sans-serif; max-width: 800px; margin: 2em auto; background: #111; color: #ddd; }
    h1   { font-size: 1.5em; border-bottom: 1px solid #444; padding-bottom: 0.3em; }
    h2   { font-size: 1.15em; margin-top: 1.8em; color: #4db87a; }
    ul   { padding-left: 1.2em; }
    li   { margin: 1em 0; }
    a    { font-weight: bold; color: #4db87a; }
    p.desc { margin: 0.3em 0 0 0; color: #999; font-size: 0.95em; }
    nav.index-nav { font-size: 0.85em; margin-bottom: 1em; color: #666; }
    nav.index-nav a { font-weight: normal; color: #4db87a; margin-right: 0.8em; }
  </style>
</head>
<body>
  <nav class="index-nav">
    <a href="demo1.html">demo1</a>
    <a href="demo2.html">demo2</a>
    <a href="demo4mcm.html">demo4mcm</a>
    <a href="demo3benchmark.html">demo3benchmark</a>
  </nav>
  <h1>BIFROST MCM propagation benchmarks</h1>""")

        seen_groups = String[]
        for (g, _, _, _) in entries
            g in seen_groups || push!(seen_groups, g)
        end
        for g in seen_groups
            println(io, "  <h2>Speed benchmarks</h2>")
            println(io, "  <ul>")
            for (eg, title, path, desc) in entries
                eg == g || continue
                println(io, "    <li>")
                println(io, "      <a href=\"$(path)\">$(title)</a>")
                println(io, "      <p class=\"desc\">$(desc)</p>")
                println(io, "    </li>")
            end
            println(io, "  </ul>")
        end
        println(io, """</body>
</html>""")
    end

    println("Wrote demo3benchmark index to: ", index_output)
    return index_output
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo3benchmark_all()
end
