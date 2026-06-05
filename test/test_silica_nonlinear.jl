"""
test_silica_nonlinear.jl — Automated regression tests for the `Nonlinear`
submodule (silica Raman + silica Brillouin).

CI-runnable companion to `test/human/test_nonlinear.ipynb`. The notebook
stays for visual judgement; this file fails the build if any of the physics
anchors below break.

# Test taxonomy
- T-PHYSICS  — compares model output to a published literature value.
- T-CONSISTENCY — bounds the mutual deviation of the three Raman models
  on robust summary statistics. A tight whole-spectrum residual would be
  fragile because the models genuinely disagree in the tail (`:bw`
  overestimates beyond 20 THz) and shoulder (`:hc` has the 15.2 THz
  feature). Instead, each spectrum is reduced to three Gaussian-like 
  scalars (peak centre, peak amplitude, FWHM).
- T-GUARDRAIL — runs the companion visual notebook
  (`test/human/test_nonlinear.ipynb`) end-to-end so a stale public-API
  reference in it fails the build instead of silently rotting.

# Brillouin scope
Even though one test point exercises GeO₂ doping (Niklès 1997 trend at
1319 nm), the bulk material is silica, so the suite lives in this file.
"""

using Test
using Statistics
using Printf
using Bifrost

const DEBUG_PLOTS = get(ENV, "BIFROST_TEST_DEBUG_PLOTS", "") in ("1", "true", "yes")

# Lazy-load Plots only when debug plotting was requested; CI should not need it.
const _PLOTS_OK = if DEBUG_PLOTS
    try
        @eval using Plots
        @eval gr()
        true
    catch err
        @warn "BIFROST_TEST_DEBUG_PLOTS set but Plots failed to load — debug plot skipped" exception=err
        false
    end
else
    false
end

# ── Helpers ───────────────────────────────────────────────────────────────

"""
    raman_summary_stats(freq_THz, g) -> (center, peak, width)

Reduce a Raman gain spectrum to three Gaussian-like summary statistics:
- `center`: argmax frequency, THz
- `peak`:   max amplitude
- `width`:  FWHM (linear-interpolated half-max crossings), THz

FWHM walks outward from the global max until `g` drops below `peak/2`,
so a secondary shoulder (e.g. HC at 15.2 THz) does not extend the
right-hand half-max past the main peak.
"""
function raman_summary_stats(freq_THz::AbstractVector, g::AbstractVector)
    i_peak = argmax(g)
    center = freq_THz[i_peak]
    peak_amp = g[i_peak]
    half = peak_amp / 2

    i_left = i_peak
    while i_left > 1 && g[i_left] >= half
        i_left -= 1
    end
    i_right = i_peak
    while i_right < length(g) && g[i_right] >= half
        i_right += 1
    end

    f_left = if i_left < i_peak && g[i_left+1] != g[i_left]
        freq_THz[i_left] + (half - g[i_left]) /
                           (g[i_left+1] - g[i_left]) *
                           (freq_THz[i_left+1] - freq_THz[i_left])
    else
        freq_THz[i_left]
    end
    f_right = if i_right > i_peak && g[i_right-1] != g[i_right]
        freq_THz[i_right-1] + (g[i_right-1] - half) /
                              (g[i_right-1] - g[i_right]) *
                              (freq_THz[i_right] - freq_THz[i_right-1])
    else
        freq_THz[i_right]
    end

    return (center = center, peak = peak_amp, width = f_right - f_left)
end

"SMF-28-like reference cross section, identical to the one in the notebook."
function _smf28_ref_fiber()
    return StepIndexCrossSection(
        SilicaGermaniaGlass(0.036),   # 3.6 mol % GeO₂ core
        SilicaGermaniaGlass(0.0),     # pure silica cladding
        8.2e-6, 125e-6;
        manufacturer = "Corning",
        model_number = "SMF-28-like",
    )
end

# ═════════════════════════════════════════════════════════════════════════
# Tests
# ═════════════════════════════════════════════════════════════════════════

@testset "silica_nonlinear" begin
    XS      = _smf28_ref_fiber()
    λ_ref   = 1550e-9
    T_ref   = 297.15
    A_eff   = effective_mode_area(XS, λ_ref, T_ref)
    γ_SMF   = (2π / λ_ref) * core_nonlinear_refractive_index(XS, λ_ref, T_ref) / A_eff
    m_GeO2  = XS.core_material.x_ge

    # Compute every Raman spectrum once; reuse for the consistency block
    # and the optional debug plot.
    freq_THz = collect(0.1:0.05:30)
    Ω        = 2π .* freq_THz .* 1e12
    models   = (:bw, :hc, :tabulated)
    spectra  = Dict{Symbol, Vector{Float64}}()
    stats    = Dict{Symbol, NamedTuple{(:center, :peak, :width)}}()
    for m in models
        g = Float64[Float64(g_R(o, γ_SMF; config = RamanConfig(model = m))) for o in Ω]
        spectra[m] = g
        stats[m]   = raman_summary_stats(freq_THz, g)
    end

    # ─── Raman T-PHYSICS — Consistency Table A (Raman peak at 13.2 THz) ─
    @testset "Raman T-PHYSICS  (peak at 13.2 THz Stokes shift)" begin
        # Each model's gain maximum must sit near the published 13.2 THz
        # silica Raman peak. rtol = 0.1 follows the reviewer's spec.
        for m in models
            @test stats[m].center ≈ 13.2 rtol = 0.1
        end
    end

    # ─── Raman T-CONSISTENCY — Gaussian-stat spread across the 3 models ─
    @testset "Raman T-CONSISTENCY  (Gaussian-stat spread across models)" begin
        centers = [stats[m].center for m in models]
        peaks   = [stats[m].peak   for m in models]
        widths  = [stats[m].width  for m in models]
        center_spread = (maximum(centers) - minimum(centers)) / mean(centers)
        peak_spread   = (maximum(peaks)   - minimum(peaks))   / mean(peaks)
        width_spread  = (maximum(widths)  - minimum(widths))  / mean(widths)

        # Peak frequency: all three models agree to < 5 %. Measured spread
        # is ~1 %; the bound exists to catch a model that drifts.
        @test center_spread < 0.05
        # Peak amplitude: BW slightly overestimates (Lorentzian tail), HC
        # slightly underestimates. Measured spread is ~14 %.
        @test peak_spread   < 0.20
        # FWHM: HC's shoulder narrows its FWHM, BW is broader. Measured
        # spread is ~30 %.
        @test width_spread  < 0.40
    end

    # ─── Brillouin T-PHYSICS — ν_B at 1550 nm SMF-28 ─────────────────────
    @testset "Brillouin T-PHYSICS  (ν_B near 11 GHz at 1550 nm SMF-28)" begin
        # Agrawal NLFO Ch. 9 reports ν_B ≈ 11.0 GHz for SMF-28 at 1550 nm.
        ν_B = brillouin_freq_shift(λ_ref; n_eff = 1.4447, m_GeO2 = m_GeO2) / 1e9
        @test ν_B ≈ 11.0 rtol = 0.05
    end

    # ─── Brillouin T-PHYSICS — Niklès 1997 doping trend at 1319 nm ──────
    @testset "Brillouin T-PHYSICS  (Niklès 1997 GeO₂ doping trend)" begin
        # Niklès, Thévenaz & Robert, JLT 15, 1842 (1997), Table II.
        # Code uses bulk acoustic velocity, so absolute values land ~15 %
        # above lit (documented in silica_brillouin.jl). The *trend* in
        # ν_B with doping is what this test pins down.
        λ_nik = 1319e-9
        n_eff = 1.4488
        ν_pure = Float64(brillouin_freq_shift(λ_nik; n_eff = n_eff, m_GeO2 = 0.0))   / 1e9
        ν_smf  = Float64(brillouin_freq_shift(λ_nik; n_eff = n_eff, m_GeO2 = 0.036)) / 1e9
        ν_dsf  = Float64(brillouin_freq_shift(λ_nik; n_eff = n_eff, m_GeO2 = 0.075)) / 1e9

        # Monotonic decrease with doping.
        @test ν_pure > ν_smf > ν_dsf

        # Magnitude: Niklès Δν_B(0 → 7.5 %) = 11.32 − 10.61 = 0.71 GHz.
        # Model gives ~0.69 GHz — agreement within 15 %.
        Δ_model = ν_pure - ν_dsf
        Δ_lit   = 11.32 - 10.61
        @test Δ_model ≈ Δ_lit rtol = 0.15
    end

    # ─── Brillouin T-PHYSICS — Kobyakov 2010 threshold scale + ordering ─
    @testset "Brillouin T-PHYSICS  (SBS threshold ordering and scale)" begin
        # Kobyakov 2010 Table 3 SMF-28 measurements:
        #   25 km ≈ 5.5 mW, 50 km ≈ 4.8 mW, 100 km ≈ 4.4 mW.
        # P_th ∝ 1/g_B, so model predictions land 30–40 % below Kobyakov.
        # This test pins the *order of magnitude* and the L-dependence.
        g_B = g_B_peak_GeO2(m_GeO2)
        # Pass α explicitly to dodge the SMF-28 empirical-loss warning that
        # would otherwise fire on every test run.
        α_1550 = 0.20 * log(10) / (10 * 1e3)    # 0.20 dB/km → 1/m
        Ps = [brillouin_threshold(A_eff, L * 1e3;
                                  λ_pump   = λ_ref,
                                  g_B_peak = g_B,
                                  α        = α_1550).P_threshold_mW
              for L in (25.0, 50.0, 100.0)]

        # All positive, finite, on the order of a few mW.
        @test all(p -> 1.0 < p < 50.0, Ps)
        # Monotone decrease with length (effective length saturation).
        @test Ps[1] > Ps[2] > Ps[3]
        # Within a factor of 2 of each Kobyakov measurement.
        for (p_model, p_lit) in zip(Ps, (5.5, 4.8, 4.4))
            @test 0.5 < p_model / p_lit < 2.0
        end
    end

    # ─── Optional debug plot ────────────────────────────────────────────
    if DEBUG_PLOTS && _PLOTS_OK
        try
            plt = Main.plot(
                xlabel  = "Frequency shift Ω/2π  (THz)",
                ylabel  = "g_R  (×10⁻³ W⁻¹ m⁻¹)",
                title   = "Raman gain — measured summary stats (debug)",
                size    = (900, 500),
                legend  = :topright,
            )
            colors = Dict(:bw => :firebrick,
                          :hc => :seagreen,
                          :tabulated => :dodgerblue)
            for m in models
                Main.plot!(plt, freq_THz, spectra[m] .* 1e3;
                           label = string(m), lw = 2, color = colors[m])
                Main.vline!(plt, [stats[m].center];
                            color = colors[m], ls = :dash, alpha = 0.6,
                            label = "")
            end
            Main.vline!(plt, [13.2]; ls = :dot, lc = :gray,
                        label = "13.2 THz (published)")
            Main.vline!(plt, [15.2]; ls = :dot, lc = :gray, alpha = 0.5,
                        label = "15.2 THz (HC shoulder)")
            outdir = joinpath(dirname(@__DIR__), "output", "test")
            mkpath(outdir)
            outfile = joinpath(outdir, "debug_raman_summary_stats.png")
            Main.savefig(plt, outfile)
            @info "Debug plot written" file=outfile
            for m in models
                @info @sprintf("%-10s  center=%.2f THz   peak=%.3e   FWHM=%.2f THz",
                               string(m), stats[m].center,
                               stats[m].peak, stats[m].width)
            end
        catch err
            @warn "Debug plot failed to render" exception=err
        end
    end
end

# ═════════════════════════════════════════════════════════════════════════
# Companion-notebook guardrail
# ═════════════════════════════════════════════════════════════════════════
#
# The visual notebook test/human/test_nonlinear.ipynb is not run by CI on its
# own, so a rename or removal in the public API can leave it calling a symbol
# that no longer exists (exactly what happened after the cross-section layer was
# split). The guardrail below extracts the notebook's Julia code cells and runs
# them in a throwaway module, asserting the whole notebook executes without
# error. Plotting verbs are stubbed so it needs no Plots dependency, matching
# this suite's no-Plots-in-CI policy; the guard is on the Bifrost API calls the
# notebook makes, which is where the regression risk lives.

const _JSON_ESCAPES = Dict('n' => '\n', 't' => '\t', 'r' => '\r', 'b' => '\b',
                           'f' => '\f', '"' => '"', '\\' => '\\', '/' => '/')

"""
    _read_json_string(text, q) -> (String, Int)

Read one JSON string from `text` whose opening quote is at index `q`. Returns the
unescaped contents and the index just past the closing quote.
"""
function _read_json_string(text::String, q::Int)
    buf = IOBuffer()
    i = nextind(text, q)
    while true
        c = text[i]
        if c == '"'
            return String(take!(buf)), nextind(text, i)
        elseif c == '\\'
            i = nextind(text, i)
            e = text[i]
            if e == 'u'
                hex = IOBuffer()
                for _ in 1:4
                    i = nextind(text, i)
                    write(hex, text[i])
                end
                write(buf, Char(parse(UInt16, String(take!(hex)); base = 16)))
            else
                write(buf, get(_JSON_ESCAPES, e, e))
            end
        else
            write(buf, c)
        end
        i = nextind(text, i)
    end
end

"""
    _read_json_string_array(text, lb) -> (String, Int)

Read a JSON array of strings whose opening bracket is at index `lb`. Returns the
concatenated contents and the index just past the closing bracket.
"""
function _read_json_string_array(text::String, lb::Int)
    parts = String[]
    i = nextind(text, lb)
    while true
        c = text[i]
        if c == '"'
            s, i = _read_json_string(text, i)
            push!(parts, s)
        elseif c == ']'
            return join(parts), nextind(text, i)
        else
            i = nextind(text, i)
        end
    end
end

"""
    _notebook_code_source(path) -> String

Concatenate the Julia source of every code cell in the `.ipynb` at `path`. Relies
only on the JSON invariant that an unescaped double quote never appears inside a
string, so `"cell_type"` and `"source"` reliably mark object keys.
"""
function _notebook_code_source(path::AbstractString)
    text = read(path, String)
    cells = String[]
    pos = firstindex(text)
    while true
        ct = findnext("\"cell_type\"", text, pos)
        ct === nothing && break
        vq = findnext('"', text, nextind(text, last(ct)))
        kind, after = _read_json_string(text, vq)
        if kind == "code"
            sk = findnext("\"source\"", text, after)
            sk === nothing && break
            lb = findnext('[', text, nextind(text, last(sk)))
            src, after = _read_json_string_array(text, lb)
            push!(cells, src)
        end
        pos = after
    end
    return join(cells, "\n\n")
end

@testset "test_nonlinear.ipynb  (T-GUARDRAIL: notebook executes)" begin
    notebook = joinpath(@__DIR__, "human", "test_nonlinear.ipynb")
    @test isfile(notebook)

    code = _notebook_code_source(notebook)
    # Extraction sanity: first and last code cells must be present.
    @test occursin("StepIndexCrossSection", code)
    @test occursin("Relative-shift error", code)

    # Drop the Plots import (verbs stubbed below) and neutralise the one display
    # call so the notebook runs headless without a Plots dependency.
    code = replace(code, "using Plots, Printf, Statistics" => "using Printf, Statistics")
    code = replace(code, "display(" => "identity(")
    @test !occursin("using Plots", code)

    sandbox = Module(:NotebookGuardSandbox)
    Core.eval(sandbox, quote
        _nbstub(args...; kwargs...) = nothing
        const plot      = _nbstub
        const plot!     = _nbstub
        const vline!    = _nbstub
        const vspan!    = _nbstub
        const scatter!  = _nbstub
        const annotate! = _nbstub
        const gr        = _nbstub
        const Plots     = (mm = 1.0,)   # only Plots.mm is referenced
    end)

    ran = try
        redirect_stdio(stdout = devnull, stderr = devnull) do
            Base.include_string(sandbox, code, "test_nonlinear.ipynb")
        end
        true
    catch err
        bt = catch_backtrace()
        @error "test_nonlinear.ipynb failed to execute" exception = (err, bt)
        false
    end
    @test ran
end
