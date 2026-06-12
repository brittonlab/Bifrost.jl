# Demo Intent

This file records the intent of every visual demo that historically lived in
`test/human/demo*.jl`, distilled by running all demos at HEAD, inspecting the generated
`output/*.html` artifacts (including headless-browser snapshots), and reading the demo
sources. It is the authoritative narrative from which the consolidated notebook
`bifrost-demos.ipynb` and its support file `demo-helper.jl` are authored. Each class
below tells the story of what the demos teach, which edge case or invariant they
codify, and what a reader should look for in the output.

Bifrost API calls appear only where the precise geometry is the point of the demo.

## Class index

| # | Class | Legacy source | Atomic demos |
| --- | --- | --- | --- |
| 0 | Smallest example | `demo-smallest.jl` | 1 |
| 1 | Path geometry | `demo-path-geometry.jl` | 7 |
| 2 | Modify (field-level MCM meta) | `demo1.jl` | 12 |
| 3 | Material spin | `demo1.jl` | 1 |
| 4 | Adaptive step-doubling | `demo1.jl` | 1 |
| 5 | JumpBy / JumpTo — 2D sweeps | `demo2.jl` | 5 |
| 6 | Meta × JumpTo interplay — 2D | `demo2.jl` | 3 |
| 7 | JumpBy / JumpTo — 3D scenes | `demo2.jl` | 12 |
| 8 | MCM temperature PTF | `demo3mcm.jl` | 2 |
| 9 | MCM speed benchmarks | `demo3benchmark.jl` | 2 |

## Class 0 — Smallest example

The minimum ceremony from nothing to a propagated Jones matrix:

```julia
xs = StepIndexCrossSection(SilicaGermaniaGlass(0.036), SilicaGermaniaGlass(0.0),
                           8.2e-6, 125e-6)
sb = SubpathBuilder(); start!(sb)
straight!(sb; length = 0.5); bend!(sb; radius = 0.05, angle = π/2)
straight!(sb; length = 0.5); seal!(sb)
fiber = Fiber(build(sb); cross_section = xs, T_ref_K = 297.15)
J, stats = propagate_fiber(fiber; λ_m = 1550e-9)
```

Lesson: the full layer stack (material → cross-section → geometry → fiber →
propagation) in ~10 lines. The 90° bend at R = 0.05 m introduces bend birefringence
with a fixed axis, so `J` comes out diagonal `diag(e^{-iφ}, e^{+iφ})` with
φ ≈ 0.033 rad — a pure retarder, |det J| = 1 (lossless SU(2)). `stats` shows the
propagation decomposed into 3 intervals: the propagator never steps across the
straight/bend segment boundaries (breakpoint invariant).

## Class 1 — Path geometry (geometry layer only, no optics)

Human-in-the-loop visual debugging for `path-geometry*.jl`. Every demo builds a path
and renders the interactive 3D inspector (see Visual techniques, §V1).

- **simple** — one Subpath: straight · 90° bend · straight · catenary (`a = 0.03`) ·
  60° bend · straight, sealed at the natural exit. Lesson: all four curve primitives
  compose with tangent (G1) continuity at every join; prints `path_length` and
  `writhe` (0 for this planar path) as scalar geometry queries.
- **segment_labels** — same spirit plus a helix, every segment authored with a
  `Nickname` meta (`lead-in`, `90° bend`, `spacer`, `sag`, `spin section`,
  `lead-out`). Lesson: `Nickname` meta flows from authoring through the build to
  presentation; names live with the segment definition.
- **helix_0 / helix_pi_3 / helix_2pi_3** — straight · 2-turn helix
  (R = 0.03, pitch = 0.02) · straight, with `axis_angle` ∈ {0, π/3, 2π/3}.
  Lesson/edge case: `axis_angle` rotates the helix axis about the incoming tangent;
  entry stays tangent-continuous for every value, the exit direction (and where the
  trailing straight goes) changes, and the arc length is invariant (0.4791 m in all
  three). Comparing the three outputs side by side is the point.
- **jumps_min_radius** — the “paddle”: a `PathBuilt` of 5 Subpaths. Vertical 1 m
  straights at x = 0, 1, 2, 3 alternate up/down and are joined by terminal `jumpto!`
  connectors with anti-parallel incoming tangents (180° U-turns), each with its own
  `min_bend_radius` (0.4, 0.1, 0.05). Subpath 4 also carries an *interior* `jumpby!`
  (local-frame `delta = (-1, 0, 0)`, reversing tangent) before its sealing `jumpto!`.
  Subpaths 2–5 use `start!(sb, :inherit)` so start states flow from the previous
  subpath. Lessons: PathBuilt assembly under a shared global arc length; `:inherit`;
  interior vs terminal connectors; per-connector `min_bend_radius`; G2 quintic
  U-turns.
- **pathbuilt** — three Subpaths (straight sealed by `jumpto!` at (0,0,0.2);
  `:inherit` 90° bend sealed by `jumpto!` at its analytic exit; `:inherit` helix
  sealed with `seal!`). Lesson: a `jumpto!` seal pins exact handoff coordinates, the
  conformity check validates each boundary, and per-Subpath `Nickname` meta labels
  whole subpaths.

## Class 2 — Modify (field-level MCM meta on segment fields)

Twelve structurally identical experiments. Baseline geometry is a 3-segment
inverted-U (straight L = 1 · π-bend R = 0.5 · straight L = 1) or, for helix targets, a
4-segment variant with a helix (R = 0.15, pitch = 0.25, turns = 1.5) inserted after
the bend. Exactly one segment carries one `MCMadd(:field, δ)` or `MCMmul(:field, f)`
with a plain `Float64` perturbation — so the “MCM” machinery is exercised with
deterministic offsets and the geometric meaning of each field is visible to the eye.
`Fiber(builder; ...)` applies the meta during its single build.

Variants render side by side (Visual techniques §V2): perturbed segment red, others
green, terminal connector faint gray.

| Demo | Target field | Variants | What it shows |
| --- | --- | --- | --- |
| straight :length (add) | first straight | −0.4 | downstream translates rigidly |
| bend :radius (add) | π-bend | −0.25, +0.50 | U widens/narrows, exit shifts |
| bend :angle (add) | π-bend | −π/2, +π/4, +π | +π closes a full circle |
| straight :length (mul) | first straight | ×(−0.4), ×0.5 | ×negative walks backward |
| bend :radius (mul) | π-bend | ×0.5, ×2.0 | scale vs offset semantics |
| bend :angle (mul) | π-bend | ×0.5, ×1.25 | — |
| helix :radius (add/mul) | helix | ±δ, ×0.5/×2 | coil fattens, pitch fixed |
| helix :pitch (add/mul) | helix | ±δ, ×0.5/×2 | coil stretches axially |
| helix :turns (add/mul) | helix | ±0.5, ×0.67/×1.5 | exit phase/direction changes |

Lessons: which geometry field each knob actually perturbs; add vs mul semantics
(including the sign-flip edge case `MCMmul(:length, -0.4)` → backward straight);
without an anchor, everything downstream of the perturbed segment rides along
rigidly — the contrast with Class 6’s anchored paths.

## Class 3 — Material spin

Geometry-only: `start!(sb; spin_rate = 2π)` on straight (1 m) · helix (R = 0.5,
pitch = 0.05, 4 turns) · straight (1 m), sealed at the natural exit. In the 3D
inspector the status box reads τ_spin = 6.2832 rad/m everywhere, ∫τ_spin ds
accumulates linearly with s, and the red spin arrow rotates in the normal–binormal
plane as the cursor moves — while T̂/N̂/B̂ follow the geometric frame. Lesson:
material spin set at `start!` is resolved by the geometry layer into path-coordinate
spin, carried independently of the helix’s geometric torsion.

## Class 4 — Adaptive step-doubling diagnostic

A solver diagnostic on a smooth, noncommuting generator (not a fiber):

    K(s) = α·i·σx·cos(πs) + β·i·σz·sin(2πs),  α = 1.2, β = 0.9, s ∈ [0, 2]

`collect_adaptive_steps` (a recording mirror of the production controller in
`path-integral.jl`; the solver itself is untouched) returns every accepted and
rejected trial; rendering is the three-panel diagnostic of §V7. With rtol = 1e-6 the
run accepts ~370 steps and rejects ~15. What to look for: step size locks to the local
error budget, growing where ‖K(s)‖ dips; each over-grown step triggers a
rejection cascade (red ✕ column) right after the ‖K‖ minima; the err/tol panel hugs
the threshold from below (controller efficiency); accepted-step err/tol never exceeds
1. Lesson: the controller is working as designed — effort concentrates where the
generator varies fastest.

## Class 5 — JumpBy / JumpTo, 2D sweeps (x–z plane)

All paths lie in y = 0; rendering is the small-multiples SVG row of §V4. Each demo
sweeps exactly one connector degree of freedom across an otherwise identical
scene — straight (1 m) · connector · straight (1 m), connector drawn red.

- **jumpby tangent_out** — `jumpby!(delta = (0.4, 0, 0.4), tangent = t)` for
  t ∈ {(+1,0,1)/√2, (0,0,1), (−1,0,1)/√2}. Lesson: a JumpBy’s outgoing tangent is
  expressed in the **local** frame; the trailing straight exits in that direction.
- **jumpto point** — `jumpto!(point = (x, 0, 1.5))` for x ∈ {0, 0.3, 0.6, 1.0}.
  Lesson: a JumpTo target is a **global** waypoint; growing transverse offset bends
  the connector harder.
- **jumpto tangent_global** — fixed point, `incoming_tangent` ∈ {+x, +z, (−1,0,1)/√2}.
  Lesson: the landing direction is a global-frame constraint; the connector arrives
  at the same point from three different headings, and the trailing straight leaves
  along each.
- **jumpto routing** — three successive `jumpto!`s with anti-parallel incoming
  tangents (each (0,0,∓1) against the previous exit) building a serpentine.
  Lesson: composite routing through waypoints; each `jumpto!` seals a Subpath, so the
  result is multi-Subpath under the hood.
- **jumpto min_radius** — the canonical infeasibility case (transverse unit chord,
  anti-parallel tangents): sweep `min_bend_radius` ∈ {0.10, 0.30, 0.49, 0.51, 0.70}.
  Feasible variants build; infeasible ones throw at build time, are trapped, and
  render as partial paths labeled “(infeasible — n/3 built)” (§V5). Lesson:
  `min_bend_radius` is a hard build-time guardrail, and the sweep locates the
  empirical threshold. *Observed at HEAD:* 0.10 and 0.30 build; 0.49, 0.51, 0.70 are
  all infeasible — i.e. the real threshold lies between 0.3 and 0.49, not at the
  ≈ 0.5 m the legacy docstring claimed (0.5 = chord/2 is the circular-arc ideal; the
  quintic connector’s peak curvature is necessarily higher).

## Class 6 — Meta × JumpTo interplay (2D overlays)

Three demos that answer: *when geometry upstream is perturbed, what happens to
everything downstream?* Rendered as baseline-vs-modified overlays (§V3) with total
path lengths and Δ% in the legend.

1. **jumpby drift** — S-curve (two opposed 90° bends) · `jumpby!(delta = (0,0,0.8))` ·
   straight, with `MCMmul(:radius, 1.25)` on both bends. The JumpBy delta is
   fiber-relative, so there is no anchor: the whole downstream trajectory drifts and
   the endpoint visibly separates (~10% length growth).
2. **jumpto anchor** — same S-curve, but sealed by `jumpto!` to a fixed lab-frame
   point (placed so total baseline length is 4.000 m). With `MCMmul(:radius, 1.5)`
   the interior swings wide, yet the endpoint stays pinned — the terminal connector
   absorbs the slack.
3. **jumpto anchor + thermal** — straight (0.5 m) and its sealing `jumpto!` *both*
   carry `MCMadd(:T_K, ΔT)` with ΔT = 0.05/α_lin (a 5% thermal expansion). The
   terminal connector thermally expands — its arc length scales by τ — while still
   landing at the fixed point (the length-constrained connector resolve). The extra
   length shows up as visible connector curvature.

Lesson trio: unanchored perturbations propagate downstream; a JumpTo anchor confines
them; `:T_K` on a `jumpto!` seal stretches the connector itself without moving the
anchor. (`:T_K` interpretation happens in `Fiber`, never in the geometry layer.)

## Class 7 — JumpBy / JumpTo, 3D scenes

The 2D lessons restaged in 3D (§V2 rendering: connector red, fixed segments green,
variants offset along x), plus 3D-only G2-continuity content. Scaffold is again
straight (1 m) · connector · straight (1 m) unless noted.

JumpBy group (local frame):

- **delta axial** — `delta = (0,0,d)`, d ∈ {0.3, 0.6, 1.0}: straight-through gap.
- **delta transverse** — `delta = (d,0,0.5)`, d ∈ {0, 0.2, 0.5}: lateral dogleg.
- **tangent_out** — 3D version of the 2D sweep.
- **curvature_out** — `curvature_out` ∈ {0, (0,±2,0)} with fixed diagonal tangent:
  the G2 exit knob; the connector leaves with prescribed curvature, bowing out of
  plane.
- **after_bend** — a 90° bend precedes the JumpBy; `delta` ∈ {(0,0,0.5), (±0.3,0,0.5)}.
  Lesson: the delta lives in the **rotated** local frame the bend left behind.
- **g2_inheritance** — bend (R ∈ {0.30, 0.50, 0.80}, π/4) directly into
  `jumpby!(delta = (0.3,0,0.3))`. Lesson: with no explicit `curvature_out` /
  incoming spec, the connector inherits the bend’s exit curvature κ = 1/R (G2
  continuity by default); smaller R visibly launches the connector on a tighter arc.

JumpTo group (global frame):

- **point (axial)** — `point = (0,0,z)`, z ∈ {1.3, 1.6, 2.0}.
- **point (transverse)** — `point = (x,0,1.5)`, x ∈ {0, 0.3, 0.6}.
- **tangent_global** — `incoming_tangent` ∈ {+x, +z, (−1,0,1)/√2} at a fixed point.
- **curvature_global** — `incoming_curvature` ∈ {0, (±10,0,0)}: the G2 landing knob;
  the connector arrives flat or pre-curled.
- **after_bend** — 90° bend then `jumpto!` to global destinations. Lesson: unlike
  JumpBy-after-bend, the target does not rotate with the local frame.
- **routing** — the serpentine composite in 3D.

## Class 8 — MCM temperature PTF (end-to-end uncertainty)

The flagship physics demo. Fiber (cross-section: 3.6% Ge-doped core / pure-silica
cladding, a = 8.2 µm, λ = 1550 nm, T_ref = 303.15 K):

    straight 5 m
    helix R = 0.025 m, pitch = 0.05 m, turns = 10001.892…   ← MCMadd(:T_K, ΔT)
    straight 5 m
    helix R = 0.025 m, pitch = 0.05 m, turns = 10000        (reference)
    straight 5 m

with a sinusoidal material spin `spin_rate = s -> sin(2πs/100)` over the whole
Subpath. Physics: bend birefringence scales as 1/R², so at R = 2.5 cm the first helix
accumulates |Δβ|·L ≈ 1775·2π rad of retardation; the fractional turn count is tuned
so that at 30 °C it sits exactly at **mid-fringe** (mod(Γ, 2π) = π), the operating
point of maximum temperature sensitivity (at a crossing, dPTF/dT = 0). ΔT enters
twice, deliberately: as `:T_K` meta on the helix (geometric thermal expansion via the
cladding CTE, baked in by `Fiber`) and as `T_ref_K = T_K_particles` on a second
`Fiber` binding (temperature-dependent material indices). Over ±5 °C the retardation
swings by ≈ ±0.75π.

- **ptf** — `StaticParticles(50)`, T ~ N(30 °C, 5 °C), **one** propagation. Each
  particle’s Jones matrix is extracted and converted to Stokes parameters for input
  state H. Output: polarization angle vs T, and S1/S2/S3/DLP vs T (§V6). The output
  stays nearly linearly polarized (DLP ≈ 1, S3 ≈ 0): the helix birefringence rotates
  the linear polarization angle without adding ellipticity, so the sensitive
  observable is the **angle**, not DLP.
- **scatter** — `Particles(500)` version (the legacy demo used 2000; reduced so the
  notebook executes in minutes — the ensemble story is unchanged); angle vs T plus
  the Poincaré equatorial projection (S1–S2), where the ensemble traces an arc
  parameterized by temperature.

Lessons: a full MCM pipeline (uncertain T → geometry meta + material binding →
ensemble propagation in one pass → per-particle post-processing); mid-fringe design;
why ensembles beat nominal-value runs for PTF questions.

## Class 9 — MCM speed benchmarks

Same fiber as Class 8. Cases: `Float64`, `Particles(N)`, `StaticParticles(50)`,
`StaticParticles(75)`. For each, two timings: **first-call** (JIT + one run, via
`time_ns`) and **steady-state** (a warm re-run). Two scenarios: `propagate_fiber`
alone on a pre-built fiber, and build + propagate (the full per-sample pipeline
including `:T_K` thermal geometry). Output: log-scale grouped bars plus a ratio
table (§V8). Lessons: the ensemble overhead of MCM relative to Float64;
`StaticParticles` SIMD sweet spot at small N; JIT cost vs marginal cost. Numbers are
environment-dependent — the demo is informational, never asserted.

Notebook scaling: the cells are guarded by a `RUN_BENCHMARKS` flag (default `false`)
so casual runs skip them; the stored outputs come from one real run with
`Particles(500)` standing in for the legacy `Particles(2000)` and a single warm
timing instead of BenchmarkTools' repeated minimum (the legacy methodology, noted in
the cell, takes hours on modest hardware).

## Legacy orchestration (replaced by notebook structure)

`demo-index.jl` + `demo-index-helpers.jl` maintained per-file registries of
`(group, fn, desc)` entries; each demo returned its HTML path(s) and optional inline
description, and a monolithic `output/demo-index.html` grouped them under section and
group headings. In the notebook this machinery disappears: sections are the classes
above, the markdown cell before each demo carries its description, and plots render
inline.

## V — Visual presentation techniques

The catalog below is sufficient to reproduce the visual storytelling without the
legacy source.

**V1. Interactive 3D path inspector** (library: `write_path_geometry_plot3d`, kept).
One path in a 3D scene, equal aspect (`aspectmode: data`): faint gray centerline;
dot markers along the path color-graded by arc length (light → dark); open-circle
markers at segment boundaries; green start dot, red end dot; segment `Nickname`
labels floated as 3D text. A cursor driven by horizontal mouse position scrubs arc
length: at the cursor sit a black dot, a translucent square spanning the local
normal–binormal plane, a short axis triad T̂ (orange), N̂ (blue), B̂ (green), and a
red arrow in the N–B plane pointing along the accumulated spin phase ∫τ_spin ds. A
status box updates live with s; x, y, z; κ; τ_geom; τ_spin; ∫τ_spin ds (rad and
deg). A `fidelity` knob scales sampling density. A footer explains the mouse
interaction.

**V2. Side-by-side variant row (3D)** — the workhorse comparison layout. N variants
of one scenario, differing in exactly one parameter, each built independently and
offset along +x by a `variant_spacing`; per-variant floating text label above the
scene; semantic line coloring: red = the element under study (connector or perturbed
segment), green = fixed context segments, faint gray = terminal connector when it is
incidental; green start dots; dark background; camera looking from −y. Driven
parametrically by a list of `(label, build_fn)` pairs — this is the parametric
tooling that makes the sweeps cheap to author.

**V3. Baseline-vs-modified overlay (2D, x–z)** — two paths drawn in the same axes:
baseline black, modified red, on a light background with an equal-aspect mapping,
padded bounds, and a unit grid. Open circles mark the start and end of every placed
segment (shared boundary markers overlay). A legend reports L_baseline, L_modified,
and Δ% — the quantitative readout of thermal/parametric growth. Legend corner is
selectable so it never covers the interesting part.

**V4. Small-multiples sweep (2D, x–z)** — like V2 but flat: per-variant x-offset
columns, red/green semantic coloring, white dots at segment ends, per-variant labels
at the top, dark background. (Legacy: hand-emitted SVG; the notebook re-renders the
same composition with the standard plotting stack.)

**V5. Failure-tolerant sweep** — a sweep variant where some parameter values are
*expected to fail to build*: each variant’s builder calls are replayed one at a time,
the first failing step is trapped, the longest successfully built prefix is rendered
in place, and the variant label is annotated “(infeasible — n/m built)”. Codifies
that `min_bend_radius` infeasibility throws at build time and shows where the
threshold sits.

**V6. Ensemble scatter colored by the uncertain parameter** — per-particle scatter
plots where marker color encodes the uncertain input (temperature, Viridis): observable
vs T (angle, Stokes components, DLP), and a Poincaré **equatorial projection**
(S1–S2, axes locked to [−1, 1], unit aspect) where the ensemble traces an arc. Hover
shows per-particle values. An explanatory paragraph above the plots states the
physical reading (DLP ≈ 1, angle is the sensitive observable).

**V7. Stacked solver-diagnostic panels** — three aligned panels over s: (top) step
size h on log-y as green accepted dots and dark-red ✕ rejected trials, with scaled
‖K(s)‖ as a shaded background band; (middle) err/tol per trial on log-y with a
dashed threshold at 1, green below / red above; (bottom) the generator’s component
coefficients as labeled line traces. Together they let the eye correlate
solver effort with generator behavior.

**V8. Benchmark bars + table** — grouped bars (first-call vs steady-state) on log-y,
one group per case, plus an HTML table of the raw ms numbers and ratios against the
Float64 baseline. Bars give shape; the table gives the numbers.

**Cross-cutting conventions**: dark backgrounds for 3D scenes and sweeps, light
backgrounds for overlay comparisons; red = subject, green = context everywhere;
equal-aspect for any geometric plot; every figure has a one-line title carrying the
swept parameter; variant labels carry the parameter value (and feasibility
annotation when relevant); prints accompanying a figure are limited to scalar
geometry/solver facts (arc length, writhe, accepted/rejected step counts).
