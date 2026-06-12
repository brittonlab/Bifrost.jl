# Demo Intent

This file records the intent of every visual demo in `bifrost-demos.ipynb`: the
lesson, edge case, or invariant each one codifies, the geometry that produces it
(as Bifrost API calls — precise enough to reconstruct the demo), and what a reader
should look for in the output. Section numbers match the notebook. All plotting
machinery lives in `demo-helper.jl` (`dh_*` functions; see Visual techniques at the
end).

## Class index

| § | Class | Atomic demos |
| --- | --- | --- |
| 1 | Smallest runnable examples | 2 |
| 2 | Path geometry | 5 |
| 3 | Modify (field-level MCM meta) | 12 |
| 4 | Material spin | 1 |
| 5 | Adaptive step-doubling | 1 |
| 6 | JumpBy / JumpTo — 2D sweeps | 5 |
| 7 | Meta × JumpTo interplay — 2D | 3 |
| 8 | JumpBy / JumpTo — 3D scenes | 12 |
| 9 | MCM temperature PTF | 2 |
| 10 | MCM speed benchmarks | 2 |

## §1 — Smallest runnable examples

### 1.0 · Deterministic

```julia
xs = StepIndexCrossSection(SilicaGermaniaGlass(0.036), SilicaGermaniaGlass(0.0),
                           8.2e-6, 125e-6)
sb = SubpathBuilder(); start!(sb)
straight!(sb; length = 0.5); bend!(sb; radius = 0.05, angle = π/2)
straight!(sb; length = 0.5); seal!(sb)
fiber = Fiber(build(sb); cross_section = xs, T_ref_K = 297.15)
J, stats = propagate_fiber(fiber; λ_m = 1550e-9)
```

Lesson: the full layer stack in ~10 lines. The bend's fixed-axis birefringence makes
`J` diagonal `diag(e^{-iφ}, e^{+iφ})`, φ ≈ 0.033 rad — a pure retarder, |det J| = 1.
`stats` has 3 intervals: the propagator never steps across segment boundaries.

### 1.1 · Smallest MCM example

Same geometry; the only change is an uncertain operating temperature:

```julia
fiber = Fiber(build(sb); cross_section = xs, T_ref_K = 297.15 ± 2.0)
J_p, _ = propagate_fiber(fiber; λ_m = 1550e-9)   # entries are Particles
```

Lesson: one uncertain input is all it takes — the whole ensemble propagates in one
call and `J_p`'s entries display as mean ± std. MCM blocks wrap in
`unsafe_comparisons(true)`; ensembles are seeded for notebook reproducibility.

## §2 — Path geometry (geometry layer only, no optics)

Human-in-the-loop visual debugging for `path-geometry*.jl`. Every demo renders the
interactive 3D inspector (§V1).

### 2.1 · simple

```julia
sb = SubpathBuilder(); start!(sb)
straight!(sb; length = 0.10);  bend!(sb; radius = 0.05, angle = π/2)
straight!(sb; length = 0.12);  catenary!(sb; a = 0.03, length = 0.10, axis_angle = 0.0)
bend!(sb; radius = 0.06, angle = π/3);  straight!(sb; length = 0.08)
seal!(sb); path = build(sb)
```

Lesson: all four curve primitives compose with tangent (G1) continuity at every
join; `path_length` and `writhe` (0 — planar path) are scalar geometry queries.

### 2.2 · segment nicknames

Same shape plus `helix!(sb; radius = 0.025, pitch = 0.015, turns = 1.2,
axis_angle = 0.0)`, every segment authored with a descriptive
`meta = [Nickname("lead-in")]` (`90° bend`, `spacer`, `sag`, `spin section`,
`lead-out`). Lesson: `Nickname` flows from authoring through build to presentation;
names live with the segment definition.

### 2.3 · helix axis_angle

```julia
straight!(sb; length = 0.05)
helix!(sb; radius = 0.03, pitch = 0.02, turns = 2.0, axis_angle = aa)  # aa ∈ {0, π/3, 2π/3}
straight!(sb; length = 0.05)
```

Lesson/edge case: `axis_angle` rotates the helix axis about the incoming tangent;
entry stays tangent-continuous for every value, the exit direction changes, and arc
length is invariant (printed: 0.4791 m for all three). Rendered as a §V2 comparison
row plus one full inspector on the `2π/3` variant.

### 2.4 · jumps_min_radius (the paddle)

A `PathBuilt` of 5 Subpaths; vertical 1 m straights joined by 180° U-turn seals:

```julia
sb1: straight!(1.0); jumpto!(point = (1,0,1), incoming_tangent = (0,0,-1), min_bend_radius = 0.4)
sb2: start!(sb2, :inherit); straight!(1.0); jumpto!((2,0,0), (0,0,1),  min_bend_radius = 0.1)
sb3: start!(sb3, :inherit); straight!(1.0); jumpto!((3,0,1), (0,0,-1), min_bend_radius = 0.05)
sb4: start!(sb4, :inherit); straight!(1.0)
     jumpby!(delta = (-1,0,0), tangent = (0,0,-1), min_bend_radius = 0.1)
     jumpto!((2,0,0), (0,0,1))
sb5: start!(sb5, :inherit); straight!(1.0); seal!(sb5)
build([sb1, sb2, sb3, sb4, sb5])
```

Lessons: PathBuilt under a shared global arc length; `:inherit` start states;
interior (`jumpby!`) vs terminal (`jumpto!`) connectors; per-connector
`min_bend_radius`; G2 quintic U-turns.

### 2.5 · pathbuilt (exact handoffs)

```julia
sb1: straight!(0.2); jumpto!(point = (0,0,0.2), incoming_tangent = (0,0,1))
sb2: start!(sb2, :inherit); bend!(radius = 0.05, angle = π/2)
     jumpto!((0.05, 0, 0.25), (1,0,0))      # the bend's analytic exit
sb3: start!(sb3, :inherit); helix!(radius = 0.025, pitch = 0.02, turns = 1.5,
     axis_angle = 0.0); seal!(sb3)
build([sb1, sb2, sb3])
```

Lesson: a `jumpto!` seal pins exact handoff coordinates, the conformity check
validates each boundary, and per-Subpath `Nickname` (on the `SubpathBuilder`
constructor) labels whole subpaths.

## §3 — Modify (field-level MCM meta on segment fields)

Twelve structurally identical experiments. Baselines:

```julia
# inverted-U (targets 1–3):
straight!(sb; length = 1.0)
bend!(sb; radius = 0.5, angle = π, axis_angle = 0.0)
straight!(sb; length = 1.0); seal!(sb)

# helix variant (targets 1–4): same with
helix!(sb; radius = 0.15, pitch = 0.25, turns = 1.5, axis_angle = 0.0)
# inserted between bend and final straight
```

Exactly one segment carries one `MCMadd(:field, δ)` or `MCMmul(:field, f)` with a
plain `Float64` value; `Fiber(sb; ...)` applies it during its single build. The
mutated segment per demo:

| Demo | Mutated segment | Meta | Variants |
| --- | --- | --- | --- |
| 3.1a | first `straight!` | `MCMadd(:length, ·)` | −0.4 |
| 3.1b | `bend!` | `MCMadd(:radius, ·)` | −0.25, +0.50 |
| 3.1c | `bend!` | `MCMadd(:angle, ·)` | −π/2, +π/4, +π |
| 3.2a | first `straight!` | `MCMmul(:length, ·)` | −0.4, +0.5 |
| 3.2b | `bend!` | `MCMmul(:radius, ·)` | 0.5, 2.0 |
| 3.2c | `bend!` | `MCMmul(:angle, ·)` | 0.5, 1.25 |
| 3.3a/d | `helix!` | `MCMadd`/`MCMmul` `(:radius, ·)` | −0.05/+0.10; ×0.5/×2 |
| 3.3b/e | `helix!` | `MCMadd`/`MCMmul` `(:pitch, ·)` | −0.10/+0.20; ×0.5/×2 |
| 3.3c/f | `helix!` | `MCMadd`/`MCMmul` `(:turns, ·)` | ±0.5; ×0.67/×1.5 |

Lessons: which geometry field each knob perturbs; add vs mul semantics (edge case:
`MCMmul(:length, -0.4)` walks the straight backward); with no anchor everything
downstream rides along rigidly — the contrast with §7's anchored paths.

## §4 — Material spin

```julia
sb = SubpathBuilder(); start!(sb; spin_rate = 2π)
straight!(sb; length = 1.0)
helix!(sb; radius = 0.5, pitch = 0.05, turns = 4.0, axis_angle = 0.0)
straight!(sb; length = 1.0); seal!(sb)
```

In the inspector the readout shows τ_spin = 6.2832 rad/m everywhere, ∫τ_spin ds
linear in s, and the red spin arrow rotating in the normal–binormal plane while
T̂/N̂/B̂ follow the geometric frame. Lesson: material spin set at `start!` is
resolved into path-coordinate spin, carried independently of geometric torsion.

## §5 — Adaptive step-doubling

A solver diagnostic on a smooth, noncommuting generator (not a fiber):

    K(s) = α·i·σx·cos(πs) + β·i·σz·sin(2πs),  α = 1.2, β = 0.9, s ∈ [0, 2]

The notebook cell retains the original demo code: `collect_adaptive_steps` (the
recording mirror of the production controller in `path-integral.jl`; the solver is
untouched) with `rtol = 1e-6, atol = 1e-9`, rendered in the original three-panel
diagnostic format (§V7). With these tolerances the run accepts ~370 steps and
rejects ~15. Look for: step size locking to the local error budget and growing where
‖K(s)‖ dips; rejection cascades (red ✕) right after the ‖K‖ minima; err/tol hugging
the threshold from below with no accepted step above 1.

## §6 — JumpBy / JumpTo — 2D sweeps

Scene for every sweep: `straight!(1.0) · connector · straight!(1.0)`, y = 0 plane,
x–z projection (§V3); the connector is the red segment. One degree of freedom per
sweep.

- **6.1 tangent_out** — `jumpby!(sb; delta = (0.4, 0, 0.4), tangent = t)`,
  `t ∈ {(1,0,1)/√2, (0,0,1), (−1,0,1)/√2}` (local frame).
- **6.2 point** — seal `jumpto!(sb1; point = (x, 0, 1.5))`,
  `x ∈ {0, 0.3, 0.6, 1.0}`; continuation `start!(sb2, :inherit)`.
- **6.3 incoming tangent** — `jumpto!(sb1; point = (0.5, 0, 1.5),
  incoming_tangent = t)`, `t ∈ {(1,0,0), (0,0,1), (−1,0,1)/√2}` (global frame).
- **6.4 routing** — three sealed Subpaths chained by `:inherit`:
  `jumpto!((1,0,1), (0,0,−1))`, `jumpto!((2,0,0), (0,0,1))`,
  `jumpto!((3,0,1), (0,0,−1))`, each with `min_bend_radius = 0.1`.
- **6.5 min_bend_radius** — `jumpto!(sb1; point = (1, 0, 1),
  incoming_tangent = (0, 0, −1), min_bend_radius = mbr)`,
  `mbr ∈ {0.10, 0.30, 0.49, 0.51, 0.70}`; infeasible variants throw at `build`
  (trapped via `dh_try_build`, drawn as the surviving lead-in). Lesson: the
  guardrail is build-time and hard. Observed: 0.49 is also infeasible — the quintic
  connector's peak curvature exceeds the circular-arc ideal chord/2 = 0.5, so the
  true threshold for this U-turn lies between 0.30 and 0.49.

## §7 — Meta × JumpTo interplay — 2D overlays

Overlays (§V4): baseline black, modified red, lengths and Δ% in the legend.

- **7.1 drift** — S-curve with fiber-relative jump; perturbation on both bends:

  ```julia
  straight!(0.3)
  bend!(radius = 0.5, angle = π/2, axis_angle = 0.0)   # + MCMmul(:radius, 1.25)
  bend!(radius = 0.5, angle = π/2, axis_angle = π)     # + MCMmul(:radius, 1.25)
  straight!(0.3); jumpby!(delta = (0, 0, 0.8)); straight!(1.0); seal!
  ```

  Lesson: no anchor — the tail shifts and the endpoint separates (~10% longer).
- **7.2 anchor** — same S-curve but sealed with `jumpto!(point = anchor)` where
  `anchor` puts the baseline at 4.000 m (computed from a probe build's `end_point`);
  perturbation `MCMmul(:radius, 1.5)` on both bends. Lesson: the interior swings
  wide, the endpoint stays pinned, the connector absorbs the slack.
- **7.3 thermal seal** — `straight!(0.5)` and its sealing
  `jumpto!(point = (1, 0, 0.5), incoming_tangent = (1, 0, 0))` both carry
  `MCMadd(:T_K, ΔT)` with `ΔT = 0.05 / cte(cladding, T_ref)` (5% expansion).
  Lesson: the length-constrained connector resolve — connector arc length scales by
  τ yet lands on the fixed point; extra length becomes curvature. `:T_K` is
  interpreted by `Fiber`, never by the geometry layer.

## §8 — JumpBy / JumpTo — 3D scenes

§6 restaged in 3D (§V2 rows) plus the G2 content that needs the third dimension.
Base scene as in §6 unless noted; connector red.

- **8.1a axial delta** — `jumpby!(delta = (0, 0, d))`, `d ∈ {0.3, 0.6, 1.0}`.
- **8.1b transverse delta** — `jumpby!(delta = (d, 0, 0.5))`, `d ∈ {0, 0.2, 0.5}`.
- **8.2a tangent_out** — as 6.1.
- **8.2b curvature_out** — `jumpby!(delta = (0.5, 0, 0.5),
  tangent = (1,0,1)/√2, curvature_out = κ)`, `κ ∈ {(0,0,0), (0,2,0), (0,−2,0)}`;
  the ±y variants bow out of plane (G2 exit knob).
- **8.3a/b point sweeps** — `jumpto!(point = (0,0,z))`, `z ∈ {1.3, 1.6, 2.0}`;
  `jumpto!(point = (x,0,1.5))`, `x ∈ {0, 0.3, 0.6}`.
- **8.3c incoming tangent** — as 6.3.
- **8.3d incoming curvature** — `jumpto!(point = (0.5, 0, 1.5),
  incoming_tangent = (0,0,1), incoming_curvature = κ)`,
  `κ ∈ {(0,0,0), (10,0,0), (−10,0,0)}` (G2 landing knob).
- **8.4a jumpby after bend** — `straight!(1.0) · bend!(radius = 0.5, angle = π/2) ·
  jumpby!(delta = d) · straight!(1.0)`, `d ∈ {(0,0,0.5), (0.3,0,0.5), (−0.3,0,0.5)}`.
  Lesson: the delta lives in the rotated local frame.
- **8.4b jumpto after bend** — same scene, `jumpto!(point = p)`,
  `p ∈ {(1,0,1.5), (1.3,0,1.5), (0.5,0,2.0)}`. Lesson: the target stays global.
- **8.5 G2 inheritance** — `bend!(radius = R, angle = π/4) ·
  jumpby!(delta = (0.3, 0, 0.3)) · straight!(1.0)`, `R ∈ {0.3, 0.5, 0.8}`. Lesson:
  with no explicit spec the connector inherits the bend's exit curvature 1/R.
- **8.6 routing 3D** — the 6.4 chain in a 3D scene.

## §9 — MCM temperature PTF

The notebook cells import the legacy demo code largely verbatim (constants,
variable names, flow); only the path spec is authored with the public builder API:

```julia
sb = SubpathBuilder(); start!(sb; spin_rate = s -> sin(2π * s / 100.0))
straight!(sb; length = 5.0)
helix!(sb; radius = 0.025, pitch = 0.05, turns = 10001.892069208387,
       axis_angle = 0.0, meta = AbstractMeta[MCMadd(:T_K, ΔT_K)])
straight!(sb; length = 5.0)
helix!(sb; radius = 0.025, pitch = 0.05, turns = 10000.0, axis_angle = 0.0)
straight!(sb; length = 5.0); seal!(sb)
Fiber(sb; cross_section = MCM_DEMO_XS, T_ref_K = 303.15)
```

Physics: bend birefringence ∝ 1/R²; at R = 2.5 cm the first helix accumulates
|Δβ|·L ≈ 1775·2π rad. The fractional turn count puts 30 °C exactly at mid-fringe
(mod(Γ, 2π) = π), the point of maximum dPTF/dT. Temperature enters twice —
`MCMadd(:T_K, ΔT_K)` geometry meta (thermal expansion via the cladding CTE, baked in
by `Fiber`) and `T_ref_K = T_K_particles` on the final binding (material indices).
Per-particle Jones matrices are sliced into Stokes observables for input state H
(`dh_stokes_ensemble`).

- **9.1 ptf** — `T_C = StaticParticles(50, Normal(30, 5))`; legacy two-panel format
  (§V6): angle vs T | S1/S2/S3/DLP vs T. Look for: angle swinging through the
  mid-fringe, DLP ≈ 1 and S3 ≈ 0 (helix birefringence rotates the linear
  polarization angle without adding ellipticity).
- **9.2 scatter** — `Particles(500, Normal(30, 5))` (legacy used 2000; 500 keeps
  runtime in minutes); two-panel: angle vs T | S1–S2 Poincaré equatorial projection
  with unit circle. Look for: the temperature-parameterized arc on the equator.

## §10 — MCM speed benchmarks

The notebook cells import the legacy benchmark structure verbatim: the same four
cases (`Float64`, `Particles(2000)`, `StaticParticles(50)`, `StaticParticles(75)`),
the same two scenarios (propagate-only on a pre-built fiber; build + propagate
including `:T_K` interpretation), first-call (JIT + run, `time_ns`) vs steady-state
(`@belapsed` minimum, 3 samples), and the §V8 presentation (grouped log-scale bars +
table). Numbers are machine-dependent; nothing is asserted.

Notebook scaling: the benchmark fiber keeps the 5-segment structure but scales the
helices to 100 turns (the legacy ~10000-turn, mid-fringe-tuned fiber exists for
sensitivity, which is irrelevant to a speed comparison) so the repeated timed runs —
in particular the `Particles(2000)` case — stay tractable; a note in the cell
records the reduction.

## Visual techniques (implemented in `demo-helper.jl`)

- **V1 path inspector** (`dh_path_inspector`) — interactive 3D scene: centerline
  with arc-length-graded markers, open circles at segment joins, green start / red
  end dots, `Nickname` labels floated above segment midpoints; a slider scrubs arc
  length carrying a translucent normal–binormal plane, the T̂/N̂/B̂ triad
  (orange/blue/green), and a red in-plane spin-phase arrow; a live readout shows
  s; x, y, z; κ; τ_geom; τ_spin; ∫τ_spin ds. Axes are fixed to the bounding box of
  all content (path + cursor extents) so they never re-range while scrubbing.
- **V2 variant rows 3D** (`dh_variant_row`) — N variants of a path side by side
  with constant x-offsets; subject segments red, context green, terminal connectors
  optionally faint gray; per-variant start/end dots and labels.
- **V3 variant rows 2D** (`dh_variant_row_2d`) — the same for planar scenes in the
  x–z projection, equal-aspect axes.
- **V4 overlays** (`dh_overlay_compare`) — baseline (black) vs modified (red) in
  x–z, with segment-boundary ticks and a legend carrying lengths and Δ%.
- **V5 try-build** (`dh_try_build`) — traps expected build-time errors so
  infeasible variants render as labeled placeholders instead of aborting a cell.
- **V6 MCM temperature rows** (`dh_temperature_ptf_row`,
  `dh_temperature_scatter_row`) — the legacy demo3mcm presentation: two side-by-side
  panels; angle vs T with a shared Viridis temperature colorbar, and either
  S1/S2/S3/DLP vs T (Reds/Greens/Blues/Oranges per series, DLP as diamonds) or the
  S1–S2 equatorial projection with unit circle at locked [−1, 1] unit aspect.
  `dh_jones_to_stokes` / `dh_stokes_ensemble` convert per-particle Jones matrices to
  Stokes observables.
- **V7 adaptive panels** (`dh_adaptive_panels`) — the legacy three-panel
  diagnostic: step size h(s) with accepted dots / rejected ✕ over a shaded ‖K(s)‖
  profile; err/tol per trial on a log axis with the acceptance threshold dashed;
  the generator's component coefficients.
- **V8 benchmark presentation** (`dh_benchmark_chart`, `dh_benchmark_table`) — the
  legacy format: grouped first-call vs steady-state bars on a log-time axis, and a
  compact Markdown table (label, N particles, first, steady, steady-per-particle).

**Cross-cutting conventions**: dark backgrounds throughout; red = subject,
green = context everywhere; equal/fixed aspect for any geometric plot; every figure
title carries the swept parameter; variant labels carry the parameter value (and
feasibility annotation when relevant); prints accompanying a figure are limited to
scalar geometry/solver facts (arc length, writhe, step counts); MCM ensembles are
seeded (`Random.seed!`) so the notebook is deterministic.
