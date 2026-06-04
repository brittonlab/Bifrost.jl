# path-geometry.jl — concepts and integral quantities

## The sliding frame

Every segment is expressed in a **local coordinate frame** whose z-axis is the
incoming tangent direction.  When segments are chained, the exit frame of one
segment becomes the entry frame of the next, so tangent continuity is guaranteed
by construction.  There is no need to specify absolute orientations; only the
relative geometry of each segment matters.

The Frenet–Serret frame (tangent T̂, principal normal N̂, binormal B̂) is
carried along the centerline.  The Frenet–Serret equations relate the frame's
rate of change to two scalar fields:

    dT̂/ds =  κ N̂
    dN̂/ds = −κ T̂ + τ B̂
    dB̂/ds =        −τ N̂

- **Curvature** κ(s): the magnitude of dT̂/ds, in 1/m — how fast the tangent
  direction turns.  A straight segment has κ = 0.  A circular bend of radius R
  has κ = 1/R.  This is the quantity integrated by `total_turning_angle`.
- **Geometric torsion** τ_geom(s): the rate at which the osculating plane (and
  thus N̂ and B̂) rotates about the tangent T̂, in rad/m.  A planar curve — any
  straight segment or circular bend — stays in one plane and has τ_geom = 0.  A
  circular helix is the canonical out-of-plane curve and has constant nonzero
  τ_geom (see the example below).
- **Spinning** Ω(s): rate at which the fiber cross-section rotates relative
  to the propagation frame, in rad/m.  This is *manufacturing spin* — the fiber
  rotated while molten during the draw, freezing a rotating axis into relaxed
  glass.  It is specified once per Subpath via the `start!(; spin_rate=…)`
  keyword (`nothing`, a constant rate, a callable of Subpath-local arc length,
  or `:inherit`), **not** as per-segment metadata.  Spinning rotates the
  intrinsic linear birefringence axes but introduces no circular birefringence.
- **Mechanical twist** τ_m(s): rate at which the *solid* fiber is twisted about
  its axis after manufacture (rotating the ends / spooling), in rad/m.  It is a
  **per-segment geometric parameter** — the optional `twist` keyword on every
  building call (`straight!`, `bend!`, `helix!`, `catenary!`, `jumpby!`,
  `jumpto!`, `seal!`), accepting `nothing`, a constant rate, or a callable of
  segment-local arc length — and joins the `AbstractPathSegment` interface as
  `twist_rate(seg, s)`.  Mechanical twist does not change the centerline shape
  (so contributes no τ_geom); the fiber layer converts it into circular
  birefringence and co-rotates the intrinsic linear axes with the cross section.

The propagation frame used for polarization is the **parallel-transport
(Bishop) frame**, not the Frenet frame: it does not spin with the curve's
torsion, giving a well-defined transverse basis even where κ = 0.  The
geometric frame rotation rate is the sum of torsion and spin:

    dψ/ds = τ_geom(s) + Ω(s)

Per-arc-length accumulated phases are available for orienting birefringence
axes (consumed by the fiber generators): `spin_phase(path, s) = ∫₀ˢ Ω`,
`twist_phase(path, s) = ∫₀ˢ τ_m`, and `torsion_phase(path, s) = ∫₀ˢ τ_geom`.
These propagate MCM `Particles` (they are not nominalized), unlike the
breakpoint-oriented `total_*` integrals below.

---

## Integral quantities compared

Four scalar integrals are available.  They measure related but distinct things:

### `total_turning_angle(path)`

    ∫ κ(s) ds

The integrated curvature κ(s) (defined above), in radians.  Measures how much
the tangent direction T̂ has turned — the "winding" of the path in 3D.
Independent of torsion and spinning.  For a closed planar loop this equals 2π.

### `total_torsion(path)`

    ∫ τ_geom(s) ds

The integrated geometric torsion, in radians.  Zero for paths made entirely of
straight segments and circular bends (both have τ_geom = 0).  Nonzero for
helices and other out-of-plane curves.  This is a property of the centerline
shape alone; it does not depend on the fiber material or applied torques.

### `total_spinning(path; s_start, s_end)`

    ∫ Ω(s) ds

The integrated applied material spinning, in radians.  Counts only the Subpath's
`spin_rate` (set at `start!`).  Knows nothing about the geometry of the
centerline — a helix with no spin contributes zero here, even though the fiber
reference frame does rotate as it traverses the helix.  Mechanical twist τ_m is
a separate, per-segment quantity and is **not** included in `total_spinning`.

### `total_frame_rotation(path; s_start, s_end)`

    ∫ [τ_geom(s) + Ω(s)] ds

The total rotation of the polarization reference frame, in radians.  This is
the physically meaningful quantity for polarization propagation: it captures
both the frame rotation due to the shape of the path (geometric torsion) and
any additional rotation applied to the fiber material (material spinning).  Use
this when you want the net polarization-axis rotation over a segment or the
whole path.

---

## Quick reference

| Function | Integrand | Shape | Spinning |
|---|---|---|---|
| `total_turning_angle` | κ(s) | ✓ | — |
| `total_torsion` | τ_geom(s) | ✓ | — |
| `total_spinning` | Ω(s) | — | ✓ |
| `total_frame_rotation` | τ_geom(s) + Ω(s) | ✓ | ✓ |

---

## Example: helix vs. circular loops

A helix with radius R, pitch-parameter h = pitch/(2π), and arc length L has:

    τ_geom = h / (R² + h²)   [constant along the helix]
    total_torsion         = τ_geom · L
    total_frame_rotation  = τ_geom · L + total_spinning

A sequence of circular `bend!` segments forming a coil has τ_geom = 0
everywhere, so `total_torsion = 0` and `total_frame_rotation` reduces to
`total_spinning` alone — even though the fiber winds around in 3D.  The
out-of-plane geometry of a coil does not, by itself, rotate the Frenet frame;
only a true helical path (nonzero torsion) does.
