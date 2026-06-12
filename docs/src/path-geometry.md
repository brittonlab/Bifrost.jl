# path-geometry.jl — concepts and integral quantities

## The sliding (construction) frame

Every segment is expressed in a **local coordinate frame** whose z-axis is the
incoming tangent direction.  When segments are chained, the exit frame of one
segment becomes the entry frame of the next, so tangent continuity is guaranteed
by construction.  There is no need to specify absolute orientations; only the
relative geometry of each segment matters.  This chained frame is internal to
`build`: it interprets each segment's `axis_angle` and the connector boundary
data, and is not exposed by the query API.

## The transported (Bishop) frame

The transverse frame returned by `bishop_e1`/`bishop_e2`/`frame` is the
**parallel-transported (Bishop) pair** (e1, e2 = T̂ × e1):

    de1/ds = −(k⃗ · e1) T̂,        k⃗ ≡ dT̂/ds

It has zero twist about the tangent, is continuous along the whole path —
through straights, inflections, and joints where the curvature direction jumps
— and is the gauge in which the fiber layer expresses every birefringence
axis.  It is anchored at s = 0 by a static lab-frame rule (Gram–Schmidt of the
world axis least aligned with the launch tangent) and is continuous across
Subpath boundaries in a `PathBuilt`.  See
[The transverse frame and the birefringence gauge](frame-and-gauge.md) for the
physics.

Scalar fields along the path:

- **Curvature** κ(s): the magnitude of dT̂/ds, in 1/m — how fast the tangent
  direction turns.  A straight segment has κ = 0.  A circular bend of radius R
  has κ = 1/R.  This is the quantity integrated by `total_turning_angle`.  The
  full **curvature vector** k⃗(s) (global coordinates, zero on straights) is
  exposed by `curvature_vector(path, s)`; the bend birefringence axis is its
  projection onto (e1, e2).
- **Geometric torsion** τ_geom(s): the rate at which the osculating plane
  rotates about the tangent, in rad/m.  A planar curve has τ_geom = 0; a helix
  has constant nonzero τ_geom.  This is a **shape diagnostic only** — in the
  transported gauge it never enters polarization propagation.
- **Spin** Ω(s): rate at which the fiber cross-section's frozen-in anisotropy
  axes rotate (imparted during draw), in rad/m, specified per Subpath via the
  `start!(; spin_rate=…)` keyword.  The rate may be constant or a callable
  function of run-local arc length.  Spin is a material property; it does not
  rotate the transported frame.

---

## Integral quantities compared

Three scalar integrals are available.  They measure related but distinct
things:

### `total_turning_angle(path)`

    ∫ κ(s) ds

The integrated curvature κ(s) (defined above), in radians.  Measures how much
the tangent direction T̂ has turned — the "winding" of the path in 3D.
Independent of torsion and spin.  For a closed planar loop this equals 2π.

### `total_torsion(path)`

    ∫ τ_geom(s) ds

The integrated geometric torsion, in radians.  Zero for paths made entirely of
straight segments and circular bends (both have τ_geom = 0).  Nonzero for
helices and other out-of-plane curves.  This is a property of the centerline
shape alone; it does not depend on the fiber material or applied torques.

### `total_spin(path; s_start, s_end)`

    ∫ Ω(s) ds

The integrated applied material spin, in radians.  Only counts contributions
from the Subpath's `spin_rate` (set at `start!`).  Knows nothing about the
geometry of the centerline — a helix with no `spin_rate` contributes
zero here.  This is the rotation of the intrinsic birefringence axes in the
transported frame (together with `twist_phase` for elastic twist).

---

## Quick reference

| Function | Integrand | Shape | Spin |
|---|---|---|---|
| `total_turning_angle` | κ(s) | ✓ | — |
| `total_torsion` | τ_geom(s) | ✓ | — |
| `total_spin` | Ω(s) | — | ✓ |

---

## Example: helix vs. circular loops

A helix with radius R, pitch-parameter h = pitch/(2π), and arc length L has:

    τ_geom = h / (R² + h²)   [constant along the helix]
    total_torsion = τ_geom · L

A sequence of circular `bend!` segments forming a coil has τ_geom = 0
everywhere, so `total_torsion = 0` — even though the fiber winds around in 3D.
In both cases the transported frame (and with it the polarization gauge) never
rotates about the tangent; what differs between the helix and the coil is how
the **curvature direction** moves relative to that frame: continuously at rate
τ_geom inside a helix, and by discrete jumps at the joints of a coil.  Both
enter the optics through the projection of `curvature_vector` onto the frame.
