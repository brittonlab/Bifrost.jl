# path-geometry.jl — concepts and integral quantities

## The sliding frame

Every segment is expressed in a **local coordinate frame** whose z-axis is the
incoming tangent direction.  When segments are chained, the exit frame of one
segment becomes the entry frame of the next, so tangent continuity is guaranteed
by construction.  There is no need to specify absolute orientations; only the
relative geometry of each segment matters.

The Frenet–Serret frame (tangent T̂, principal normal N̂, binormal B̂) is
carried along the centerline.  Two scalar fields characterize how fast the frame
rotates:

- **Geometric torsion** τ_geom(s): rate at which B̂ rotates about T̂, in rad/m.
  A pure circular arc has τ_geom = 0.  A helix has constant nonzero τ_geom.
- **Material twist** Ω(s): rate at which the fiber cross-section rotates relative
  to the Frenet frame, in rad/m.  This arises from applied torque or spinning
  during fiber lay-down, and is specified via `twist!` overlays.

The total rotation rate of the polarization reference frame is their sum:

    dψ/ds = τ_geom(s) + Ω(s)

---

## Integral quantities compared

Four scalar integrals are available.  They measure related but distinct things:

### `total_turning_angle(path)`

    ∫ κ(s) ds

The total curvature, in radians.  Measures how much the tangent direction T̂
has turned — the "winding" of the path in 3D.  Independent of torsion and
twist.  For a closed planar loop this equals 2π.

### `total_torsion(path)`

    ∫ τ_geom(s) ds

The integrated geometric torsion, in radians.  Zero for paths made entirely of
straight segments and circular bends (both have τ_geom = 0).  Nonzero for
helices and other out-of-plane curves.  This is a property of the centerline
shape alone; it does not depend on the fiber material or applied torques.

### `total_material_twist(path; s_start, s_end)`

    ∫ Ω(s) ds

The integrated applied material twist, in radians.  Only counts contributions
from explicit `TwistOverlay`s added via `twist!`.  Knows nothing about the
geometry of the centerline — a helix with no `twist!` overlay contributes zero
here, even though the fiber reference frame does rotate as it traverses the
helix.

### `total_frame_rotation(path; s_start, s_end)`

    ∫ [τ_geom(s) + Ω(s)] ds

The total rotation of the polarization reference frame, in radians.  This is
the physically meaningful quantity for polarization propagation: it captures
both the frame rotation due to the shape of the path (geometric torsion) and
any additional rotation applied to the fiber material (material twist).  Use
this when you want the net polarization-axis rotation over a segment or the
whole path.

---

## Quick reference

| Function | Integrand | Shape | Material twist |
|---|---|---|---|
| `total_turning_angle` | κ(s) | ✓ | — |
| `total_torsion` | τ_geom(s) | ✓ | — |
| `total_material_twist` | Ω(s) | — | ✓ |
| `total_frame_rotation` | τ_geom(s) + Ω(s) | ✓ | ✓ |

---

## Example: helix vs. circular loops

A helix with radius R, pitch-parameter h = pitch/(2π), and arc length L has:

    τ_geom = h / (R² + h²)   [constant along the helix]
    total_torsion         = τ_geom · L
    total_frame_rotation  = τ_geom · L + total_material_twist

A sequence of circular `bend!` segments forming a coil has τ_geom = 0
everywhere, so `total_torsion = 0` and `total_frame_rotation` reduces to
`total_material_twist` alone — even though the fiber winds around in 3D.  The
out-of-plane geometry of a coil does not, by itself, rotate the Frenet frame;
only a true helical path (nonzero torsion) does.
