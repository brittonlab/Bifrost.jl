# The transverse frame and the birefringence gauge

This page derives, from first principles, the frame convention in which BIFROST
expresses polarization: the **parallel-transported (Bishop) frame**. It fixes the
meaning of every angle that enters the generators, explains why geometric torsion
appears nowhere in the optics, and states the anchor and Subpath-continuity
conventions. It is the companion to the implementation in
`src/geometry/path-geometry.jl` (frame transport) and `src/fiber/fiber-path.jl`
(generator assembly), and to the visual demo
`test/human/demo-frame-birefringence.jl`.

## 1. The propagation equation and what "gauge" means here

BIFROST integrates the Jones equation along the fiber arc length $s$:

```math
\frac{dJ(s)}{ds} = K(s)\,J(s), \qquad J(0) = \mathbb{1},
```

with $K$ a traceless anti-Hermitian-times-$i$ generator (lossless propagation;
see [Deriving the generators](generators.md)). The two components of the Jones
vector are the field amplitudes on **some** pair of transverse unit vectors
$(e_1(s), e_2(s))$ with $e_2 = \hat{T}\times e_1$. That pair is a *gauge choice*:
nothing physical depends on it, but every axis angle inside $K$ — the bend axis,
the ellipse axis — is measured *in* it, so the equations are only as trustworthy
as the frame convention is consistent. The solver (`src/path-integral.jl`) never
sees the frame; the gauge lives entirely in how `fiber-path.jl` builds $K(s)$
from geometry queries.

## 2. Free propagation parallel-transports the field

For a transparent, locally isotropic medium guiding a transverse wave along a
gently curved path (radii large compared to the wavelength and the cladding,
which is the regime BIFROST is restricted to), Maxwell's equations force the
electric field to stay transverse while exchanging direction with the tangent
only at the minimal rate required by $E \perp \hat{T}$:

```math
\frac{d\vec{E}}{ds} = -\left(\vec{k}\cdot\vec{E}\right)\hat{T},
\qquad \vec{k} \equiv \frac{d\hat{T}}{ds},
```

where $\vec{k}(s)$ is the **curvature vector** (magnitude $\kappa$, direction
toward the local center of curvature). The field never rotates *about* the
tangent. This is Rytov's law, observed in fiber as the Berry-phase rotation of
polarization along helical paths (Tomita and Chiao,
doi:10.1103/PhysRevLett.57.937; Berry, doi:10.1038/326277a0). The glass itself
obeys the same rule: an unspun, untwisted cross-section is carried along the
path without rotating about the tangent, so its frozen-in anisotropy axes are
parallel-transported too (Ulrich and Simon, doi:10.1364/AO.18.002241).

## 3. The Bishop (relatively-parallel) frame

Define the transverse pair by the same transport law (Bishop,
doi:10.2307/2319846):

```math
\frac{de_1}{ds} = -\left(\vec{k}\cdot e_1\right)\hat{T},
\qquad e_2 = \hat{T}\times e_1 .
```

Properties that matter here:

- **Existence and smoothness everywhere.** The transport needs only $\hat T(s)$;
  it is defined on straights, through inflections, and across joints where the
  curvature direction jumps. (The Frenet–Serret normal $\hat N = \vec k/\kappa$
  needs $\kappa \neq 0$ and flips at inflections.)
- **Zero twist:** $\langle de_1/ds,\, e_2\rangle = 0$ by construction — the
  defining property, enforced as a T-PHYSICS test in
  `test/test_bishop_frame.jl`.
- **Gauge freedom is one constant angle.** If $(e_1, e_2)$ is relatively
  parallel, so is the pair rotated by any *constant* angle about $\hat T$. The
  frame is unique once $e_1(0)$ is chosen (Section 7).
- **Relation to Frenet–Serret.** Where the FS frame exists,
  $e_1 = \hat N\cos\theta - \hat B\sin\theta$ with
  $\theta(s) = \theta(0) + \int_0^s \tau_{\mathrm{geom}}\,ds'$: the FS pair
  rotates about the tangent at the geometric-torsion rate relative to parallel
  transport. This identity is the closed-form helix transport in
  `_parallel_transport_local(::HelixSegment, …)`.

**The "M1 vs M2" question dissolves.** A common objection to Bishop frames is
that one cannot tell how the osculating plane lies relative to the two frame
legs. That question never needs an answer: the frame is *gauge*, not physics.
Physical directions — the curvature vector, the ellipse major axis — are
*projected onto* the pair, and angles like
$\theta_b = \operatorname{atan2}(\vec k\cdot e_2,\ \vec k\cdot e_1)$ are
continuous functions of the geometry wherever the physics is continuous, and
jump exactly where the physics jumps (a corner between bend planes), always on
an existing breakpoint.

## 4. Why the generator contains no geometric-torsion term

Expand the field on the transported pair,
$\vec E = a_1 e_1 + a_2 e_2$. Substituting the transport law for both $\vec E$
(Section 2) and $(e_1, e_2)$ (Section 3) cancels the tangential exchange terms
identically, leaving

```math
\frac{d}{ds}\begin{pmatrix}a_1\\ a_2\end{pmatrix} = K_{\mathrm{aniso}}(s)
\begin{pmatrix}a_1\\ a_2\end{pmatrix},
```

where $K_{\mathrm{aniso}}$ contains **only material anisotropy**. Free
propagation — straight, bent, or helical — contributes nothing. Had we expanded
on the FS pair instead, the relative rotation of Section 3 would inject a
connection term $\tau_{\mathrm{geom}}\,(i\sigma_2)$ that must then be carried,
spliced, and (as issues #88/#62/#24 showed) mis-carried across joints,
inflections, and near-straight regions. In the transported gauge
$\tau_{\mathrm{geom}}$ exits the optics entirely; it survives only as the shape
diagnostic `geometric_torsion`/`total_torsion`, with no optical role.

## 5. The birefringence sources in this gauge

With the gauge fixed, each mechanism enters $K$ as a linear retarder
$\tfrac{i}{2}\Delta\beta\,(\cos 2\varphi\,\sigma_3 + \sin 2\varphi\,\sigma_1)$
or a circular rotator $\tfrac{1}{2}\Delta\beta_c\,(-i\sigma_2)$, with axis
$\varphi$ measured from $e_1$:

| Source | Axis $\varphi(s)$ | Magnitude | Code |
| --- | --- | --- | --- |
| Bend (curvature stress) | $\theta_b = \operatorname{atan2}(\vec k\cdot e_2, \vec k\cdot e_1)$ | $\Delta\beta_b \propto \kappa^2$ | `bend_generator_K` |
| Axial tension on a bend | $\theta_b$ (same eigen-axis) | $\propto F\kappa$ | `tension_generator_K` |
| Core ellipticity + asymmetric thermal stress | $\theta_{\mathrm{int}} = \theta_{\mathrm{int}}(0) + \phi_{\mathrm{spin}}(s) + \phi_{\mathrm{twist}}(s)$ | from cross-section | `ellipticity_generator_K` |
| Mechanical twist | — (circular) | $\Delta\beta_c = g\,\tau_m$ | `twist_generator_K` |

- $\vec k\cdot e_{1,2}$ are `bend_components` (`fiber-path.jl`): the projection
  of `curvature_vector(path, s)` onto `bishop_e1`/`bishop_e2`. The double-angle
  form $(\cos 2\theta_b, \sin 2\theta_b) = ((k_x^2-k_y^2)/k^2,\ 2k_xk_y/k^2)$
  is normalization- and branch-free (`_bend_axis_c2s2`).
- $\phi_{\mathrm{spin}} = \int_0^s \xi\,ds'$ (`spin_phase`) is the frozen-in
  rotation of the glass imparted during draw: it carries the intrinsic axes but
  creates no stress and no circular birefringence.
- $\phi_{\mathrm{twist}} = \int_0^s \tau_m\,ds'$ (`twist_phase`) is elastic
  twist: it both co-rotates the intrinsic axes *and* produces photoelastic
  circular birefringence $g\,\tau_m$ with $g \approx 0.14$–$0.16$ in silica
  (Ulrich and Simon, doi:10.1364/AO.18.002241).
- Bend-stress magnitudes follow Ulrich, Rashleigh, and Eickhoff
  (doi:10.1364/OL.5.000273); the spun-bent formulation parallels
  Przhiyalkovsky et al. (doi:10.1109/JLT.2020.3017795).

Three rotations that must never be conflated: **geometric torsion** (gauge —
absent here), **spin** (material rotation of the anisotropy axes, no stress),
**twist** (material rotation *plus* circular birefringence). Only the last two
are physical inputs; both are integrals of authored rates and are continuous
across Subpath boundaries by construction.

## 6. DGD and gauge invariance

The sensitivity system propagates $G = \partial_\omega J$ alongside $J$, and
the differential group delay is the eigenvalue spread of
$H = -i\,J^{-1}G$ (`output_dgd`, `output_dgd_2x2`). A change of anchor or any
fixed rotation of the gauge conjugates $J \mapsto R\,J\,R^{\mathsf T}$ and $G$
likewise, which conjugates $H$ and leaves its spectrum — hence the DGD and all
polarization observables — unchanged. This is why the pre-refactor helix
reference value is reproduced exactly (`test/test_bishop_frame.jl`,
T-SIM-REGRESSION) even though $J$ itself is reported in a different gauge than
the old code's.

## 7. The anchor: $e_1(0)$

Parallel transport determines the frame up to one constant: $e_1(0)$. BIFROST
uses a **static lab-frame rule** (`_initial_frame_from_tangent`): take the
world axis least aligned with the launch tangent and Gram–Schmidt it,

```math
e_1(0) = \frac{\hat a - (\hat a\cdot\hat T)\hat T}{\lVert \cdots \rVert},
\qquad \hat a = \operatorname*{arg\,min}_{\hat x,\hat y,\hat z}\ |\hat T\cdot \hat a|.
```

For the common $\hat T(0) = \hat z$ launch this is simply $e_1(0) = \hat x$.
The rule is deliberately blind to curvature: it depends only on the launch
direction, is reproducible by hand, never moves when downstream geometry is
edited, and has no special cases. Alternatives considered — "first curvature
direction" (undefined on a leading straight, unstable under small geometry
edits) and a user-supplied vector (API surface for a pure gauge choice) — buy
nothing physical, because observables are anchor-independent (Section 6). To
report $J$ on specific laboratory axes, conjugate by the constant rotations
between $(e_1, e_2)$ and those axes at the two fiber ends — both available from
`bishop_e1`/`bishop_e2` at $s = 0$ and $s = L$.

## 8. Subpath boundaries

Each `SubpathBuilt` is built standalone with its own lab anchor. Stacking
subpaths into a `PathBuilt` must not let the gauge jump mid-fiber, so
`build(::Vector{SubpathBuilt})` runs `_resolve_bishop_gauge`: for each later
Subpath it stores the one constant angle (`_bishop_gauge_at_s0`) that rotates
its anchored frame onto the predecessor's transported frame at the boundary.
Because a constant transverse rotation of a relatively-parallel field is itself
relatively parallel, the correction is **exact** — no re-transport, no
approximation. Consequently the optical gauge is continuous across every
interior boundary (hand-loaded or `:inherit`), and a path built as one Subpath
or as several produces identical generators and identical $J$
(`test/test_bishop_frame.jl`, issue #89).

Authoring caveat, distinct from the gauge: segment parameters that reference
the *construction* frame (`axis_angle`, `jumpby!` deltas) are interpreted in a
frame re-derived from the boundary tangent, not continued from the predecessor
— so frame-dependent authoring after a split can change the 3D shape itself.
That is a property of the authoring DSL owned by the subpath-concatenation work
(issues #51, #32), not of the optical gauge.

## 9. Why the Frenet gauge failed (the defect inventory)

The pre-refactor code oriented the bend axis by
$\varphi = \int_0^s \tau_{\mathrm{geom}}\,ds'$ — correct only while the FS
normal evolves smoothly. Documented failures, each reproduced visually in
`test/human/demo-frame-birefringence.jl`:

1. **Discrete curvature-direction jumps dropped** (#88): `axis_angle` never
   reached the generators; perpendicular-plane corners and S-bends produced the
   same Jones axes as their planar counterparts.
2. **Connector torsion spikes**: FS torsion's $\sim\kappa^2$ denominator made
   a near-straight quintic connector inject large spurious axis phase.
3. **Internal contradiction**: `total_frame_rotation` skipped connector
   torsion while `torsion_phase` integrated it by quadrature.
4. **Subpath gauge scrambling** (#89): each Subpath re-anchored its frame, so
   the relative bend-vs-intrinsic angle jumped at boundaries.
5. **Helix display jumps** (#24): the τ-inclusive "frame rate" diagnostic
   jumped at helix entry/exit.
6. The one configuration the hybrid got right — helix followed by
   `bend!(…; axis_angle = 0)`, because chaining FS end frames absorbs exactly
   $\int\tau$ — is preserved as a regression case (DGD equality with the old
   code, Section 6).

The S-shape conformity rejection (#62) is the same Frenet-flip pathology
appearing in the *authoring* layer's endpoint checks; its fix is to compare
curvature **vectors** (which are continuous through an S-joint) rather than
unit normals, and belongs to the conformity check, not the gauge.

## 10. Implementation map

| Concept | Where |
| --- | --- |
| Per-segment closed-form transport | `_parallel_transport_local` (straight/bend/catenary/helix), `path-geometry.jl` |
| Connector discrete transport (double reflection, doi:10.1145/1330511.1330513) | `e1_table` + `_parallel_transport_local(::QuinticConnector, …)`, `path-geometry-connector.jl` |
| Transported frame queries | `bishop_e1`, `bishop_e2`, `bishop_frame` |
| Curvature vector | `curvature_vector`, `_curvature_vector_local` |
| Anchor | `_initial_frame_from_tangent` |
| Subpath gauge continuity | `_resolve_bishop_gauge`, `_bishop_gauge_at_s0` |
| Bend-axis projection | `bend_components`, `_bend_axis_c2s2`, `fiber-path.jl` |
| Material phases | `spin_phase`, `twist_phase` |
| Shape diagnostics (no optical role) | `geometric_torsion`, `total_torsion`, `writhe` |
