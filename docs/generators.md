# Deriving the Generators

BIFROST tries to solve the differential equation

$$ \frac{\partial J(z)}{\partial z} = K(z) J(z) $$

where $J(z)$ is the Jones matrix for our system (which also depends on e.g. $\omega$, the angular frequency of the light), $z$ is the position coordinate along the fiber, and $K(z)$ is a "generator." The generator contains all of the fundamental physics in the problem, so this document shows how to derive it.

Our formulation is very similar to that of Przhiyalkovsky et al., "Polarization Dynamics of Light Propagating in Bent Spun Birefringent Fiber," *Journal of Lightwave Technology* 24(38), 6879-6885 (2020), DOI: 10.1109/jlt.2020.3017795.

## Formulating the Equation

Suppose our fiber is laid along the $z$-axis. Propagation of monochromatic light in the fiber is dictated by the wave equation

$$ \frac{\partial^2\vec{E}}{\partial z^2} = -k_0^2 \epsilon \vec{E}. $$

Here $k_0 = \omega/c$ and $\epsilon$ is the permittivity tensor. This is what contains all of the interesting physics; we'll come back to it in a moment.

We approximate the waves as transverse, so $\vec{E} = (E_x, E_y)$ (i.e. only two components). There is a common phase gained by both polarization modes of the form $e^{-i k_0 \bar{n} z}$ where $\bar{n}$ is an average refractive index. Let's transform this away by defining $E_{x,y} = \mathcal{E}_{x,y} e^{-ik_0 \bar{n}z}$. If the anisotropy is low ($\epsilon_i, \epsilon_c \ll \bar{\epsilon}$), then this acts like a slowly-varying envelope approximation and, upon substitution into our wave equation, we can ignore the second-order derivatives of $\mathcal{E}_{x,y}$, which turns our wave equation into a Schrödinger-like equation:

$$ 2i\bar{n} \frac{\partial \mathcal{E}}{\partial z} = k_0(\epsilon - \bar{n}^2)\mathcal{E}. $$

So the "generator" we're interested in is

$$ K(z) = -\frac{i}{2} \frac{k_0}{\bar{n}} (\epsilon - \bar{n}^2). $$

## The Permittivity Tensor

The permittivity tensor has contributions from every birefringence mechanism available. It explicitly describes how the two modes gain phase, so we can use our physical intuition to write it down. We must pick a basis for doing so; let's choose the $N,B$ basis where $N$ and $B$ are the normal and binormal of the Frenet-Serret frame. Then, as a starting point, we can write

$$ \epsilon = \bar{\epsilon} + i\epsilon_c \left( \begin{array}{cc} 0 & -1 \\ 1 & 0 \end{array} \right) + \epsilon_i \left( \begin{array}{cc} \cos 2\xi z & \sin 2\xi z \\ \sin 2\xi z & -\cos 2\xi z \end{array} \right). $$

Here $\bar{\epsilon} = \bar{n}^2$ is the average permittivity, $\epsilon_c$ is the circular birefringence, and $\epsilon_i$ is the intrinsic linear birefringence, which is rotating with the spinning of the fiber at spin rate $\xi$. Writing $\xi = 2\pi/L_s$ with a single spin pitch $L_s$ is the constant-rate special case used here for clarity. In general the spin rate is an arbitrary function of arc length, $\xi = \xi(z)$, and the implementation treats it that way: a Subpath's `spin_rate` (set at `start!(; spin_rate=…)`) in [`path-geometry.jl`](../src/geometry/path-geometry.jl) is a `Union{Nothing, Float64, Function}`, and the run integrator (`_integrate_rate`) takes the analytic branch for a constant rate and adaptive Gauss–Kronrod quadrature for a function rate. Replace $2\xi z$ below with the accumulated spin phase $2\,\phi_\xi(z) = 2\int_0^z \xi(z')\,dz'$ (queryable as `spin_phase(path, z)`) for the function-valued case.

This is only a starting point because the choice of the $N,B$ frame doesn't really mean anything. Without the curvature of bending, the direction of the normal vector is arbitrary. In this case, the normal vector direction is evidently set by the fact that the first component of the Jones vector here is along the slow axis and the second component is along the fast axis. Thus the normal vector is set along the slow axis at $z=0$ (and doesn't rotate with the spinning). That's a loose constraint because we sort of don't care about these overall rotations for our problems where we'll be sending in light of all polarizations.

## Three rotations: geometric torsion, spin, and mechanical twist

Once the centerline bends and twists, three *physically distinct* rotations enter, and the implementation (issue #8) keeps them separate. The propagation frame is the **parallel-transport (Bishop) frame** of the centerline — not the Frenet frame — because the Bishop frame does not spin with the curve's torsion and so gives a well-defined transverse basis even where the curvature vanishes.

| Source | Origin | Effect on $\epsilon$ |
| --- | --- | --- |
| **Geometric torsion** $\tau_g(z)$ | The 3D *shape* of the centerline (Frenet–Serret torsion). Zero for any planar curve; constant $h/(R^2+h^2)$ for a helix. | Rotates the **curvature-direction (bend) linear axis** relative to the Bishop frame at rate $\tau_g$. No circular birefringence. |
| **Manufacturing spin** $\xi(z)$ | The fiber rotated about its axis *while molten* during the draw, freezing a rotating axis orientation into relaxed glass. | Rotates the **intrinsic linear axis** (ellipticity + thermal stress) at rate $\xi$. No circular birefringence. |
| **Mechanical twist** $\tau_m(z)$ | The *solid* fiber twisted about its axis after manufacture (torsional shear stress). | **Circular** birefringence via the photoelastic effect (rate $g\,\tau_m$), and co-rotates the intrinsic linear axis with the cross section. |

So the full generator is built as three additive contributions in the Bishop frame:

* **Bending / axial tension** — linear birefringence of magnitude $\Delta\beta_\text{bend}$ at orientation $\varphi_\text{bend}(z) = \int_0^z \tau_g\,dz'$ (queryable as `torsion_phase(path, z)`), reducing to the fixed $(\kappa, 0)$ axis on any planar path.
* **Core ellipticity + asymmetric thermal stress** — linear birefringence of summed magnitude $\Delta\beta_\text{ellip} + \Delta\beta_\text{stress}$ at orientation $\varphi_\text{int}(z) = \xi_0 + \phi_\xi(z) + \phi_{\tau_m}(z)$, where $\xi_0$ is the frozen ellipse angle (`ellipticity_axis_angle`) and $\phi_{\tau_m}$ is the accumulated mechanical-twist phase (`twist_phase(path, z)`).
* **Mechanical twist** — circular birefringence $\Delta\beta_c = g(\lambda, T)\,\tau_m$, the only source of the antisymmetric (optical-activity) term.

Writing these out, with $c_b = \cos 2\varphi_\text{bend}$, $s_b = \sin 2\varphi_\text{bend}$, and $c_i, s_i$ likewise for $\varphi_\text{int}$,

$$ K = \tfrac{i}{2}\Delta\beta_\text{bend}\!\left( \begin{array}{cc} c_b & s_b \\ s_b & -c_b \end{array} \right) + \tfrac{i}{2}(\Delta\beta_\text{ellip}+\Delta\beta_\text{stress})\!\left( \begin{array}{cc} c_i & s_i \\ s_i & -c_i \end{array} \right) + \tfrac{1}{2}\Delta\beta_c\!\left( \begin{array}{cc} 0 & -1 \\ 1 & 0 \end{array} \right). $$

These map onto `linear_birefringence_generator` and `circular_birefringence_generator` in [`fiber-path.jl`](../src/fiber/fiber-path.jl). The cross-section returns each birefringence *magnitude* (see [`step-index.jl`](../src/fiber-cross-section/step-index.jl)); the fiber layer supplies the orientation from the three accumulated phases. Each birefringence is evaluated at the segment-local temperature `local_temperature(fiber, z)`, so a `:T_K` excursion shifts the optics (notably $\Delta\beta_\text{stress} \propto |T_\text{soft} - T|$), not only the geometry.

If there were no spinning, twist, or torsion (a planar, unspun fiber with its ellipse along $N$), this collapses to

$$ \epsilon = \bar{\epsilon} + \epsilon_i \left( \begin{array}{cc} 1 & 0 \\ 0 & -1 \end{array} \right), $$

i.e. a fixed linear retarder, with the circular term present only under mechanical twist.