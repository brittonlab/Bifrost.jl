# Deriving the Generators

BIFROST tries to solve the differential equation

```math
\frac{\partial J(z)}{\partial z} = K(z) J(z)
```

where $J(z)$ is the Jones matrix for our system (which also depends on e.g. $\omega$, the angular frequency of the light), $z$ is the position coordinate along the fiber, and $K(z)$ is a "generator." The generator contains all of the fundamental physics in the problem, so this document shows how to derive it.

Our formulation is very similar to that of Przhiyalkovsky et al., "Polarization Dynamics of Light Propagating in Bent Spun Birefringent Fiber," *Journal of Lightwave Technology* 24(38), 6879-6885 (2020), DOI: 10.1109/jlt.2020.3017795.

## Formulating the Equation

Suppose our fiber is laid along the $z$-axis. Propagation of monochromatic light in the fiber is dictated by the wave equation

```math
\frac{\partial^2\vec{E}}{\partial z^2} = -k_0^2 \epsilon \vec{E}.
```

Here $k_0 = \omega/c$ and $\epsilon$ is the permittivity tensor. This is what contains all of the interesting physics; we'll come back to it in a moment.

We approximate the waves as transverse, so $\vec{E} = (E_x, E_y)$ (i.e. only two components). There is a common phase gained by both polarization modes of the form $e^{-i k_0 \bar{n} z}$ where $\bar{n}$ is an average refractive index. Let's transform this away by defining $E_{x,y} = \mathcal{E}_{x,y} e^{-ik_0 \bar{n}z}$. If the anisotropy is low ($\epsilon_i, \epsilon_c \ll \bar{\epsilon}$), then this acts like a slowly-varying envelope approximation and, upon substitution into our wave equation, we can ignore the second-order derivatives of $\mathcal{E}_{x,y}$, which turns our wave equation into a Schrödinger-like equation:

```math
2i\bar{n} \frac{\partial \mathcal{E}}{\partial z} = k_0(\epsilon - \bar{n}^2)\mathcal{E}.
```

So the "generator" we're interested in is

```math
K(z) = -\frac{i}{2} \frac{k_0}{\bar{n}} (\epsilon - \bar{n}^2).
```

## The Permittivity Tensor

The permittivity tensor has contributions from every birefringence mechanism available. It explicitly describes how the two modes gain phase, so we can use our physical intuition to write it down. We must pick a basis for doing so; let's choose the $N,B$ basis where $N$ and $B$ are the normal and binormal of the Frenet-Serret frame. Then, as a starting point, we can write

```math
\epsilon = \bar{\epsilon} + i\epsilon_c \left( \begin{array}{cc} 0 & -1 \\ 1 & 0 \end{array} \right) + \epsilon_i \left( \begin{array}{cc} \cos 2\xi z & \sin 2\xi z \\ \sin 2\xi z & -\cos 2\xi z \end{array} \right).
```

Here $\bar{\epsilon} = \bar{n}^2$ is the average permittivity, $\epsilon_c$ is the circular birefringence, and $\epsilon_i$ is the intrinsic linear birefringence, which is rotating with the spinning of the fiber at spin rate $\xi$. Writing $\xi = 2\pi/L_s$ with a single spin pitch $L_s$ is the constant-rate special case used here for clarity. In general the spin rate is an arbitrary function of arc length, $\xi = \xi(z)$, and the implementation treats it that way: `ResolvedSpinningRate.rate` in `src/geometry/path-geometry.jl` is a `Union{Float64, Function}`, and the run integrator (`_integrate_rate`) takes the analytic branch for a constant rate and adaptive Gauss–Kronrod quadrature for a function rate. Replace $2\xi z$ below with $2\int_0^z \xi(z')\,dz'$ for the function-valued case.

This is a only starting point because the choice of the $N,B$ frame doesn't really mean anything. Without the curvature of bending, the direction of the normal vector is arbitrary. In this case, the normal vector direction is evidently set by the fact that the first component of the Jones vector here is along the slow axis and the second component is along the fast axis. Thus the normal vector is set along the slow axis at $z=0$ (and doesn't rotate with the spinning). That's a loose constraint because we sort of don't care about these overall rotations for our problems where we'll be sending in light of all polarizations. 

Where things get really complicated is when we introduce bending, because now we have to keep track of the linear birefringence direction and the spinning relative to the bend plane, and the curvature means the normal vector is rotating for actual physical reasons. **To be continued...**

If there were no spinning, then we would say 

```math
\epsilon = \bar{\epsilon} + i\epsilon_c \left( \begin{array}{cc} 0 & -1 \\ 1 & 0 \end{array} \right) + \epsilon_i \left( \begin{array}{cc} 1 & 0 \\ 0 & -1 \end{array} \right).
```

## Assembling the Fiber Generator

While the full bent-spun derivation above is still being worked out, the
implementation assembles the generator from independently-modeled mechanisms.
From here on we write the propagation coordinate as the arc length $s$ along the
centerline (the same role $z$ plays above).

The original BIFROST-style sliced approach calculates

```math
J_{\mathrm{total}}=\prod_i J_i,
```

with matrix order matching the order light encounters along the fiber. That
approach is simple, but it is difficult to attach a meaningful error bound when
linear birefringence, spinning, and other non-commuting terms vary along the
path.

The Julia implementation instead assembles a local generator:

```math
K(s,\omega)=K_{\mathrm{bend}}(s,\omega)+K_{\mathrm{spin}}(s,\omega).
```

The bending contribution comes from path curvature. For a local bend radius
`R(s)`, the implemented perturbation uses the bending birefringence response
from `src/fiber/fiber-cross-section.jl`; in the simplest stress model the
magnitude scales like `1/R(s)^2`.

The spinning contribution uses the total frame rotation rate:

```math
\tau_{\mathrm{path}}(s)=\tau_{\mathrm{geom}}(s)+\Omega(s).
```

Here `geometric_torsion(path, s)` comes from the centerline, while
`spinning_rate(path, s)` comes from resolved `Spinning` metadata.

The same decomposition exists for the frequency derivative:

```math
K_\omega(s,\omega)
=K_{\mathrm{bend},\omega}(s,\omega)+K_{\mathrm{spin},\omega}(s,\omega).
```

That keeps ordinary Jones propagation and DGD sensitivity propagation aligned:
both use the same `Fiber`, wavelength, breakpoint partition, and adaptive
integration strategy.

## Numerical Propagation

The exponential midpoint step is

```math
J_{n+1}=\exp\!\left(hK(s_n+h/2)\right)J_n.
```

It is useful here because the solution of a constant-coefficient matrix ODE is
exactly an exponential. Step by step, the method preserves the multiplicative
structure of Jones propagation. Under the lossless assumption, after removing
common phase, the Jones matrices live in `SU(2)`.

The adaptive controller uses step doubling:

- take one full step of size `h`,
- take two half steps of size `h/2`,
- compare the two results using `phase_insensitive_error`,
- accept or reject the step,
- update `h` with a cubic-root controller because the estimate scales as
  `O(h^3)`.

Path breakpoints come from the built path:

```julia
fiber_breakpoints(fiber) = breakpoints(fiber.path)
```

Those breakpoints include path segment boundaries and resolved spinning-run
boundaries. `propagate_fiber` calls `propagate_piecewise`, which integrates
independently over each smooth interval.

The 2x2 Jones exponential uses a closed form based on Cayley-Hamilton. The
implementation also factors out small numerical trace drift before applying the
traceless formula:

```math
\exp(A)=\exp(\operatorname{tr}(A)/2)
\left[\cosh(\mu)I+\operatorname{sinhc}(\mu)\tilde A\right],
\qquad \mu^2=-\det(\tilde A).
```

Here `sinhc(mu) = sinh(mu) / mu`, with a Taylor branch near zero. The full
derivation is given in the [appendix](@ref cayley-hamilton).

## DGD Sensitivity Propagation

The finite-difference DGD estimate used by the legacy implementation has the
form

```math
\partial_\omega J
\approx \frac{J(\omega+\Delta\omega)-J(\omega)}{\Delta\omega}.
```

The Julia propagator instead integrates the sensitivity matrix
`G = partial_omega J` directly:

```math
\frac{dJ}{ds}=KJ,\qquad J(s_0)=I,
```

```math
\frac{dG}{ds}=K_\omega J+KG,\qquad G(s_0)=0.
```

At the output, the PMD generator is

```math
H_{\mathrm{PMD}}=-iJ^{-1}G.
```

`output_dgd(J, G)` returns the eigenvalue spread of that generator. For 2x2
MCM-valued matrices, `output_dgd_2x2(J, G)` computes the same spread using a
closed-form Hermitian 2x2 formula, avoiding `LinearAlgebra.eigvals`.

The coupled sensitivity step is implemented using a closed-form Frechet
derivative of the 2x2 exponential:

```math
\exp\!\left(h\begin{bmatrix}K & K_\omega\\0 & K\end{bmatrix}\right)
=
\begin{bmatrix}E & F\\0 & E\end{bmatrix}.
```

This is implemented by `exp_block_upper_triangular_2x2`. Avoiding generic 4x4
`LinearAlgebra.exp` is important for MCM compatibility and is also faster for
ordinary `Float64` cases.

## [Appendix: Cayley-Hamilton 2x2 Exponential](@id cayley-hamilton)

For any 2x2 matrix `A`, Cayley-Hamilton gives

```math
A^2-\operatorname{tr}(A)A+\det(A)I=0.
```

If `A` is traceless, then

```math
A^2=-\det(A)I.
```

Define `mu^2 = -det(A)`. Then

```math
A^2=\mu^2 I,\qquad A^3=\mu^2 A,\qquad A^4=\mu^4 I,\ldots
```

The exponential series splits into even and odd powers:

```math
e^A
=\sum_{k=0}^{\infty}\frac{A^{2k}}{(2k)!}
 +\sum_{k=0}^{\infty}\frac{A^{2k+1}}{(2k+1)!}.
```

Using `A^(2k) = mu^(2k) I` and `A^(2k+1) = mu^(2k) A`:

```math
e^A
=
\left(\sum_{k=0}^{\infty}\frac{\mu^{2k}}{(2k)!}\right)I
+
\left(\sum_{k=0}^{\infty}\frac{\mu^{2k}}{(2k+1)!}\right)A.
```

Therefore

```math
e^A=\cosh(\mu)I+\frac{\sinh(\mu)}{\mu}A,
\qquad \mu^2=-\det(A).
```

At `mu = 0`, interpret `sinh(mu) / mu -> 1`, so `e^A = I + A`.
