# Deriving the Generators

BIFROST tries to solve the differential equation

$$ \frac{\partial J(z)}{\partial z} = K(z) J(z) $$

where $J(z)$ is the Jones matrix for our system (which also depends on e.g. $\omega$, the angular frequency of the light), $z$ is the position coordinate along the fiber, and $K(z)$ is a "generator." The generator contains all of the fundamental physics in the problem, so this document shows how to derive it.

Our formulation is very similar to that of Przhiyalkovsky et al., "Polarization Dynamics of Light Propagating in Bent Spun Birefringent Fiber," *Journal of Lightwave Technology* 24(38), 6879-6885 (2020), DOI: 10.1109/jlt.2020.3017795.

## Formulating the Equation

Suppose our fiber is laid along the $z$-axis. Propagation of monochromatic light in the fiber is dictated by the wave equation

$$ \frac{\partial^2\vec{E}}{\partial z^2} = -k_0^2 \epsilon \vec{E}. $$

Here $k_0 = \omega/c$ and $\epsilon$ is the permittivity tensor. This is what contains all of the interesting physics; we'll come back to it in a moment.

We approximate the waves as transverse, so $\vec{E} = (E_x, E_y)$ (i.e. only two components). There is a common phase gained by both polarization modes of the form $e^{-i k_0 \bar{n} z}$ where $\bar{n}$ is an average refractive index. Let's transform this away by defining $E_{x,y} = \mathcal{E}_{x,y} e^{-ik_0 \bar{n}z}$. If the anisotropy is low ($\epsilon_i, \epsilon_c \ll \bar{\epsilon}$), then this acts like a slowly-varying envelope approximation and, upon substitution into our wave equation, we can ignore the second-order derivatives of $\mathcal{E}_{x,y}$, which turns our wave equation into a Schr\"odinger-like equation:

$$ 2i\bar{n} \frac{\partial \mathcal{E}}{\partial z} = k_0(\epsilon - \bar{n}^2)\mathcal{E}. $$

So the "generator" we're interested in is

$$ K(z) = -\frac{i}{2} \frac{k_0}{\bar{n}} (\epsilon - \bar{n}^2). $$

## The Permittivity Tensor

The permittivity tensor has contributions from every birefringence mechanism available. It explicitly describes how the two modes gain phase, so we can use our physical intuition to write it down. We must pick a basis for doing so; let's choose the $N,B$ basis where $N$ and $B$ are the normal and binormal of the Frenet-Serret frame. Then, as a starting point, we can write

$$ \epsilon = \bar{\epsilon} + i\epsilon_c \left( \begin{array}{cc} 0 & -1 \\ 1 & 0 \end{array} \right) + \epsilon_i \left( \begin{array}{cc} \cos 2\xi z & \sin 2\xi z \\ \sin 2\xi z & -\cos 2\xi z \end{array} \right). $$

Here $\bar{\epsilon} = \bar{n}^2$ is the average permittivity, $\epsilon_c$ is the circular birefringence, and $\epsilon_i$ is the intrinsic linear birefringence, which is rotating with the spinning of the fiber at spin rate $\xi = 2\pi/L_s$, with $L_s$ the spin pitch.

This is a only starting point because the choice of the $N,B$ frame doesn't really mean anything. Without the curvature of bending, the direction of the normal vector is arbitrary. In this case, the normal vector direction is evidently set by the fact that the first component of the Jones vector here is along the slow axis and the second component is along the fast axis. Thus the normal vector is set along the slow axis at $z=0$ (and doesn't rotate with the spinning). That's a loose constraint because we sort of don't care about these overall rotations for our problems where we'll be sending in light of all polarizations. 

Where things get really complicated is when we introduce bending, because now we have to keep track of the linear birefringence direction and the spinning relative to the bend plane, and the curvature means the normal vector is rotating for actual physical reasons. **To be continued...**

If there were no spinning, then we would say 

$$ \epsilon = \bar{\epsilon} + i\epsilon_c \left( \begin{array}{cc} 0 & -1 \\ 1 & 0 \end{array} \right) + \epsilon_i \left( \begin{array}{cc} 1 & 0 \\ 0 & -1 \end{array} \right). $$