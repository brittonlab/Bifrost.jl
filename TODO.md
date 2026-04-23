This list of tasks helps keep humans organized but also provides agents with an idea of
upcoming features. Agents should not start work on these without explicit authorization.

- [ ] add struct in path-integral.jl that reflects simulation parameters (eg `rtol`,
      `atol`, step controls)

- [x] Update AGENT.md to provide guidance on how to structure unit tests. Focus on
      clearly deliniating high-value physics motivated tests, test involving validation
      data, tests that accompany select simulation test runs and finally tests that help
      keep agents from wandering

- [x] Is there a julia convention for building a 3D path one bit at a time? I guess this
      will likey involve the Frenet-Serret frame. I'm open to the idea that there's
      already a convention for this. In that case I'd use multiple dispatch to express
      the detailed intent for my fiber based application.

- [x] implement for T in (JumpBy, JumpTo)
  - [x] properly implement sample_path for these

# quality of life changes for path-geometry-plot.jl
- [x] Make the scale factor for x,y,z the same in path-geometry-plot.jl 
- [x] Remove the "Local Frenet section" subplot in path-geometry.plot.jl.
- [x] Make the left-right scrub linear in s not linear in the sample number. 
- [x] Add an option to give a string-based nickname to a Segment. This would be
an optional argument for eg straight!. By default there is no nickname. If 
there is a nickname add a text label to the segment 
rendered by path-geometry-plot.jl. If possible the label should lie in the same plane as the segment. 

- [x] Add minimum radius of curvature optional parameter for jumpto and jumpby in path-geometry.jl.
- [x] Let's think if there's a more straightforward check for violation of min_bend_radius. 1) Are you using _hc_peak_curvature() to check? Based on its name it may be relevant. 2) Since HermiteConnector implements only a cubic spline does that bound the number of points that have to be sampled to conclusively determine that there is no violation of min_bend_radius? 
  - This didn't help. 


- [ ] In light of the recent addition of ARTHITECTURE.md and AGENT.md and README.md is
      there any refactoring that should take place? Are there any consequential
      contradictions between these new files and the existing organization, logic,
      represented physics or physics modeling in the codebase?

- [ ] Move the portion of fiber-path-plot.jl that relates to the 3D geometry of the fiber
      to to path-geometry-plot.jl. The 3D plotting should be informed by the geometries
      and geometric calculations in path-geometry.jl. It should include a movable plane
      that is moved by the mouse along the path length from start to finish. At each
      intermediate point the fernet frame coordinates should be graphically illustrated
      in a cross- section plot. Nothing in path-geometry-plot.jl should relate to optical
      fiber. Create some simple example paths in demo.jl that illustrate
      path-geometry-plot.jl.

  - [ ] test the refactoring from Sunday night related to above

- [ ] Move the portion of fiber-path-plot.jl related to plotting the Poincare sphere and
      move it into a new file called poincare-sphere-plot.jl.

- [ ] Restructure fiber-path-plot.jl to draw upon on path-geometry-plot.jl and
      poincare-sphere-plot.jl while still retaining its overall goals. That is,
      illustrating the transformation of an input polarization state as a function of
      distance s along the length of an optical fiber.



- [ ] Create a path-geometry.md that documents how it works. Add specific 
  illustratings for important features of path-geometry.jl. Each is described in
  .md and accompanied by code-generated .png that illustrate the points. The
  file defining these illustrations is called path-geometry-illustrated.jl. 
    - [ ] Show that the reference frame for adding additional segments does not 
    rotate with fiber twist.
    - [ ] Illustrate how the axis_angle option is defined and how it
    changes segment addition.
    - [ ] Show how the orientation of the prior segment influences the orientation of a helix and the helix exit path.
    - [ ] Devise examples that illustrate how segments respond to shrinkage and
    contrast it with 

- [ ] Do I want to change the implementation of path-geometry.jl to 
permit all the parametric parameters to be Functions? This would parallel
how it's done in fiber-path.jl (eg in BendSegment)
  - This could be motivated by my desire to use the MonteCarloMeasurements.jl 
  type Particles or StaticParticles for the following parameters temperature, path geometry parameters, core_diameter, cladding_diameter, core_noncircularity, axial tension. The MonteCarloMeasurements.jl docs discusses supporting new 
  functions here https://baggepinnen.github.io/MonteCarloMeasurements.jl/stable/overloading/

- [ ] change the implementation of twist! bend! etc in fiber-path.jl to leverage
  the same in path-geometry.pl. Note that the contents of path-geometry.jl ought 
  not be  loaded into the global namespace of fiber-path as fiber-path 
  implements its own twist! bend! etc variants

- [ ] The following relates to path-geometry.jl, fiber-cross-section.jl and fiber-path.jl. I may need to refactor some things due to rewriting path-geometry.jl. 
I need to decide where in the stack to store temperature as a function of s
  - temperature directly relates to some macroscopic degrees of freedom
  represented in path-geometry. 
    - shrinkage depends on the material property poisson_ratio(T) 
    - core_noncircularity birefringence 
    - axial_tension birefringence
    - twisting biregringence  
    - asymmetric thermal stress
  - temperature influences some intrinsic material properies (that don't 
  depend on geometry) for example 
    - core_refractive_index
    - cladding_refractive_index
    - mode_terms
    - chromatic_dispersion_parameter



- [ ] fiber-path.jl depends on path-geometry.jl and r-path.jl to fiber-