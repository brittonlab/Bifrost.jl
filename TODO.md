This list of tasks helps keep humans organized but also provides agents with an idea of
upcoming features. Agents should not start work on these without explicit authorization.

- [x] Update AGENT.md to provide guidance on how to structure unit tests. Focus on
      clearly deliniating high-value physics motivated tests, test involving validation
      data, tests that accompany select simulation test runs and finally tests that help
      keep agents from wandering

- [x] Is there a julia convention for building a 3D path one bit at a time? I guess this
      will likey involve the Frenet-Serret frame. I'm open to the idea that there's
      already a convention for this. In that case I'd use multiple dispatch to express
      the detailed intent for my fiber based application.

- [x] The method uncovered_intervals in path-geometry.jl seems redundant since fiber segments 
  are built piecewise. My construction it's not possible for there to be gaps. Don't do this.... they are needed. 


- [x] Move the portion of fiber-path-plot.jl that relates to the 3D geometry of the fiber
      to to path-geometry-plot.jl. The 3D plotting should be informed by the geometries
      and geometric calculations in path-geometry.jl. It should include a movable plane
      that is moved by the mouse along the path length from start to finish. At each
      intermediate point the fernet frame coordinates should be graphically illustrated
      in a cross- section plot. Nothing in path-geometry-plot.jl should relate to optical
      fiber. Create some simple example paths in demo.jl that illustrate
      path-geometry-plot.jl.

- [ ] Update demo.jl to use physical birefringences if possible.

- [ ] TODO fix the MCM demo in demo.jl 4/28 task

##################
### HIGH LEVEL ###
##################

- [x] Ensure compatability of this codebase with Windows, Linux and MacOS. JWB See linter created 
for brittonlab/stablepoint

- [ ] Pythonic interface to the julia code.
   - auto-generated and follows the julia API almost verbatum 

- [ ] new folder structure... proposal
      ```bash
      src 
      src/nonlinear
      src/geometry
      src/fiber
      test
      test/legacy-python
      test/human
      doc
      Manifest.toml
      Project.toml
      ```

- [ ] Confirm that spinning and twisting are separated/named correctly everywhere in the code, e.g.
path-geometry defines a Twist meta that seems like it should be a Spin meta.

- [ ] Read and edit all docstrings as appropriate; for instance path-geometry.jl's intro docs
reference a fiber-path-shrinkage.jl that doesn't exist.

#################
#### TESTING ####
#################

- [ ] PB tests that confirm julia port reproduces representative tests cases of legacy-python

- [ ] PB port python unit tests to julia (eg reproducing fiber-properties)

- [ ] Add unit tests for fiber paddles that reproduce Thorlabs website. 

#######################
#### path-integral ####
#######################

- [ ] add struct in path-integral.jl that reflects simulation parameters (eg `rtol`,
      `atol`, step controls)

- [ ] Resume work on branch benchmark-integrate. 

- [ ] Make it possible for users to use adaptive-step-doubling graphical interface for 
understanding/debugging. Do we even want this type of feature? 

#########################
#### path-geometry ######
#########################
JWB focuses on these features/bugs. 

- [ ] high level reframing to make first and last point of each Subpath be absolute (jumpto like)

- [ ] in `demo1.jl` ... `helix-mcm-spinning.html` discontinuous jump in int \tau_spin ds

- [ ] in `demo2.jl` `modify-jumpby-drift-2d.html` too small radius... should be caught, right?

- [ ] in `demo2.jl` confirm that length preservation is physical

- [ ] Remove the "Local Frenet section" subplot in path-geometry.plot.jl.

- [ ] Move the portion of fiber-path-plot.jl related to plotting the Poincare sphere and
      move it into a new file called poincare-sphere-plot.jl.

- [ ] Restructure fiber-path-plot.jl to draw upon on path-geometry-plot.jl and
      poincare-sphere-plot.jl while still retaining its overall goals. That is,
      illustrating the transformation of an input polarization state as a function of
      distance s along the length of an optical fiber.

- [ ] Is this what we want? Piecewise bend! loops don't accumulate geometric spinning in 
      total_spinning — but they shouldn't, because a BendSegment has geometric_torsion = 0 
      (circular arcs have zero torsion). A helix does have nonzero geometric_torsion, but 
      that's captured in geometric_torsion(seg, s), not in total_spinning.

- [ ] Verify the `Spinning` `is_continuous = true` carry-over semantics are what we
      want: when a `Spinning` meta has `is_continuous = true`, the resolver in
      path-geometry.jl computes its `phi_0` as `prev.phi_0 + ∫_0^{prev_run_length}
      prev.rate(s_local) ds_local` — i.e. the absolute phase at the start of the
      new run equals the accumulated phase at the end of the prior run. Confirm
      this is the intended physical meaning before any downstream consumer
      (polarization propagator, plotting, etc.) starts relying on `phi_0`.

- [ ] Add a `resolved::Bool` flag (or similar mechanism) to `AbstractMeta`
      so the system can certify that all meta on a `SubpathBuilt` /
      `PathBuilt` has been fully processed. Some meta is interpreted by code
      internal to `path-geometry*.jl` (e.g. `Spinning`); other meta is
      interpreted by external code (e.g. `MCMadd`/`MCMmul` in
      `fiber-path-modify.jl`). A `resolved` flag would let downstream
      consumers (or a `complete(::SubpathBuilt)` / `complete(::PathBuilt)`
      check) verify that no meta has been silently ignored before a built
      structure is treated as authoritative.

- [ ] Create a path-geometry.md that documents how it works. Add specific 
  illustratings for important features of path-geometry.jl. Each is described in
  .md and accompanied by code-generated .png that illustrate the points. The
  file defining these illustrations is called path-geometry-illustrated.jl. 
    - [ ] Show that the reference frame for adding additional segments does not 
    rotate with fiber spinning.
    - [ ] Illustrate how the axis_angle option is defined and how it
    changes segment addition.
    - [ ] Show how the orientation of the prior segment influences the orientation of a helix 
    and the helix exit path.
    - [ ] Devise examples that illustrate how segments respond to shrinkage and
    contrast it with 

- [ ] JWB Joe misunderstands that there is an additional source of birefringence due to physical 
    twisting (beyond factory-based spinning). Ask PB about this. 

- [ ] I want to make some modifications that focus on path-geometry.jl. Let's worry about the 
      downstream consequences of these changes later. 

      Currently SpinningOverlay is specified by s_start and length. I want to make changes so that 
      the start and end of each spinning is defined with respect to segment boundaries. There are 
      several ways  I can think of implementing this. 

      OPTION 1 :: use meta
      Create a new struct Spinning <: AbstractMeta with members
            rate::Function
            \phi_0::Real
            is_continuous::Bool
      In this approach Spinning meta is associated with the segment where a particular spinning 
      rate commences and continues until the end of the Path or until another segment has an 
      associated Spinning meta. The bool is_continuous specifies if the spinning phase remains 
      continuous with the prior Spinning specification. If is_continuous is True then \phi_0 
      shouldn't be specified. If is_continuous is False than \phi_0 must be specified as this is 
      the starting phase. The rate::Function must accept \phi as a parameter. 

      OPTION 2 :: Use zero-length Segment
      Create a new struct Spinning <: AbstractPathSegment with members
            rate::Function
            \phi_0::Real
            is_continuous::Bool
      In this approach, a Spinning segment is inserted into Path placed_segment as a zero-length 
      segment that demarks the start of a particular spinning specified by its member data. The 
      same spinning rate  continues until the end of the Path or until another Spinning segment 
      is added to the Path. The bool is_continuous specifies if the spinning phase remains 
      continuous with the prior Spinning specification. If is_continuous is True then \phi_0 
      shouldn't be specified. If is_continuous is False than \phi_0 must be specified as this 
      is the starting phase. The rate::Function must accept \phi as a parameter. 

      One detail common to both approaches is that if is_continuous is True the length of prior 
      segments is important in calculating the phase continuity at the boundary between old and new 
      spinning rates. 

      While it's outside the context of the current refactoring it's important to note that 
      properties of individual path segments can be modified using meta

      Please help me think which is the right approach of if there is another that's better. 


- [ ] TODO20260505 PlacedSegment.origin and PlacedSegment.frame are stored in 
  global coordinates, evolving forward as segments are added.   
  The "first and last point anchored globally" part of your 
  vision is correct, and "coordinate frame evolving per segment"
   is correct — the frame on segment N is the global frame at
  the start of segment N. But the storage is global, not
  Subpath-local. The path just happens to start at global
  (0,0,0) because there's no Subpath start_point concept yet.

  What the new plan currently says                              
   
  Same as today, except pos starts at collect(sub.start_point)  
  (still global, just anchored at the Subpath's start_point 
  instead of (0,0,0)). PlacedSegments still hold global         
  origin/frame.                                             

  What "true Subpath-local storage" would mean

  PlacedSegment.origin would be the offset from start_point     
  (Subpath-local), and PlacedSegment.frame would be a rotation
  matrix whose columns are local N/B/T axes (with local +z =    
  start_outgoing_tangent). Global queries would compute     
  start_point + R_start · ps.origin + R_start · ps.frame ·
  v_segment_local, where R_start is the rotation taking
  Subpath-local axes to global.

  This is genuinely stronger independence: a SubpathBuilt could 
  be "relocated" by mutating its Subpath's start_point without
  rebuilding, and the internal data is invariant to where in    
  space the Subpath sits. The cost is more careful build() math
  (rotating jumpto_point into local before solving the
  connector) and a query layer that does R_start · ...
  everywhere.

#################################
### High-Level Future Changes ###
#################################

- [ ] Eventually, web-based documentation like that of NumPy. A quick search suggests Documenter.jl
exists for this purpose and can automatically generate web-based documentation using GitHub Actions,
though I suspect we might want to write some of our own.

- [ ] PB: I want to eventually change the philosophy of the materials/birefringences part of the
code to be more accessible to users who want to add new fiber types, new birefringence responses,
etc. This may somewhat merge with what Prakriti is doing, but I'd like the material data to come
from JSON files or something clearly separated as "data", and I'd like the fiber-cross-section Julia
file to be well-documented enough that a user who hasn't learned Julia doesn't really need to in
order to add a new birefringence rersponse if they want to. Right now it's a disorganized morass of
methods and some of them won't work for all fiber types. I particularly have in mind graded-index
fiber, which should be a priority for us to add in the future since SMF-28 is likely graded-index
fiber.

- [ ] PB: Related to the above, I'd like to consider in the future turning a lot of hard-coded
things into somewhat more dynamic lists. For instance right now, the generator_K() method directly
adds bending plus twisting, and if someone wants the generator to also have, say, birefringence due
to axial magnetic fields, well then they have to go in and hard-modify the generator_K() function.
And a user might want to add this for only their fiber type or something. In the past I proposed
a BirefringenceRegistry and MaterialRegistry to deal with this. 

- [ ] PB: Remove fluorinated silica glass. I put it in BIFROST v0 to do some testing and it has
so many unsupported pieces to it (and I didn't carefully research what I did put in), so I don't
feel comfortable having it in. It could be a separate material if we separate fiber materials later.
