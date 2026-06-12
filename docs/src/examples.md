# Examples

## Julia

From the repository root, start Julia with `julia --project=.`:

```julia
using Bifrost

xs = StepIndexCrossSection(
    GermaniaSilicaGlass(0.036),
    GermaniaSilicaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like",
)

spec = PathSpecBuilder()
straight!(spec; length = 0.5, meta = [Nickname("lead-in")])
bend!(spec; radius = 0.05, angle = pi / 2, meta = [Nickname("90 deg bend")])
straight!(spec; length = 0.5, meta = [Nickname("lead-out")])

path = build(spec)
fiber = Fiber(path; cross_section = xs, T_ref_K = 297.15)

J, stats = propagate_fiber(fiber; λ_m = 1550e-9, rtol = 1e-9, verbose = false)
```

For DGD:

```julia
J, G, stats = propagate_fiber_sensitivity(
    fiber;
    λ_m = 1550e-9,
    rtol = 1e-9,
    verbose = false,
)

dgd = output_dgd(J, G)
```

For Monte Carlo Measurements (MCM) paths, prefer `output_dgd_2x2(J, G)` over
`output_dgd(J, G)` because it avoids `eigvals`. See
[MCM compatibility](@ref mcm-compatibility) for the surrounding conventions.

## Python

The `wrapper.py` shim boots juliacall against this project and returns a
tab-completable view of any Julia module. Names listed in `dir(Bf)` are exactly
what that module `export`s — none of the ~1000 names a juliacall module
otherwise inherits from `Base`. Julia names ending in `!` use the `_b` suffix.

```python
import wrapper

Bf = wrapper.wrap("Bifrost")

xs = Bf.StepIndexCrossSection(
    Bf.GermaniaSilicaGlass(0.036),
    Bf.GermaniaSilicaGlass(0.0),
    8.2e-6,
    125e-6,
)

spec = Bf.PathSpecBuilder()
Bf.straight_b(spec, length=0.1)
fiber = Bf.Fiber(Bf.build(spec), cross_section=xs)

J, stats = Bf.propagate_fiber(fiber, λ_m=1550e-9, verbose=False)
```

For a fuller juliacall example, see
[`docs/juliacall-demo.py`](https://github.com/brittonlab/BIFROST/blob/main/docs/juliacall-demo.py)
and its Julia-side module
[`docs/juliacall-demo.jl`](https://github.com/brittonlab/BIFROST/blob/main/docs/juliacall-demo.jl).
To run python code in the uv environment use
`uv run python docs/juliacall-demo.py`.
