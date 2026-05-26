import bifrost as bf

xs = bf.FiberCrossSection(
    bf.GermaniaSilicaGlass(0.036),
    bf.GermaniaSilicaGlass(0.0),
    8.2e-6,
    125e-6,
)

print(xs.core_diameter)

# spec = jl.PathSpecBuilder()
# jl.straight_b(spec, length=0.1)
# fiber = jl.Fiber(jl.build(spec), cross_section=xs)

# J, stats = jl.propagate_fiber(
#     fiber,
#     **{"λ_m": 1550e-9, "verbose": False},
# )