import wrapper as bf

bf.info()    # prints what got loaded

xs = bf.FiberCrossSection(
    bf.GermaniaSilicaGlass(0.036),
    bf.GermaniaSilicaGlass(0.0),
    8.2e-6, 125e-6,
)

spec = bf.PathSpecBuilder()
bf.straight_b(spec, length=0.5)
bf.bend_b(spec, radius=0.05, angle=3.14159 / 2)
bf.straight_b(spec, length=0.5)

fiber = bf.Fiber(bf.build(spec), cross_section=xs, T_ref_K=297.15)
J, stats = bf.propagate_fiber(fiber, λ_m=1550e-9)

print("Jones matrix:")
print(J)
print(f"\nintervals: {len(stats)}")