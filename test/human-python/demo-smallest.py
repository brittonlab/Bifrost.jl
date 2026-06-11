import bifrost_py as bf
#bf.info()

xs = bf.StepIndexCrossSection(
    bf.SilicaGermaniaGlass(0.036),
    bf.SilicaGermaniaGlass(0.0),
    8.2e-6,
    125.0e-6,
    manufacturer = "Corning",
    model_number = "SMF-like"
)
print(xs.manufacturer)

spec = bf.SubpathBuilder(); bf.start_b(spec)
bf.straight_b(spec, length = 0.5, meta=[bf.Nickname("lead-in")])
bf.bend_b(spec, radius = 0.05, angle = bf.pi / 2, meta=[bf.Nickname("90 deg bend")])
bf.straight_b(spec, length = 0.5, meta=[bf.Nickname("lead-out")])
bf.seal_b(spec)

seg = spec.segments[1]  # Access the SECOND segment... Python indexing, not Julia
print("The bend's nickname is ", bf.segment_nickname(seg))

fiber = bf.Fiber(bf.build(spec), cross_section=xs, T_ref_K=297.15)
J, stats = bf.propagate_fiber(fiber, lambda_m=1550e-9)

print("Jones matrix:")
print(J)
print(f"\nintervals: {len(stats)}")
print(f"First interval: {stats[0]}")