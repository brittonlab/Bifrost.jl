using Bifrost

xs = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like",
)

# Single Subpath: lead-in straight + 90° bend (axis_angle=0) + lead-out straight.
# After the quarter bend on +z at R=0.05, the local frame's +z aligns with the
# global +x axis, so the lead-out straight runs along +x.
# seal! ends the Subpath at its natural exit with no terminal connector bending.
sb = SubpathBuilder(); start!(sb)
straight!(sb; length = 0.5, meta = [Nickname("lead-in")])
bend!(sb;     radius = 0.05, angle = π / 2, meta = [Nickname("90 deg bend")])
straight!(sb; length = 0.5, meta = [Nickname("lead-out")])
seal!(sb)

fiber = Fiber(build(sb); cross_section = xs, T_ref_K = 297.15)

J, stats = propagate_fiber(fiber; λ_m = 1550e-9)

println("J =")
display(J)
println()
println("intervals = ", length(stats))
