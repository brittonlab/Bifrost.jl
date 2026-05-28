import math
import sys
from pathlib import Path

# Make `src/wrapper.py` importable without requiring PYTHONPATH=src on the command
# line, so `uv run python docs/juliacall-demo.py` works from a clean checkout.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import wrapper  # noqa: E402

# Python is the workspace; Julia modules are guest libraries. Tab-completion on
# `Zoo` lists only Zoo's exports (Animal, describe, ...), not the ~1000 names a
# raw juliacall module inherits from Base. Likewise for `Bifrost`.
Zoo = wrapper.wrap("Zoo", "docs/juliacall-demo.jl")
Bf = wrapper.wrap("Bifrost")


def main():
    print("\n \n")
    use_zoo()
    print("\n \n")
    use_bifrost()


def use_zoo():
    print("This method uses the Julia Zoo module. Welcome to the zoo! \n")
    # Build a few animals
    nala = Zoo.Animal("Nala", "lion", "F", 180.0, 7.5)
    kibo = Zoo.Animal("Kibo", "elephant", "M", 5400.0, 150.0)
    mei = Zoo.Animal("Mei", "giant panda", "F", 95.0, 12.0)

    print(Zoo.describe(3.14), "\n")    # dispatch on Number
    print(Zoo.describe(nala))    # dispatch on Animal
    print(Zoo.describe(kibo))
    print(Zoo.describe(mei), "\n")

    print(f"Nala food fraction = {Zoo.daily_food_fraction(nala):.4f}")
    print(f"Kibo food fraction = {Zoo.daily_food_fraction(kibo):.4f}", "\n")

    # Python list of animals -> Julia iterates over it directly; no raw jl handle needed.
    animals = [nala, kibo, mei]
    print("Weekly intake (kg):", np.array(Zoo.weekly_intake(animals)))
    print(f"Herd total mass = {Zoo.herd_total_mass(animals)} kg")

    intake = np.array(
        [
            [7.5, 7.4, 7.6, 7.5, 7.5, 7.5, 7.5],                # Nala
            [150.0, 148.0, 152.0, 149.0, 151.0, 150.0, 150.0],  # Kibo
            [12.0, 11.5, 12.5, 12.0, 12.0, 12.0, 12.0],         # Mei
        ]
    )
    scaled = Zoo.scale_intake(intake, 1.1)  # +10% feeding plan
    print("Scaled intake (kg/day):")
    print(np.array(scaled))


def use_bifrost():
    print("This method uses the Julia Bifrost module. \n")
    xs = Bf.FiberCrossSection(
        Bf.GermaniaSilicaGlass(0.036),
        Bf.GermaniaSilicaGlass(0.0),
        8.2e-6,
        125e-6,
    )

    # `straight!` and `bend!` end in `!` (Julia mutator convention), which Python
    # syntax can't spell directly — getattr is the standard juliacall escape. The
    # `meta=` kwarg expects `AbstractVector{<:AbstractMeta}`; juliacall converts a
    # Python list to `PyList{Any}` which fails the type constraint, so the demo
    # omits per-segment Nicknames here.
    straight = getattr(Bf, "straight!")
    bend = getattr(Bf, "bend!")

    spec = Bf.PathSpecBuilder()
    straight(spec, length=0.08)
    bend(spec, radius=0.06, angle=math.pi / 2)
    straight(spec, length=0.06)

    path = Bf.build(spec)
    fiber = Bf.Fiber(path, cross_section=xs)

    J = Bf.propagate_fiber(fiber, λ_m=1550e-9, verbose=False)[0]
    print("Final Jones matrix:")
    print(np.array(J))

    # Sample the path centerline and pull it into numpy so the plot is a pure
    # Python-workspace artifact (no Julia plotting backend, no HTML output).
    path_sample = Bf.sample_path(path, path.spec.s_start, path.s_end, fidelity=3.0)
    xyz = np.array([[float(s.position[i]) for i in range(3)] for s in path_sample.samples])

    # Popup window via matplotlib's default GUI backend (TkAgg/QtAgg). Not a browser.
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], color="tab:blue", linewidth=1.5)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("z (m)")
    ax.set_title("Bifrost fiber path: lead-in / 90° bend / spacer")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
