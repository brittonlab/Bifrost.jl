from pathlib import Path

import numpy as np
from bifrost import start

jl = start()

# Load Julia file
jl.include(str(Path(__file__).with_suffix(".jl")))
Zoo = jl.Zoo


def main():
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

    # Vector of animals -> weekly intake (kg)
    animals = jl.Vector[Zoo.Animal]([nala, kibo, mei])
    print("Weekly intake (kg):", np.array(Zoo.weekly_intake(animals)))
    print(f"Herd total mass = {Zoo.herd_total_mass(animals)} kg")

    # Numpy matrix of daily intakes (rows = animals, cols = days) scaled in Julia
    intake = np.array(
        [
            [7.5, 7.4, 7.6, 7.5, 7.5, 7.5, 7.5],     # Nala
            [150.0, 148.0, 152.0, 149.0, 151.0, 150.0, 150.0],  # Kibo
            [12.0, 11.5, 12.5, 12.0, 12.0, 12.0, 12.0],         # Mei
        ]
    )
    scaled = Zoo.scale_intake(intake, 1.1)  # +10% feeding plan
    print("Scaled intake (kg/day):")
    print(np.array(scaled))


if __name__ == "__main__":
    main()
