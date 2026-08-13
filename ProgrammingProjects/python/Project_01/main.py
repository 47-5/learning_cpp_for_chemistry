import argparse
from pathlib import Path

from chemistry import (
    GeometryAnalyzer,
    Molecule,
    RotationalAnalyzer,
)


DEFAULT_INPUT = (
    Path(__file__).resolve().parents[2]
    / "resources"
    / "Project_01"
    / "input"
    / "allene.dat"
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Project 1 molecular analysis"
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        type=Path,
        default=DEFAULT_INPUT,
        help="geometry file (default: allene.dat)",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    molecule = Molecule.from_file(arguments.input_file)
    geometry = GeometryAnalyzer(molecule)
    rotation = RotationalAnalyzer(molecule)

    print(f"Number of atoms: {molecule.atom_count}")
    print("Input Cartesian coordinates:")
    for atom in molecule.atoms:
        x, y, z = atom.position
        print(
            f"{atom.atomic_number} "
            f"{x:.12f} {y:.12f} {z:.12f}"
        )

    print("Interatomic distances (bohr):")
    for i in range(molecule.atom_count):
        for j in range(i):
            print(f"{i} {j} {geometry.distance(i, j)}")

    print("Bond angles:")
    for i in range(molecule.atom_count):
        for j in range(i + 1, molecule.atom_count):
            for k in range(j + 1, molecule.atom_count):
                if (
                    geometry.distance(i, j) <= 4.0
                    and geometry.distance(j, k) <= 4.0
                ):
                    print(
                        f"{i} {j} {k} "
                        f"{geometry.angle_degrees(i, j, k):.6f}"
                    )

    print("Out-of-plane angles:")
    for i in range(molecule.atom_count):
        for j in range(molecule.atom_count):
            for k in range(molecule.atom_count):
                for l in range(molecule.atom_count):
                    if len({i, j, k, l}) < 4 or j <= l:
                        continue

                    if (
                        geometry.distance(i, k) < 4.0
                        and geometry.distance(j, k) < 4.0
                        and geometry.distance(k, l) < 4.0
                    ):
                        angle = (
                            geometry.out_of_plane_angle_degrees(
                                i, j, k, l
                            )
                        )
                        print(f"{i} {j} {k} {l} {angle:.6f}")

    print("Torsional angles:")
    for i in range(molecule.atom_count):
        for j in range(i):
            for k in range(j):
                for l in range(k):
                    if (
                        geometry.distance(i, j) < 4.0
                        and geometry.distance(j, k) < 4.0
                        and geometry.distance(k, l) < 4.0
                    ):
                        angle = (
                            geometry.torsional_angle_degrees(
                                i, j, k, l
                            )
                        )
                        print(f"{i}-{j}-{k}-{l} {angle:.6f}")

    print("Center of mass:")
    print(rotation.center_of_mass())
    print("Moment of inertia tensor:")
    print(rotation.inertia_tensor())
    print("Principal moments of inertia:")
    print(rotation.principal_moments())
    print("Rotor type:")
    print(rotation.rotor_type())
    print("Principal moments of inertia (amu * AA^2):")
    print(rotation.principal_moments_amu_angstrom2())
    print("Principal moments of inertia (g * cm^2):")
    print(rotation.principal_moments_g_cm2())
    print("Rotational constants (MHz):")
    print(rotation.rotational_constants_mhz())
    print("Rotational constants (cm^-1):")
    print(rotation.rotational_constants_cm_inverse())


if __name__ == "__main__":
    main()
