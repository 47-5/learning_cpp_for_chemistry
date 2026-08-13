import argparse
from pathlib import Path

from chemistry import Molecule


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Project 2 vibrational analysis"
    )
    parser.add_argument(
        "geometry_file",
        type=Path,
    )
    parser.add_argument(
        "hessian_file",
        type=Path,
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    molecule = Molecule.from_file(
        arguments.geometry_file
    )

    if not arguments.hessian_file.is_file():
        raise ValueError(
            "Cannot open Hessian file: "
            f"{arguments.hessian_file}"
        )

    dimension = 3 * molecule.atom_count
    print(f"Number of atoms: {molecule.atom_count}")
    print(f"Hessian dimension: {dimension} x {dimension}")

    # Step 2: Read the Cartesian Hessian into a NumPy array.


if __name__ == "__main__":
    main()
