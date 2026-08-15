"""Load a Project 3 case and prepare an RHF calculation."""

from pathlib import Path

from chemistry.electronic_structure import (
    RHFCalculator,
    SCFOptions,
    load_ao_integrals,
)


INPUT_DIRECTORY = (
    Path(__file__).resolve().parents[2]
    / "resources"
    / "Project_03"
    / "input"
    / "h2o"
    / "STO-3G"
)


def main(input_directory: Path = INPUT_DIRECTORY) -> None:
    integrals = load_ao_integrals(input_directory)
    calculator = RHFCalculator(
        integrals,
        SCFOptions(
            energy_tolerance=1.0e-10,
            density_tolerance=1.0e-8,
            max_iterations=100,
        ),
    )

    print(f"Input directory: {input_directory}")
    print(f"Basis functions: {integrals.n_basis}")
    print(f"Electrons: {integrals.n_electrons}")
    print(f"Occupied RHF orbitals: {calculator.n_occupied}")
    print(
        "Nuclear repulsion: "
        f"{integrals.nuclear_repulsion:.12f} Eh"
    )
    result = calculator.run()

    print()
    print("SCF iterations")
    print(
        f"{'Iter':>4}"
        f" {'E(elec) / Eh':>20}"
        f" {'E(total) / Eh':>20}"
        f" {'Delta E / Eh':>14}"
        f" {'Delta D':>14}"
    )
    print("-" * 78)

    for iteration in result.history:
        print(
            f"{iteration.index:4d}"
            f" {iteration.electronic_energy:20.12f}"
            f" {iteration.total_energy:20.12f}"
            f" {iteration.energy_change:14.6e}"
            f" {iteration.density_change:14.6e}"
        )

    print()
    print(f"Converged: {result.converged}")
    print(f"SCF iterations: {result.history[-1].index}")
    print(
        "Final electronic energy: "
        f"{result.electronic_energy:.12f} Eh"
    )
    print(
        "Nuclear repulsion: "
        f"{integrals.nuclear_repulsion:.12f} Eh"
    )
    print(f"Final total energy: {result.total_energy:.12f} Eh")


if __name__ == "__main__":
    main()
