from pathlib import Path

from chemistry import Molecule
from chemistry.vibrational_analysis import VibrationalAnalysis

INPUT_DIRECTORY = (
    Path(__file__).resolve().parents[2]
    / "resources"
    / "Project_02"
    / "input"
)

GEOMETRY_FILE = INPUT_DIRECTORY / "h2o_geom.txt"
HESSIAN_FILE = INPUT_DIRECTORY / "h2o_hessian.txt"


def main(
    geometry_file: Path = GEOMETRY_FILE,
    hessian_file: Path = HESSIAN_FILE,
) -> None:
    molecule = Molecule.from_file(geometry_file)

    if not hessian_file.is_file():
        raise ValueError(
            "Cannot open Hessian file: "
            f"{hessian_file}"
        )

    dimension = 3 * molecule.atom_count
    print(f"Number of atoms: {molecule.atom_count}")
    print(f"Hessian dimension: {dimension} x {dimension}")

    # Step 2: Read the Cartesian Hessian into a NumPy array.
    vibrational_analysis = VibrationalAnalysis(molecule, hessian_file)
    print('hessian:')
    print(vibrational_analysis.hessian)

    # Step 3: 获取质权坐标下的Hessian
    mass_weighted_hessian = vibrational_analysis.get_mass_weighted_hessian()
    print('mass_weighted_hessian:')
    print(mass_weighted_hessian)

    # Step 4: 对角化F，从而得到特征值和normal_mode
    eigen_values, normal_modes = vibrational_analysis.get_normal_modes()
    print('eigen_values:')
    print(eigen_values)
    print('normal_modes:')
    print(normal_modes)
    # Step 5: Convert eigenvalues to spectroscopic wavenumbers.
    wavenumbers = (
        vibrational_analysis.get_wavenumbers_cm_inverse()
    )
    print("wavenumbers (cm^-1):")
    for mode_index, wavenumber in enumerate(wavenumbers):
        print(f"{mode_index:3d} {wavenumber:14.6f}")

if __name__ == "__main__":
    main()
