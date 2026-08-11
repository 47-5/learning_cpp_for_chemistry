from pathlib import Path

from molecule import Molecule



if __name__ == '__main__':

    data_path = Path('../input/allene.dat')
    m = Molecule.from_file(file_path=data_path)

    print(m)
    print(m.N)
    print(m.positions)
    print(m.atomic_numbers)
    m.print_distances()
    m.print_angles()
    m.print_out_of_plane_angles()
    m.print_torsional_angles()
    print(m.get_center_of_mass())
    print(m.get_moment_of_inertia_tensor())
    print(m.get_principal_moments_of_inertia())

    moments_bohr2, moments_angstrom2, moments_g_cm2 = (
        m.get_principal_moments_in_different_units()
    )

    constants_mhz, constants_cm_inverse = (
        m.get_rotational_constants()
    )

    print("amu * bohr^2:", moments_bohr2)
    print("amu * AA^2:", moments_angstrom2)
    print("g * cm^2:", moments_g_cm2)

    print("MHz:", constants_mhz)
    print("cm^-1:", constants_cm_inverse)
