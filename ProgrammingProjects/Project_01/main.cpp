#include "molecule.h"

#include <exception>
#include <iomanip>
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr
            << "Usage: project01 <input-file>\n";

        return 1;
    }

    try {
        Molecule molecule =
            Molecule::from_file(argv[1]);

        std::cout
            << "Number of atoms: "
            << molecule.atomCount()
            << '\n';

        std::cout
            << "Input Cartesian coordinates:\n";

        std::cout
            << std::fixed
            << std::setprecision(12);

        for (
            std::size_t i = 0;
            i < molecule.atomCount();
            ++i
        ) {
            const Atom& atom = molecule.atom(i);

            std::cout
                << atom.atomic_number << ' '
                << atom.position.x() << ' '
                << atom.position.y() << ' '
                << atom.position.z() << '\n';
        }

        // 打印原子间距离
        molecule.print_distances();
        // 打印键角
        molecule.print_angles();
        // 打印面外角
        molecule.print_out_of_plane_angles();
        molecule.print_torsional_angles();
        // 打印质心
        Eigen::Vector3d center_of_mass = molecule.get_center_of_mass();
        std::cout << "center of mass:\n";
        std::cout << center_of_mass << '\n';
        // 打印moment of inertia tensor
        std::cout << "moment of inertia tensor:\n" << molecule.get_moment_of_inertia_tensor() << '\n';
        //principal moments of inertia.
        std::cout << "principal moments of inertia.\n" << molecule.get_principal_moments_of_inertia() << '\n';
        // 对称性
        std::cout << "rotor_type:\n" << molecule.get_rotor_type() << '\n';

        std::cout << "Principal moments of inertia (amu * AA^2):\n"
                  << molecule.get_principal_moments_amu_angstrom2() << '\n';
        std::cout << "Principal moments of inertia (g * cm^2):\n"
                  << std::scientific << molecule.get_principal_moments_g_cm2() << '\n';
        std::cout << std::fixed << "Rotational constants (MHz):\n"
                  << molecule.get_rotational_constants_mhz() << '\n';
        std::cout << "Rotational constants (cm^-1):\n"
                  << molecule.get_rotational_constants_cm_inverse() << '\n';

    }
    catch (const std::exception& error) {
        std::cerr
            << "Error: "
            << error.what()
            << '\n';

        return 1;
    }

    return 0;
}