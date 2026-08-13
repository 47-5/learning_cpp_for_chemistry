#include <chemistry/geometry_analyzer.h>
#include <chemistry/molecule.h>
#include <chemistry/rotational_analyzer.h>

#include <exception>
#include <iomanip>
#include <iostream>
#include <set>

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: project01 <input-file>\n";
        return 1;
    }

    try {
        const chemistry::Molecule molecule =
            chemistry::Molecule::from_file(argv[1]);
        const chemistry::GeometryAnalyzer geometry{molecule};
        const chemistry::RotationalAnalyzer rotation{molecule};

        std::cout
            << "Number of atoms: "
            << molecule.atom_count()
            << '\n'
            << "Input Cartesian coordinates:\n"
            << std::fixed
            << std::setprecision(12);

        for (const chemistry::Atom& atom : molecule.atoms()) {
            std::cout
                << atom.atomic_number << ' '
                << atom.position.x() << ' '
                << atom.position.y() << ' '
                << atom.position.z() << '\n';
        }

        std::cout << "Interatomic distances (bohr):\n";
        for (std::size_t i = 0; i < molecule.atom_count(); ++i) {
            for (std::size_t j = 0; j < i; ++j) {
                std::cout
                    << i << ' ' << j << ' '
                    << geometry.distance(i, j) << '\n';
            }
        }

        std::cout << "Bond angles:\n";
        for (std::size_t i = 0; i < molecule.atom_count(); ++i) {
            for (std::size_t j = i + 1;
                 j < molecule.atom_count();
                 ++j) {
                for (std::size_t k = j + 1;
                     k < molecule.atom_count();
                     ++k) {
                    if (geometry.distance(i, j) <= 4.0 &&
                        geometry.distance(j, k) <= 4.0) {
                        std::cout
                            << i << ' ' << j << ' ' << k << ' '
                            << std::setprecision(6)
                            << geometry.angle_degrees(i, j, k)
                            << '\n';
                    }
                }
            }
        }

        std::cout << "Out-of-plane angles:\n";
        for (std::size_t i = 0; i < molecule.atom_count(); ++i) {
            for (std::size_t j = 0; j < molecule.atom_count(); ++j) {
                for (std::size_t k = 0;
                     k < molecule.atom_count();
                     ++k) {
                    for (std::size_t l = 0;
                         l < molecule.atom_count();
                         ++l) {
                        if (std::set<std::size_t>{i, j, k, l}.size() < 4 ||
                            j <= l) {
                            continue;
                        }

                        if (geometry.distance(i, k) < 4.0 &&
                            geometry.distance(j, k) < 4.0 &&
                            geometry.distance(k, l) < 4.0) {
                            std::cout
                                << i << ' ' << j << ' '
                                << k << ' ' << l << ' '
                                << geometry.out_of_plane_angle_degrees(
                                       i, j, k, l
                                   )
                                << '\n';
                        }
                    }
                }
            }
        }

        std::cout << "Torsional angles:\n";
        for (std::size_t i = 0; i < molecule.atom_count(); ++i) {
            for (std::size_t j = 0; j < i; ++j) {
                for (std::size_t k = 0; k < j; ++k) {
                    for (std::size_t l = 0; l < k; ++l) {
                        if (geometry.distance(i, j) < 4.0 &&
                            geometry.distance(j, k) < 4.0 &&
                            geometry.distance(k, l) < 4.0) {
                            std::cout
                                << std::setw(2) << i << '-'
                                << std::setw(2) << j << '-'
                                << std::setw(2) << k << '-'
                                << std::setw(2) << l << ' '
                                << std::setw(10)
                                << geometry.torsional_angle_degrees(
                                       i, j, k, l
                                   )
                                << '\n';
                        }
                    }
                }
            }
        }

        std::cout
            << "Center of mass:\n"
            << rotation.center_of_mass() << '\n'
            << "Moment of inertia tensor:\n"
            << rotation.inertia_tensor() << '\n'
            << "Principal moments of inertia:\n"
            << rotation.principal_moments() << '\n'
            << "Rotor type:\n"
            << rotation.rotor_type() << '\n'
            << "Principal moments of inertia (amu * AA^2):\n"
            << rotation.principal_moments_amu_angstrom2() << '\n'
            << "Principal moments of inertia (g * cm^2):\n"
            << std::scientific
            << rotation.principal_moments_g_cm2() << '\n'
            << std::fixed
            << "Rotational constants (MHz):\n"
            << rotation.rotational_constants_mhz() << '\n'
            << "Rotational constants (cm^-1):\n"
            << rotation.rotational_constants_cm_inverse() << '\n';
    }
    catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
