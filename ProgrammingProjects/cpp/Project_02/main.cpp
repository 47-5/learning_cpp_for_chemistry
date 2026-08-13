#include <chemistry/molecule.h>

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cerr
            << "Usage: project02 <geometry-file> <hessian-file>\n";
        return 1;
    }

    try {
        const std::filesystem::path geometry_path{argv[1]};
        const std::filesystem::path hessian_path{argv[2]};

        const chemistry::Molecule molecule =
            chemistry::Molecule::from_file(geometry_path);

        std::ifstream hessian_file{hessian_path};
        if (!hessian_file) {
            throw std::runtime_error(
                "Cannot open Hessian file: " +
                hessian_path.string()
            );
        }

        std::cout
            << "Number of atoms: "
            << molecule.atom_count() << '\n'
            << "Hessian dimension: "
            << 3 * molecule.atom_count() << " x "
            << 3 * molecule.atom_count() << '\n';

        // Step 2: Read the Cartesian Hessian into an Eigen::MatrixXd.
    }
    catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}