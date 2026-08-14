#include "chemistry/vibrational_analysis.h"

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
        chemistry::VibrationalAnalysis vib = chemistry::VibrationalAnalysis{molecule, hessian_path};
        const Eigen::MatrixXd& hessian = vib.hessian();
        std::cout << "Hessian:\n" << hessian << '\n';

        // Step3: 获得质权坐标下的Hessian
        Eigen::MatrixXd F = vib.get_mass_weighted_hessian();
        std::cout <<"mass weighted Hessian:\n" << F << '\n';

        // Step4: 对角化F，从而得到特征值和normal mode
        const chemistry::NormalModeResult result = vib.get_normal_modes();

        std::cout
            << "Eigenvalues:\n"
            << result.eigenvalues
            << '\n';

        std::cout
            << "Normal modes:\n"
            << result.normal_modes
            << '\n';

        // Step5: 转换特征值的单位，得到以cm-1为单位的谐振频率
        const Eigen::VectorXd wavenumbers = vib.get_wavenumbers_cm_inverse();

        std::cout << "Wavenumbers (cm^-1):\n"
                  << wavenumbers << '\n';

    }
    catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}