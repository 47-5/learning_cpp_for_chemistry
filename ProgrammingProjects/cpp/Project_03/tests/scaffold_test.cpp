#include "minihf/integrals.h"
#include "minihf/rhf_calculator.h"

#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: project03_scaffold_test <input-directory>\n";
        return 1;
    }

    try {
        const minihf::AOIntegrals integrals =
            minihf::AOIntegrals::from_directory(
                std::filesystem::path{argv[1]});

        require(
            integrals.electron_count() == 10,
            "H2O should contain 10 electrons");
        require(
            integrals.basis_size() == 7,
            "H2O/STO-3G should contain 7 basis functions");
        require(
            std::abs(integrals.nuclear_repulsion()
                     - 8.002367061810450)
                < 1.0e-12,
            "Unexpected nuclear repulsion energy");
        require(
            integrals.overlap().isApprox(
                integrals.overlap().transpose()),
            "Overlap matrix was not restored symmetrically");
        require(
            std::abs(integrals.eri(1, 0, 0, 0)
                     - integrals.eri(0, 1, 0, 0))
                < 1.0e-14,
            "ERI permutation symmetry was not restored");

        const minihf::RHFCalculator calculator{integrals};
        require(
            calculator.occupied_orbitals() == 5,
            "H2O RHF should contain 5 occupied spatial orbitals");
    }
    catch (const std::exception& error) {
        std::cerr << "Test failure: " << error.what() << '\n';
        return 1;
    }

    std::cout << "Project 3 C++ scaffold checks passed.\n";
    return 0;
}
