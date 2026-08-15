#include "minihf/integrals.h"
#include "minihf/rhf_calculator.h"

#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: project03 <Project-3-input-directory>\n";
        return 1;
    }

    try {
        const std::filesystem::path input_directory{argv[1]};
        const minihf::AOIntegrals integrals =
            minihf::AOIntegrals::from_directory(input_directory);
        const minihf::RHFCalculator calculator{integrals};

        std::cout
            << "Input directory: " << input_directory << '\n'
            << "Basis functions: " << integrals.basis_size() << '\n'
            << "Electrons: " << integrals.electron_count() << '\n'
            << "Occupied RHF orbitals: "
            << calculator.occupied_orbitals() << '\n'
            << std::fixed << std::setprecision(12)
            << "Nuclear repulsion: "
            << integrals.nuclear_repulsion() << " Eh\n"
            << "Integral loading is complete; RHF methods are TODOs.\n";

        // After implementing RHFCalculator::run():
        // const minihf::RHFResult result = calculator.run();
        // Print result.history here, just like the Python entry point.
    }
    catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
