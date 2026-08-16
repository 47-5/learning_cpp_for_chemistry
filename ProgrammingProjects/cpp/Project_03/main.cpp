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
            << "Input directory: " << input_directory.string() << '\n'
            << "Basis functions: " << integrals.basis_size() << '\n'
            << "Electrons: " << integrals.electron_count() << '\n'
            << "Occupied RHF orbitals: "
            << calculator.occupied_orbitals() << '\n'
            << std::fixed << std::setprecision(12)
            << "Nuclear repulsion: "
            << integrals.nuclear_repulsion() << " Eh\n";

        const minihf::RHFResult result = calculator.run();

        std::cout
            << "\nSCF iterations\n"
            << std::setw(4) << "Iter"
            << ' ' << std::setw(20) << "E(elec) / Eh"
            << ' ' << std::setw(20) << "E(total) / Eh"
            << ' ' << std::setw(14) << "Delta E / Eh"
            << ' ' << std::setw(14) << "Delta D" << '\n'
            << "------------------------------------------------------------------------------\n";

        for (const minihf::SCFIteration& iteration : result.history) {
            std::cout
                << std::setw(4) << iteration.index
                << ' ' << std::fixed << std::setprecision(12)
                << std::setw(20) << iteration.electronic_energy
                << ' ' << std::setw(20) << iteration.total_energy
                << ' ' << std::scientific << std::setprecision(6)
                << std::setw(14) << iteration.energy_change
                << ' ' << std::setw(14) << iteration.density_change
                << '\n';
        }

        std::cout
            << "\nConverged: "
            << (result.converged ? "True" : "False") << '\n'
            << "SCF iterations: " << result.iterations() << '\n'
            << std::fixed << std::setprecision(12)
            << "Final electronic energy: "
            << result.electronic_energy << " Eh\n"
            << "Nuclear repulsion: "
            << integrals.nuclear_repulsion() << " Eh\n"
            << "Final total energy: "
            << result.total_energy << " Eh\n";

    }
    catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
