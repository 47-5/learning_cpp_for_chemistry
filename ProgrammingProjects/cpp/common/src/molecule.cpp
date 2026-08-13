#include <chemistry/molecule.h>

#include <chemistry/atomic_masses.h>

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace chemistry {

Molecule::Molecule(std::vector<Atom> atoms)
    : atoms_(std::move(atoms))
{
}

Molecule Molecule::from_file(const std::filesystem::path& path)
{
    std::ifstream file{path};
    if (!file) {
        throw std::runtime_error(
            "Cannot open geometry file: " + path.string()
        );
    }

    std::size_t atom_count{};
    if (!(file >> atom_count)) {
        throw std::runtime_error(
            "Failed to read atom count from: " + path.string()
        );
    }

    std::vector<Atom> atoms;
    atoms.reserve(atom_count);

    for (std::size_t index = 0; index < atom_count; ++index) {
        double atomic_number_value{};
        double x{};
        double y{};
        double z{};

        if (!(file >> atomic_number_value >> x >> y >> z)) {
            throw std::runtime_error(
                "Failed to read atom at index " +
                std::to_string(index)
            );
        }

        if (!std::isfinite(atomic_number_value) ||
            std::floor(atomic_number_value) != atomic_number_value) {
            throw std::runtime_error(
                "Invalid atomic number at index " +
                std::to_string(index)
            );
        }

        const auto atomic_number =
            static_cast<int>(atomic_number_value);

        if (atomic_number <= 0 ||
            static_cast<std::size_t>(atomic_number) >=
                atomic_masses.size()) {
            throw std::runtime_error(
                "Atomic number out of range at index " +
                std::to_string(index)
            );
        }

        atoms.push_back(Atom{
            atomic_number,
            atomic_masses[static_cast<std::size_t>(atomic_number)],
            Eigen::Vector3d{x, y, z},
        });
    }

    return Molecule{std::move(atoms)};
}

std::size_t Molecule::atom_count() const noexcept
{
    return atoms_.size();
}

const Atom& Molecule::atom(std::size_t index) const
{
    return atoms_.at(index);
}

const std::vector<Atom>& Molecule::atoms() const noexcept
{
    return atoms_;
}

} // namespace chemistry
