#ifndef LEARNING_CPP_FOR_CHEMISTRY_MOLECULE_H
#define LEARNING_CPP_FOR_CHEMISTRY_MOLECULE_H

#include <Eigen/Core>

#include <cstddef>
#include <filesystem>
#include <vector>

namespace chemistry {

struct Atom {
    int atomic_number;
    double atomic_mass;
    Eigen::Vector3d position;
};

class Molecule {
public:
    explicit Molecule(std::vector<Atom> atoms);

    static Molecule from_file(const std::filesystem::path& path);

    [[nodiscard]] std::size_t atom_count() const noexcept;
    [[nodiscard]] const Atom& atom(std::size_t index) const;
    [[nodiscard]] const std::vector<Atom>& atoms() const noexcept;

private:
    std::vector<Atom> atoms_;
};

} // namespace chemistry

#endif // LEARNING_CPP_FOR_CHEMISTRY_MOLECULE_H
