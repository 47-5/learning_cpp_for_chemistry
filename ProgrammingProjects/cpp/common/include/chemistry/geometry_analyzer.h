#ifndef LEARNING_CPP_FOR_CHEMISTRY_GEOMETRY_ANALYZER_H
#define LEARNING_CPP_FOR_CHEMISTRY_GEOMETRY_ANALYZER_H

#include <chemistry/molecule.h>

#include <cstddef>

namespace chemistry {

class GeometryAnalyzer {
public:
    explicit GeometryAnalyzer(const Molecule& molecule);

    [[nodiscard]] double distance(std::size_t i, std::size_t j) const;

    [[nodiscard]] double angle_degrees(
        std::size_t i,
        std::size_t j,
        std::size_t k
    ) const;

    [[nodiscard]] double out_of_plane_angle_degrees(
        std::size_t i,
        std::size_t j,
        std::size_t k,
        std::size_t l
    ) const;

    [[nodiscard]] double torsional_angle_degrees(
        std::size_t i,
        std::size_t j,
        std::size_t k,
        std::size_t l
    ) const;

private:
    const Molecule& molecule_;
};

} // namespace chemistry

#endif // LEARNING_CPP_FOR_CHEMISTRY_GEOMETRY_ANALYZER_H
