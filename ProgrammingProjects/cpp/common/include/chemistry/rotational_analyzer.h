#ifndef LEARNING_CPP_FOR_CHEMISTRY_ROTATIONAL_ANALYZER_H
#define LEARNING_CPP_FOR_CHEMISTRY_ROTATIONAL_ANALYZER_H

#include <chemistry/molecule.h>

#include <Eigen/Core>

#include <string>

namespace chemistry {

class RotationalAnalyzer {
public:
    explicit RotationalAnalyzer(const Molecule& molecule);

    [[nodiscard]] Eigen::Vector3d center_of_mass() const;
    [[nodiscard]] Eigen::Matrix3d inertia_tensor() const;
    [[nodiscard]] Eigen::Vector3d principal_moments() const;
    [[nodiscard]] std::string rotor_type() const;

    [[nodiscard]] Eigen::Vector3d
    principal_moments_amu_angstrom2() const;

    [[nodiscard]] Eigen::Vector3d
    principal_moments_g_cm2() const;

    [[nodiscard]] Eigen::Vector3d
    rotational_constants_mhz() const;

    [[nodiscard]] Eigen::Vector3d
    rotational_constants_cm_inverse() const;

private:
    const Molecule& molecule_;
};

} // namespace chemistry

#endif // LEARNING_CPP_FOR_CHEMISTRY_ROTATIONAL_ANALYZER_H
