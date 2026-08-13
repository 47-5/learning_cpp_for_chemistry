#include <chemistry/rotational_analyzer.h>

#include <chemistry/physical_constants.h>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <stdexcept>

namespace chemistry {

RotationalAnalyzer::RotationalAnalyzer(const Molecule& molecule)
    : molecule_(molecule)
{
}

Eigen::Vector3d RotationalAnalyzer::center_of_mass() const
{
    if (molecule_.atom_count() == 0) {
        throw std::runtime_error(
            "Cannot calculate the center of mass of an empty molecule"
        );
    }

    Eigen::Vector3d weighted_sum = Eigen::Vector3d::Zero();
    double total_mass = 0.0;

    for (const Atom& atom : molecule_.atoms()) {
        weighted_sum += atom.atomic_mass * atom.position;
        total_mass += atom.atomic_mass;
    }

    return weighted_sum / total_mass;
}

Eigen::Matrix3d RotationalAnalyzer::inertia_tensor() const
{
    Eigen::Matrix3d inertia = Eigen::Matrix3d::Zero();
    const Eigen::Vector3d center = center_of_mass();

    for (const Atom& atom : molecule_.atoms()) {
        const Eigen::Vector3d position = atom.position - center;
        const double x = position.x();
        const double y = position.y();
        const double z = position.z();
        const double mass = atom.atomic_mass;

        inertia(0, 0) += mass * (y * y + z * z);
        inertia(1, 1) += mass * (x * x + z * z);
        inertia(2, 2) += mass * (x * x + y * y);

        inertia(0, 1) -= mass * x * y;
        inertia(0, 2) -= mass * x * z;
        inertia(1, 2) -= mass * y * z;
    }

    inertia(1, 0) = inertia(0, 1);
    inertia(2, 0) = inertia(0, 2);
    inertia(2, 1) = inertia(1, 2);

    return inertia;
}

Eigen::Vector3d RotationalAnalyzer::principal_moments() const
{
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver{
        inertia_tensor()
    };

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error(
            "Failed to diagonalize inertia tensor"
        );
    }

    return solver.eigenvalues();
}

std::string RotationalAnalyzer::rotor_type() const
{
    constexpr double tolerance = 1e-4;

    if (molecule_.atom_count() == 2) {
        return "diatomic";
    }

    const Eigen::Vector3d moments = principal_moments();

    if (std::abs(moments[0]) < tolerance) {
        return "linear polyatomic";
    }

    const bool ab_equal =
        std::abs(moments[0] - moments[1]) < tolerance;
    const bool bc_equal =
        std::abs(moments[1] - moments[2]) < tolerance;

    if (ab_equal && bc_equal) {
        return "spherical top";
    }

    if (ab_equal) {
        return "oblate symmetric top";
    }

    if (bc_equal) {
        return "prolate symmetric top";
    }

    return "asymmetric top";
}

Eigen::Vector3d
RotationalAnalyzer::principal_moments_amu_angstrom2() const
{
    return principal_moments() *
        (bohr_to_angstrom * bohr_to_angstrom);
}

Eigen::Vector3d
RotationalAnalyzer::principal_moments_g_cm2() const
{
    return principal_moments() *
        amu_to_gram *
        (bohr_to_cm * bohr_to_cm);
}

Eigen::Vector3d
RotationalAnalyzer::rotational_constants_mhz() const
{
    const Eigen::Vector3d moments =
        principal_moments_g_cm2();
    Eigen::Vector3d constants;

    for (Eigen::Index index = 0;
         index < moments.size();
         ++index) {
        if (moments[index] <= 0.0) {
            throw std::runtime_error(
                "Cannot calculate a rotational constant "
                "from a non-positive moment of inertia"
            );
        }

        const double constant_hz =
            planck_constant_erg_s /
            (8.0 * pi * pi * moments[index]);

        constants[index] = constant_hz / 1e6;
    }

    return constants;
}

Eigen::Vector3d
RotationalAnalyzer::rotational_constants_cm_inverse() const
{
    return rotational_constants_mhz() *
        1e6 /
        speed_of_light_cm_s;
}

} // namespace chemistry
