#include <chemistry/geometry_analyzer.h>

#include <chemistry/physical_constants.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace chemistry {

GeometryAnalyzer::GeometryAnalyzer(const Molecule& molecule)
    : molecule_(molecule)
{
}

double GeometryAnalyzer::distance(std::size_t i, std::size_t j) const
{
    return (
        molecule_.atom(i).position -
        molecule_.atom(j).position
    ).norm();
}

double GeometryAnalyzer::angle_degrees(
    std::size_t i,
    std::size_t j,
    std::size_t k
) const
{
    const Eigen::Vector3d ji =
        molecule_.atom(i).position -
        molecule_.atom(j).position;
    const Eigen::Vector3d jk =
        molecule_.atom(k).position -
        molecule_.atom(j).position;

    const double ji_norm = ji.norm();
    const double jk_norm = jk.norm();

    if (ji_norm == 0.0 || jk_norm == 0.0) {
        throw std::runtime_error(
            "Cannot calculate an angle with coincident atoms"
        );
    }

    const double cosine = std::clamp(
        ji.dot(jk) / (ji_norm * jk_norm),
        -1.0,
        1.0
    );

    return std::acos(cosine) * 180.0 / pi;
}

double GeometryAnalyzer::out_of_plane_angle_degrees(
    std::size_t i,
    std::size_t j,
    std::size_t k,
    std::size_t l
) const
{
    const Eigen::Vector3d ki =
        molecule_.atom(i).position -
        molecule_.atom(k).position;
    const Eigen::Vector3d kj =
        molecule_.atom(j).position -
        molecule_.atom(k).position;
    const Eigen::Vector3d kl =
        molecule_.atom(l).position -
        molecule_.atom(k).position;

    const Eigen::Vector3d plane_normal = kj.cross(kl);
    const double ki_norm = ki.norm();
    const double normal_norm = plane_normal.norm();

    if (ki_norm == 0.0) {
        throw std::runtime_error(
            "Out-of-plane atom coincides with the central atom"
        );
    }

    if (normal_norm < 1e-12) {
        return 0.0;
    }

    const double sine = std::clamp(
        ki.dot(plane_normal) / (ki_norm * normal_norm),
        -1.0,
        1.0
    );

    return std::asin(sine) * 180.0 / pi;
}

double GeometryAnalyzer::torsional_angle_degrees(
    std::size_t i,
    std::size_t j,
    std::size_t k,
    std::size_t l
) const
{
    const Eigen::Vector3d ij =
        molecule_.atom(j).position -
        molecule_.atom(i).position;
    const Eigen::Vector3d jk =
        molecule_.atom(k).position -
        molecule_.atom(j).position;
    const Eigen::Vector3d kl =
        molecule_.atom(l).position -
        molecule_.atom(k).position;

    const Eigen::Vector3d ijk_normal = ij.cross(jk);
    const Eigen::Vector3d jkl_normal = jk.cross(kl);

    const double ijk_norm = ijk_normal.norm();
    const double jkl_norm = jkl_normal.norm();

    if (ijk_norm < 1e-12 || jkl_norm < 1e-12) {
        throw std::runtime_error(
            "Torsional angle is undefined for collinear atoms"
        );
    }

    const double cosine = std::clamp(
        ijk_normal.dot(jkl_normal) / (ijk_norm * jkl_norm),
        -1.0,
        1.0
    );

    double angle = std::acos(cosine) * 180.0 / pi;
    const double sign =
        jk.dot(ijk_normal.cross(jkl_normal));

    if (sign < 0.0) {
        angle *= -1.0;
    }

    return angle;
}

} // namespace chemistry
