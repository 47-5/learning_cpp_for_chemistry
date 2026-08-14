#include "chemistry/vibrational_analysis.h"
#include "chemistry/physical_constants.h"

#include <filesystem>
#include <fstream>
#include <utility>
#include <cmath>
#include <Eigen/Eigenvalues>
#include <stdexcept>
//
// Created by lrc on 2026/8/14.
//
namespace chemistry {

VibrationalAnalysis::VibrationalAnalysis(Molecule molecule, const std::filesystem::path& hessian_path):
_molecule(std::move(molecule)), _hessian(read_hessian(hessian_path, _molecule.atom_count()))
{
}
const Eigen::MatrixXd& VibrationalAnalysis::hessian() const {
    return _hessian;
}

Eigen::MatrixXd VibrationalAnalysis::read_hessian(const std::filesystem::path& path,
                                                  std::size_t expected_atom_count) {
    std::ifstream hessian_f(path);

    std::size_t file_atom_count;
    if (!(hessian_f >> file_atom_count)) {
        throw std::runtime_error("Can not open hessian file");
    }

    if (file_atom_count != expected_atom_count) {
        throw std::runtime_error("Atom count in Hessian file does not match expected number of atoms");
    }

    const Eigen::Index hessian_dim = 3 * static_cast<Eigen::Index>(expected_atom_count);

    Eigen::MatrixXd hessian(hessian_dim, hessian_dim);

    for (Eigen::Index i = 0; i < hessian_dim; ++i) {
        for (Eigen::Index j = 0; j < hessian_dim; ++j) {
            if (!(hessian_f >> hessian(i, j))) {
                throw std::runtime_error("Not enough values in Hessian file");
            }
        }
    }
    return hessian;
}

Eigen::MatrixXd VibrationalAnalysis::get_mass_weighted_hessian() const {
    int hessian_dim = 3 * _molecule.atom_count();
    Eigen::MatrixXd sqrt_mass_matrix(hessian_dim, hessian_dim);

    for (int i = 0; i < hessian_dim; ++i) {
        for (int j = 0; j < hessian_dim; ++j) {
            sqrt_mass_matrix(i, j) = sqrt(_molecule.atoms()[i/3].atomic_mass * _molecule.atoms()[j/3].atomic_mass) ;
        }
    }

    Eigen::MatrixXd mass_weighted_hessian = _hessian.cwiseQuotient(sqrt_mass_matrix);
    return mass_weighted_hessian;
}

NormalModeResult VibrationalAnalysis::get_normal_modes() const {
    const Eigen::MatrixXd mass_weighted_hessian = get_mass_weighted_hessian();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver{mass_weighted_hessian};

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error(
            "Failed to diagonalize mass-weighted Hessian"
        );
    }

    return NormalModeResult{solver.eigenvalues(), solver.eigenvectors()};
}

Eigen::VectorXd VibrationalAnalysis::get_wavenumbers_cm_inverse(
    double zero_tolerance) const {
    if (zero_tolerance < 0.0) {
        throw std::invalid_argument(
            "zero_tolerance must be non-negative"
        );
    }

    const Eigen::VectorXd eigenvalues =
        get_normal_modes().eigenvalues;

    const double conversion_factor =
        std::sqrt(
            hartree_to_joule
            / (amu_to_kilogram
               * bohr_to_meter
               * bohr_to_meter)
        )
        / (2.0 * pi * speed_of_light_cm_s);

    Eigen::VectorXd wavenumbers(eigenvalues.size());

    for (Eigen::Index index = 0;
         index < eigenvalues.size();
         ++index) {
        const double eigenvalue = eigenvalues(index);

        if (std::abs(eigenvalue) < zero_tolerance) {
            wavenumbers(index) = 0.0;
        } else {
            wavenumbers(index) = std::copysign(
                std::sqrt(std::abs(eigenvalue))
                    * conversion_factor,
                eigenvalue
            );
        }
    }

    return wavenumbers;
}

} // namespace chemistry
