//
// Created by lrc on 2026/8/14.
//

#ifndef LEARNING_CPP_FOR_CHEMISTRY_VIBRATIONAL_ANALYSIS_H
#define LEARNING_CPP_FOR_CHEMISTRY_VIBRATIONAL_ANALYSIS_H
#include <chemistry/molecule.h>
#include <Eigen/Dense>
#include <cstddef>
#include <filesystem>
#include <vector>

namespace chemistry {

struct NormalModeResult {
    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd normal_modes;
};

class VibrationalAnalysis {
    public:
    explicit VibrationalAnalysis(Molecule molecule, const std::filesystem::path& hessian_path);
    [[nodiscard]] const Eigen::MatrixXd& hessian() const;
    [[nodiscard]] Eigen::MatrixXd get_mass_weighted_hessian() const;

    [[nodiscard]] NormalModeResult get_normal_modes() const;

    [[nodiscard]] Eigen::VectorXd get_wavenumbers_cm_inverse(
        double zero_tolerance = 1e-10) const;

    private:
    static Eigen::MatrixXd read_hessian(
        const std::filesystem::path& path,
        std::size_t expected_atom_count);

    Molecule _molecule;
    Eigen::MatrixXd _hessian;

};

}

#endif // LEARNING_CPP_FOR_CHEMISTRY_VIBRATIONAL_ANALYSIS_H
