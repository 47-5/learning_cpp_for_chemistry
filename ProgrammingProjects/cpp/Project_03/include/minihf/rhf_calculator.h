#ifndef MINIHF_RHF_CALCULATOR_H
#define MINIHF_RHF_CALCULATOR_H

#include "minihf/integrals.h"

#include <Eigen/Dense>

#include <vector>

namespace minihf {

struct SCFOptions {
    double energy_tolerance{1.0e-10};
    double density_tolerance{1.0e-8};
    int max_iterations{100};
    double overlap_eigenvalue_tolerance{1.0e-10};

    void validate() const;
};

struct SCFIteration {
    int index{};
    double electronic_energy{};
    double total_energy{};
    double energy_change{};
    double density_change{};
    Eigen::MatrixXd density;
};

struct RHFResult {
    bool converged{};
    double electronic_energy{};
    double total_energy{};
    Eigen::VectorXd orbital_energies;
    Eigen::MatrixXd coefficients;
    Eigen::MatrixXd density;
    Eigen::MatrixXd fock;
    std::vector<SCFIteration> history;

    [[nodiscard]] int iterations() const noexcept;
};

struct RoothaanSolution {
    Eigen::VectorXd orbital_energies;
    Eigen::MatrixXd coefficients;
};

class RHFCalculator {
public:
    explicit RHFCalculator(
        const AOIntegrals& integrals,
        SCFOptions options = {});

    [[nodiscard]] Eigen::Index occupied_orbitals() const noexcept;
    [[nodiscard]] RHFResult run() const;

    // These numerical building blocks are public so that each Project 3
    // step can be implemented and tested independently.
    [[nodiscard]] Eigen::MatrixXd build_orthogonalizer() const;

    [[nodiscard]] RoothaanSolution solve_roothaan(
        const Eigen::MatrixXd& fock,
        const Eigen::MatrixXd& orthogonalizer) const;

    [[nodiscard]] Eigen::MatrixXd build_density(
        const Eigen::MatrixXd& coefficients) const;

    [[nodiscard]] Eigen::MatrixXd build_fock(
        const Eigen::MatrixXd& density) const;

    [[nodiscard]] double electronic_energy(
        const Eigen::MatrixXd& density,
        const Eigen::MatrixXd& fock) const;

    [[nodiscard]] static double density_change(
        const Eigen::MatrixXd& old_density,
        const Eigen::MatrixXd& new_density);

private:
    const AOIntegrals& integrals_;
    SCFOptions options_;
    Eigen::Index occupied_orbitals_;
};

} // namespace minihf

#endif // MINIHF_RHF_CALCULATOR_H
