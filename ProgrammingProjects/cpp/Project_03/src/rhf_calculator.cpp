#include "minihf/rhf_calculator.h"

#include <stdexcept>
#include <utility>

namespace minihf {

void SCFOptions::validate() const
{
    if (energy_tolerance <= 0.0) {
        throw std::invalid_argument(
            "energy_tolerance must be positive");
    }
    if (density_tolerance <= 0.0) {
        throw std::invalid_argument(
            "density_tolerance must be positive");
    }
    if (max_iterations <= 0) {
        throw std::invalid_argument(
            "max_iterations must be positive");
    }
    if (overlap_eigenvalue_tolerance <= 0.0) {
        throw std::invalid_argument(
            "overlap_eigenvalue_tolerance must be positive");
    }
}

int RHFResult::iterations() const noexcept
{
    return history.empty() ? 0 : history.back().index;
}

RHFCalculator::RHFCalculator(
    const AOIntegrals& integrals,
    SCFOptions options)
    : integrals_{integrals},
      options_{std::move(options)},
      occupied_orbitals_{integrals.electron_count() / 2}
{
    options_.validate();

    if (integrals.electron_count() % 2 != 0) {
        throw std::invalid_argument(
            "RHFCalculator requires an even electron count");
    }
    if (occupied_orbitals_ > integrals.basis_size()) {
        throw std::invalid_argument(
            "Occupied orbitals do not fit in the AO basis");
    }
}

Eigen::Index RHFCalculator::occupied_orbitals() const noexcept
{
    return occupied_orbitals_;
}

RHFResult RHFCalculator::run() const
{
    throw std::logic_error(
        "TODO 6: assemble the RHF SCF loop");
}

Eigen::MatrixXd RHFCalculator::build_orthogonalizer() const
{
    throw std::logic_error(
        "TODO 1: diagonalize S and construct S^(-1/2)");
}

RoothaanSolution RHFCalculator::solve_roothaan(
    const Eigen::MatrixXd& /*fock*/,
    const Eigen::MatrixXd& /*orthogonalizer*/) const
{
    throw std::logic_error(
        "TODO 2: solve the transformed Fock eigenproblem");
}

Eigen::MatrixXd RHFCalculator::build_density(
    const Eigen::MatrixXd& /*coefficients*/) const
{
    throw std::logic_error(
        "TODO 3: form D from occupied MO columns without factor 2");
}

Eigen::MatrixXd RHFCalculator::build_fock(
    const Eigen::MatrixXd& /*density*/) const
{
    throw std::logic_error(
        "TODO 4: contract D with Coulomb and exchange ERIs");
}

double RHFCalculator::electronic_energy(
    const Eigen::MatrixXd& /*density*/,
    const Eigen::MatrixXd& /*fock*/) const
{
    throw std::logic_error(
        "TODO 5: evaluate sum D * (H_core + F)");
}

double RHFCalculator::density_change(
    const Eigen::MatrixXd& old_density,
    const Eigen::MatrixXd& new_density)
{
    if (old_density.rows() != new_density.rows()
        || old_density.cols() != new_density.cols()) {
        throw std::invalid_argument(
            "Density matrices must have the same shape");
    }
    return (new_density - old_density).norm();
}

} // namespace minihf
