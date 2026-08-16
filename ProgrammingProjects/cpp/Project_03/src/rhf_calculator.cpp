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
    bool converged = false;

    Eigen::MatrixXd X = build_orthogonalizer();
    Eigen::MatrixXd F = integrals_.core_hamiltonian();
    RoothaanSolution roothaan_solution = solve_roothaan(F, X);
    Eigen::MatrixXd D = build_density(roothaan_solution.coefficients);

    F = build_fock(D);
    double EE = electronic_energy(D, F);
    double E = EE + integrals_.nuclear_repulsion();
    std::vector<SCFIteration> history;
    history.push_back(SCFIteration{0, EE, E, 0.0, 0.0, D});

    for (int run_index = 1; run_index < options_.max_iterations + 1; ++run_index ) {

        RoothaanSolution new_roothaan_solution = solve_roothaan(F, X);
        Eigen::MatrixXd new_D = build_density(new_roothaan_solution.coefficients);
        Eigen::MatrixXd new_F = build_fock(new_D);
        double new_EE = electronic_energy(new_D, new_F);
        double new_E = new_EE + integrals_.nuclear_repulsion();

        const double E_change = new_E - E;
        const double D_change = density_change(new_D, D);
        history.push_back(SCFIteration{run_index, new_EE, new_E,
            E_change, D_change, new_D});

        // 将新一轮结果提升为当前状态
        roothaan_solution = std::move(new_roothaan_solution);
        D = std::move(new_D);
        F = std::move(new_F);
        EE = std::move(new_EE);
        E = std::move(new_E);

        // 判断收敛
        if ((abs(E_change) < options_.energy_tolerance) && (abs(D_change) < options_.density_tolerance)) {
            converged = true;
            break;
        }
    }
    return RHFResult{converged, EE, E,
        roothaan_solution.orbital_energies, roothaan_solution.coefficients,
        D, F, history};
}

Eigen::MatrixXd RHFCalculator::build_orthogonalizer() const
{
    const Eigen::MatrixXd S = integrals_.overlap();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(S);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error(
            "Failed to diagonalize overlap"
        );
    }
    const Eigen::MatrixXd U = solver.eigenvectors();
    const Eigen::VectorXd eigenvalues = solver.eigenvalues();

    if (
    (eigenvalues.array() <= options_.overlap_eigenvalue_tolerance).any()) {
        throw std::runtime_error(
            "Overlap matrix is singular or "
            "nearly linearly dependent"
        );
    }

    const Eigen::VectorXd inverse_sqrt_eigenvalues = eigenvalues.array().sqrt().inverse();
    const Eigen::MatrixXd inverse_sqrt_s = inverse_sqrt_eigenvalues.asDiagonal();

    const Eigen::MatrixXd X = U * inverse_sqrt_s * U.transpose();
    return X;
}

RoothaanSolution RHFCalculator::solve_roothaan(
    const Eigen::MatrixXd& fock,
    const Eigen::MatrixXd& orthogonalizer) const
{
    Eigen::MatrixXd transformed_F = orthogonalizer.transpose() * fock * orthogonalizer;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(transformed_F);

    Eigen::MatrixXd C = orthogonalizer * solver.eigenvectors();
    return RoothaanSolution{solver.eigenvalues(), C};
}

Eigen::MatrixXd RHFCalculator::build_density(
    const Eigen::MatrixXd& coefficients) const
{
    const Eigen::MatrixXd occupied_coefficients = coefficients.leftCols(occupied_orbitals_);
    const Eigen::MatrixXd D = occupied_coefficients * occupied_coefficients.transpose();
    return D;
}

Eigen::MatrixXd RHFCalculator::build_fock(
    const Eigen::MatrixXd& density) const
{
    const Eigen::Index basis_size = integrals_.basis_size();
    if (density.rows() != basis_size || density.cols() != basis_size) {
        throw std::invalid_argument("Density matrix has the wrong shape");
    }

    Eigen::MatrixXd F = integrals_.core_hamiltonian();
    for (Eigen::Index mu = 0; mu < basis_size; ++mu) {
        for (Eigen::Index nu = 0; nu < basis_size; ++nu) {
            double two_electron_contribution = 0.0;

            for (Eigen::Index lambda = 0; lambda < basis_size; ++lambda) {
                for (Eigen::Index sigma = 0; sigma < basis_size; ++sigma) {
                    const double coulomb = integrals_.eri(mu, nu, lambda, sigma);
                    const double exchange = integrals_.eri(mu, lambda, nu, sigma);
                    two_electron_contribution += density(lambda, sigma) * (2 * coulomb - exchange);
                }
            }
            F(mu, nu) += two_electron_contribution;
        }
    }
    return F;
}

double RHFCalculator::electronic_energy(
    const Eigen::MatrixXd& density,
    const Eigen::MatrixXd& fock) const
{
    //     "evaluate sum D * (H_core + F)");
    const double ee = density.cwiseProduct(integrals_.core_hamiltonian() + fock).sum();
    return ee;
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
