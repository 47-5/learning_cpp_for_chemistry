#ifndef MINIHF_INTEGRALS_H
#define MINIHF_INTEGRALS_H

#include <Eigen/Dense>

#include <cstddef>
#include <filesystem>
#include <vector>

namespace minihf {

class AOIntegrals {
public:
    AOIntegrals(
        int electron_count,
        double nuclear_repulsion,
        Eigen::MatrixXd overlap,
        Eigen::MatrixXd kinetic,
        Eigen::MatrixXd nuclear_attraction,
        std::vector<double> electron_repulsion);

    static AOIntegrals from_directory(
        const std::filesystem::path& input_directory,
        int molecular_charge = 0);

    [[nodiscard]] int electron_count() const noexcept;
    [[nodiscard]] Eigen::Index basis_size() const noexcept;
    [[nodiscard]] double nuclear_repulsion() const noexcept;

    [[nodiscard]] const Eigen::MatrixXd& overlap() const noexcept;
    [[nodiscard]] const Eigen::MatrixXd& kinetic() const noexcept;
    [[nodiscard]] const Eigen::MatrixXd& nuclear_attraction() const noexcept;
    [[nodiscard]] Eigen::MatrixXd core_hamiltonian() const;

    [[nodiscard]] double eri(
        Eigen::Index mu,
        Eigen::Index nu,
        Eigen::Index lambda,
        Eigen::Index sigma) const;

private:
    [[nodiscard]] std::size_t eri_index(
        Eigen::Index mu,
        Eigen::Index nu,
        Eigen::Index lambda,
        Eigen::Index sigma) const;

    void validate() const;

    int electron_count_;
    double nuclear_repulsion_;
    Eigen::MatrixXd overlap_;
    Eigen::MatrixXd kinetic_;
    Eigen::MatrixXd nuclear_attraction_;
    std::vector<double> electron_repulsion_;
};

} // namespace minihf

#endif // MINIHF_INTEGRALS_H
