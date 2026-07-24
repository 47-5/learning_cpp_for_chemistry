#include <Eigen/Dense>

#include <iostream>

int main() {
    Eigen::Matrix2d fock;
    fock << 1.0, 0.2, 0.2, 2.0;

    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(fock);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Eigenvalue decomposition failed.\n";
        return 1;
    }

    std::cout << "Toy Fock matrix:\n" << fock << "\n\n";
    std::cout << "Orbital energies:\n" << solver.eigenvalues() << "\n\n";
    std::cout << "Orbital coefficients (columns):\n" << solver.eigenvectors() << '\n';

    return 0;
}
