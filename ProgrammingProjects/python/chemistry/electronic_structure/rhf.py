"""Restricted Hartree--Fock scaffold for Project 3."""

import numpy as np

from .integrals import AOIntegrals
from .models import RHFResult, SCFOptions, SCFIteration


class RHFCalculator:
    """Closed-shell RHF calculator using the Project 3 density convention.

    The density matrix used here does not contain the occupation factor 2::

        D[mu, nu] = sum_i_occ C[mu, i] C[nu, i]

    Complete the private methods in their numerical order. Keeping these
    operations separate makes each formula independently testable before
    assembling the SCF loop.
    """

    def __init__(
        self,
        integrals: AOIntegrals,
        options: SCFOptions | None = None,
    ) -> None:
        if integrals.n_electrons % 2 != 0:
            raise ValueError(
                "RHFCalculator requires an even number of electrons"
            )

        n_occupied = integrals.n_electrons // 2
        if n_occupied > integrals.n_basis:
            raise ValueError(
                f"{n_occupied} occupied orbitals do not fit in "
                f"{integrals.n_basis} basis functions"
            )

        self.integrals = integrals
        self.options = options or SCFOptions()
        self.n_occupied = n_occupied

    def run(self) -> RHFResult:
        """Run the SCF calculation and return all final numerical data.

        Suggested implementation order:

        1. Build X = S^(-1/2).
        2. Diagonalize the core Hamiltonian for an initial density.
        3. Repeatedly build F[D], solve the Roothaan--Hall equation,
           and build a new density.
        4. Check both energy and density convergence.
        5. Return an RHFResult instead of printing from this class.
        """
        converged = False

        X = self._build_orthogonalizer()
        F0 = self.integrals.core_hamiltonian
        orbital_energies, orbital_coefficients = self._solve_roothaan(fock=F0, orthogonalizer=X)
        D = self._build_density(coefficients=orbital_coefficients)
        F = self._build_fock(density=D)
        EE = self._electronic_energy(density=D, fock=F)
        E = EE + self.integrals.nuclear_repulsion
        history = [SCFIteration(index=0, electronic_energy=EE, total_energy=E, energy_change=0, density_change=0, density=D)]
        for run_index in range(1, self.options.max_iterations + 1):
            orbital_energies, orbital_coefficients = self._solve_roothaan(fock=F, orthogonalizer=X)
            D = self._build_density(coefficients=orbital_coefficients)
            F = self._build_fock(density=D)
            EE = self._electronic_energy(density=D, fock=F)
            E = EE + self.integrals.nuclear_repulsion
            energy_change = E - history[-1].total_energy
            density_change = self._density_change(history[-1].density, D)
            history.append(SCFIteration(index=run_index, electronic_energy=EE, total_energy=E,
                                        energy_change=energy_change, density_change=density_change, density=D))
            if (np.abs(energy_change) < self.options.energy_tolerance) and \
                (np.abs(density_change) < self.options.density_tolerance):
                converged = True
                break

        return RHFResult(
            converged=converged,
            electronic_energy=EE, total_energy=E,
            orbital_energies=orbital_energies, coefficients=orbital_coefficients,
            density=D, fock=F,
            history=history,
        )


    def _build_orthogonalizer(self) -> np.ndarray:
        """Return symmetric orthogonalizer X = S^(-1/2)."""
        S = self.integrals.overlap
        eigen_values, eigen_vectors = np.linalg.eigh(S)
        tolerance = self.options.overlap_eigenvalue_tolerance

        if np.any(eigen_values <= tolerance):
            raise ValueError(
                "Overlap matrix is singular or nearly linearly dependent; "
                f"smallest eigenvalue: {eigen_values.min():.3e}"
            )

        inverse_s_sqrt = np.diag(1 / np.sqrt(eigen_values))
        X = eigen_vectors @ inverse_s_sqrt @ eigen_vectors.T
        return X

    def _solve_roothaan(
        self,
        fock: np.ndarray,
        orthogonalizer: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return ascending orbital energies and AO coefficients.

        Transform F to X.T @ F @ X, diagonalize it, then transform the
        eigenvectors back with C = X @ C_orthogonal.
        """
        transformed_F = orthogonalizer.T @ fock @ orthogonalizer
        orbital_energies, transformed_orbital_coefficients = np.linalg.eigh(transformed_F)
        orbital_coefficients = orthogonalizer @ transformed_orbital_coefficients
        return orbital_energies, orbital_coefficients


    def _build_density(self, coefficients: np.ndarray) -> np.ndarray:
        """Return D = C_occ @ C_occ.T without an occupation factor 2."""
        occupied_coefficients = coefficients[:, :self.n_occupied]
        D = occupied_coefficients @ occupied_coefficients.T
        return D

    def _build_fock(self, density: np.ndarray) -> np.ndarray:
        """Return F = H_core + sum D * [2 Coulomb - exchange]."""
        eri = self.integrals.electron_repulsion
        coulomb = np.einsum(
            "ls,mnls->mn",
            density,
            eri,
            optimize=True,
        )
        exchange = np.einsum(
            "ls,mlns->mn",
            density,
            eri,
            optimize=True,
        )

        F = self.integrals.core_hamiltonian + 2 * coulomb - exchange
        return F

    def _electronic_energy(
        self,
        density: np.ndarray,
        fock: np.ndarray,
    ) -> float:
        """Return E_elec = sum D * (H_core + F)."""
        energy = np.einsum(
            "mn,mn->",
            density,
            self.integrals.core_hamiltonian + fock
        )
        return float(energy)

    @staticmethod
    def _density_change(
        old_density: np.ndarray,
        new_density: np.ndarray,
    ) -> float:
        """Return the Frobenius norm used by the Project 3 notes."""
        return float(np.linalg.norm(new_density - old_density))
