import numpy as np
from pathlib import Path

from .molecule import Molecule
from .physical_constants import (
    AMU_TO_KILOGRAM,
    BOHR_TO_METER,
    HARTREE_TO_JOULE,
    SPEED_OF_LIGHT_CM_S,
)


class VibrationalAnalysis:
    def __init__(self, molecule: Molecule, hessian_path):
        self._molecule = molecule
        self._hessian = VibrationalAnalysis.read_hessian(hessian_path, self._molecule.atom_count)

    @staticmethod
    def read_hessian(hessian_path: str | Path, expected_n_atom: int):
        hessian_f = open(hessian_path).readlines()
        file_n_atom = int(hessian_f[0])
        if file_n_atom != expected_n_atom:
            raise ValueError(f"Expected {expected_n_atom} atoms, got {file_n_atom}")

        values = []
        for line in hessian_f[1:]:
            values += [float(i) for i in line.strip().split()]
        values = np.array(values)

        hessian = np.reshape(values, (3 * expected_n_atom, 3 * expected_n_atom))
        return hessian

    @property
    def hessian(self):
        return self._hessian

    def get_mass_weighted_hessian(self):
        hessian = self._hessian
        molecule = self._molecule
        atoms = molecule.atoms
        sqrt_mass_matrix = np.zeros_like(hessian)
        for i in range(3 * molecule.atom_count):
            for j in range(3 * molecule.atom_count):
                sqrt_mass_matrix[i, j] = np.sqrt(atoms[i // 3].atomic_mass * atoms[j // 3].atomic_mass)

        mass_weighted_hessian = hessian / sqrt_mass_matrix
        return mass_weighted_hessian

    @property
    def mass_weighted_hessian(self):
        return self.get_mass_weighted_hessian()

    def get_normal_modes(self):
        mass_weighted_hessian = self.get_mass_weighted_hessian()

        eigenvalues, normal_modes = np.linalg.eigh(
            mass_weighted_hessian
        )

        return eigenvalues, normal_modes
    def get_wavenumbers_cm_inverse(
        self,
        zero_tolerance: float = 1e-10,
    ) -> np.ndarray:
        """Return signed normal-mode wavenumbers in cm^-1.

        The eigenvalues have units hartree / (amu * bohr^2).
        Converting that unit to SI gives angular-frequency squared:

            omega^2 = eigenvalue * E_h / (m_u * a_0^2)

        The spectroscopic wavenumber is omega / (2*pi*c), with c in
        cm/s. Significant negative eigenvalues are returned as negative
        wavenumbers, the usual notation for imaginary modes.
        """
        if zero_tolerance < 0.0:
            raise ValueError(
                "zero_tolerance must be non-negative"
            )

        eigenvalues, _ = self.get_normal_modes()

        conversion_factor = (
            np.sqrt(
                HARTREE_TO_JOULE
                / (
                    AMU_TO_KILOGRAM
                    * BOHR_TO_METER**2
                )
            )
            / (2.0 * np.pi * SPEED_OF_LIGHT_CM_S)
        )

        wavenumbers = (
            np.sign(eigenvalues)
            * np.sqrt(np.abs(eigenvalues))
            * conversion_factor
        )

        numerical_zero = np.abs(eigenvalues) < zero_tolerance
        wavenumbers[numerical_zero] = 0.0

        return wavenumbers
