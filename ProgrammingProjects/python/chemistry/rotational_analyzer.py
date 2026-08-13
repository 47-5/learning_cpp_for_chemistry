import numpy as np

from .molecule import Molecule
from .physical_constants import (
    AMU_TO_GRAM,
    BOHR_TO_ANGSTROM,
    BOHR_TO_CM,
    PLANCK_CONSTANT_ERG_S,
    SPEED_OF_LIGHT_CM_S,
)


class RotationalAnalyzer:
    def __init__(self, molecule: Molecule) -> None:
        self._molecule = molecule

    def center_of_mass(self) -> np.ndarray:
        if self._molecule.atom_count == 0:
            raise ValueError(
                "Cannot calculate the center of mass "
                "of an empty molecule"
            )

        masses = self._molecule.masses.reshape(-1, 1)
        return np.sum(
            masses * self._molecule.positions,
            axis=0,
        ) / np.sum(masses)

    def inertia_tensor(self) -> np.ndarray:
        centered_positions = (
            self._molecule.positions
            - self.center_of_mass()
        )
        inertia = np.zeros((3, 3))

        for mass, position in zip(
            self._molecule.masses,
            centered_positions,
        ):
            x, y, z = position

            inertia[0, 0] += mass * (y**2 + z**2)
            inertia[1, 1] += mass * (x**2 + z**2)
            inertia[2, 2] += mass * (x**2 + y**2)

            inertia[0, 1] -= mass * x * y
            inertia[0, 2] -= mass * x * z
            inertia[1, 2] -= mass * y * z

        inertia[1, 0] = inertia[0, 1]
        inertia[2, 0] = inertia[0, 2]
        inertia[2, 1] = inertia[1, 2]

        return inertia

    def principal_moments(self) -> np.ndarray:
        return np.linalg.eigvalsh(self.inertia_tensor())

    def rotor_type(self) -> str:
        tolerance = 1e-4

        if self._molecule.atom_count == 2:
            return "diatomic"

        moment_a, moment_b, moment_c = (
            self.principal_moments()
        )

        if abs(moment_a) < tolerance:
            return "linear polyatomic"

        ab_equal = abs(moment_a - moment_b) < tolerance
        bc_equal = abs(moment_b - moment_c) < tolerance

        if ab_equal and bc_equal:
            return "spherical top"

        if ab_equal:
            return "oblate symmetric top"

        if bc_equal:
            return "prolate symmetric top"

        return "asymmetric top"

    def principal_moments_amu_angstrom2(
        self,
    ) -> np.ndarray:
        return (
            self.principal_moments()
            * BOHR_TO_ANGSTROM**2
        )

    def principal_moments_g_cm2(self) -> np.ndarray:
        return (
            self.principal_moments()
            * AMU_TO_GRAM
            * BOHR_TO_CM**2
        )

    def rotational_constants_mhz(self) -> np.ndarray:
        moments = self.principal_moments_g_cm2()

        if np.any(moments <= 0.0):
            raise ValueError(
                "Cannot calculate a rotational constant "
                "from a non-positive moment of inertia"
            )

        constants_hz = PLANCK_CONSTANT_ERG_S / (
            8.0 * np.pi**2 * moments
        )
        return constants_hz / 1e6

    def rotational_constants_cm_inverse(
        self,
    ) -> np.ndarray:
        return (
            self.rotational_constants_mhz()
            * 1e6
            / SPEED_OF_LIGHT_CM_S
        )
