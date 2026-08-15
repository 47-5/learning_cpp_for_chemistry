"""Configuration and result objects for SCF calculations."""

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class SCFOptions:
    energy_tolerance: float = 1.0e-10
    density_tolerance: float = 1.0e-8
    max_iterations: int = 100
    overlap_eigenvalue_tolerance: float = 1.0e-10

    def __post_init__(self) -> None:
        if self.energy_tolerance <= 0.0:
            raise ValueError("energy_tolerance must be positive")
        if self.density_tolerance <= 0.0:
            raise ValueError("density_tolerance must be positive")
        if self.max_iterations <= 0:
            raise ValueError("max_iterations must be positive")
        if self.overlap_eigenvalue_tolerance <= 0.0:
            raise ValueError(
                "overlap_eigenvalue_tolerance must be positive"
            )


@dataclass(frozen=True)
class SCFIteration:
    index: int
    electronic_energy: float
    total_energy: float
    energy_change: float
    density_change: float
    density: np.ndarray


@dataclass
class RHFResult:
    converged: bool
    electronic_energy: float
    total_energy: float
    orbital_energies: np.ndarray
    coefficients: np.ndarray
    density: np.ndarray
    fock: np.ndarray
    history: list[SCFIteration] = field(default_factory=list)

    @property
    def iterations(self) -> int:
        return len(self.history)
