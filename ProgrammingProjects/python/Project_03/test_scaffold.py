"""Smoke tests for the Project 3 data and class scaffold."""

from pathlib import Path
import unittest

import numpy as np

from chemistry.electronic_structure import (
    RHFCalculator,
    load_ao_integrals,
)


INPUT_DIRECTORY = (
    Path(__file__).resolve().parents[2]
    / "resources"
    / "Project_03"
    / "input"
    / "h2o"
    / "STO-3G"
)


class Project03ScaffoldTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.integrals = load_ao_integrals(INPUT_DIRECTORY)

    def test_loads_h2o_sto3g_dimensions(self) -> None:
        self.assertEqual(self.integrals.n_electrons, 10)
        self.assertEqual(self.integrals.n_basis, 7)
        self.assertEqual(self.integrals.overlap.shape, (7, 7))
        self.assertEqual(
            self.integrals.electron_repulsion.shape,
            (7, 7, 7, 7),
        )

    def test_restores_input_symmetry(self) -> None:
        np.testing.assert_allclose(
            self.integrals.overlap,
            self.integrals.overlap.T,
        )
        np.testing.assert_allclose(
            self.integrals.electron_repulsion,
            self.integrals.electron_repulsion.transpose(1, 0, 2, 3),
        )
        np.testing.assert_allclose(
            self.integrals.electron_repulsion,
            self.integrals.electron_repulsion.transpose(2, 3, 0, 1),
        )

    def test_constructs_closed_shell_calculator(self) -> None:
        calculator = RHFCalculator(self.integrals)
        self.assertEqual(calculator.n_occupied, 5)


if __name__ == "__main__":
    unittest.main()
