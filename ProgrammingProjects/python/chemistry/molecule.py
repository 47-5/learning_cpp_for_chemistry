from pathlib import Path

import numpy as np

from .atomic_masses import ATOMIC_MASSES


class Atom:
    def __init__(
        self,
        atomic_number: int,
        position: np.ndarray,
    ) -> None:
        if not 1 <= atomic_number < len(ATOMIC_MASSES):
            raise ValueError(
                f"Invalid atomic number: {atomic_number}"
            )

        self.atomic_number = atomic_number
        self.atomic_mass = ATOMIC_MASSES[atomic_number]
        self.position = np.asarray(position, dtype=float)


class Molecule:
    def __init__(self, atoms: list[Atom]) -> None:
        self._atoms = tuple(atoms)

    @classmethod
    def from_file(cls, path: str | Path) -> "Molecule":
        geometry_path = Path(path)

        try:
            lines = geometry_path.read_text(
                encoding="utf-8"
            ).splitlines()
        except OSError as error:
            raise ValueError(
                f"Cannot open geometry file: {geometry_path}"
            ) from error

        if not lines:
            raise ValueError(
                f"Geometry file is empty: {geometry_path}"
            )

        try:
            atom_count = int(lines[0])
        except ValueError as error:
            raise ValueError(
                "Failed to read atom count"
            ) from error

        if len(lines) - 1 < atom_count:
            raise ValueError(
                f"Expected {atom_count} atoms, "
                f"but found {len(lines) - 1}"
            )

        atoms: list[Atom] = []

        for index, line in enumerate(
            lines[1 : atom_count + 1]
        ):
            fields = line.split()
            if len(fields) != 4:
                raise ValueError(
                    f"Invalid atom record at index {index}"
                )

            try:
                atomic_number_value = float(fields[0])
                coordinates = np.array(
                    fields[1:],
                    dtype=float,
                )
            except ValueError as error:
                raise ValueError(
                    f"Invalid atom record at index {index}"
                ) from error

            if not atomic_number_value.is_integer():
                raise ValueError(
                    f"Invalid atomic number at index {index}"
                )

            atoms.append(
                Atom(
                    int(atomic_number_value),
                    coordinates,
                )
            )

        return cls(atoms)

    @property
    def atom_count(self) -> int:
        return len(self._atoms)

    @property
    def atoms(self) -> tuple[Atom, ...]:
        return self._atoms

    def atom(self, index: int) -> Atom:
        if not 0 <= index < self.atom_count:
            raise IndexError("Atom index out of range")

        return self._atoms[index]

    @property
    def positions(self) -> np.ndarray:
        return np.array(
            [atom.position for atom in self._atoms]
        )

    @property
    def atomic_numbers(self) -> np.ndarray:
        return np.array(
            [atom.atomic_number for atom in self._atoms],
            dtype=int,
        )

    @property
    def masses(self) -> np.ndarray:
        return np.array(
            [atom.atomic_mass for atom in self._atoms]
        )
