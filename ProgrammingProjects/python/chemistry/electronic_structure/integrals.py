"""AO integral data and readers for the Project 3 teaching inputs."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class AOIntegrals:
    """All fixed input data needed by an electronic-structure method.

    The electron-repulsion integrals use chemist notation and dense
    storage::

        electron_repulsion[mu, nu, lam, sig]
            == (mu nu | lam sig)

    Dense storage is intentionally used in the Python teaching version.
    Project 3 is small enough that the clarity is worth more than ERI
    compression at this stage.
    """

    n_electrons: int
    nuclear_repulsion: float
    overlap: np.ndarray
    kinetic: np.ndarray
    nuclear_attraction: np.ndarray
    electron_repulsion: np.ndarray

    def __post_init__(self) -> None:
        self.overlap = np.asarray(self.overlap, dtype=float)
        self.kinetic = np.asarray(self.kinetic, dtype=float)
        self.nuclear_attraction = np.asarray(
            self.nuclear_attraction,
            dtype=float,
        )
        self.electron_repulsion = np.asarray(
            self.electron_repulsion,
            dtype=float,
        )
        self._validate()

    @property
    def n_basis(self) -> int:
        return self.overlap.shape[0]

    @property
    def core_hamiltonian(self) -> np.ndarray:
        """Return the one-electron core Hamiltonian H_core = T + V."""
        return self.kinetic + self.nuclear_attraction

    def _validate(self) -> None:
        if self.n_electrons <= 0:
            raise ValueError("n_electrons must be positive")

        if self.overlap.ndim != 2:
            raise ValueError("overlap must be a matrix")

        n_basis, n_columns = self.overlap.shape
        if n_basis != n_columns:
            raise ValueError("overlap must be square")

        matrix_shape = (n_basis, n_basis)
        for name, matrix in (
            ("kinetic", self.kinetic),
            ("nuclear_attraction", self.nuclear_attraction),
        ):
            if matrix.shape != matrix_shape:
                raise ValueError(
                    f"{name} has shape {matrix.shape}; "
                    f"expected {matrix_shape}"
                )

        eri_shape = (n_basis, n_basis, n_basis, n_basis)
        if self.electron_repulsion.shape != eri_shape:
            raise ValueError(
                "electron_repulsion has shape "
                f"{self.electron_repulsion.shape}; "
                f"expected {eri_shape}"
            )

        for name, matrix in (
            ("overlap", self.overlap),
            ("kinetic", self.kinetic),
            ("nuclear_attraction", self.nuclear_attraction),
        ):
            if not np.allclose(matrix, matrix.T):
                raise ValueError(f"{name} must be symmetric")


def load_ao_integrals(
    input_directory: str | Path,
    molecular_charge: int = 0,
) -> AOIntegrals:
    """Load the Project 3 integral files from one input directory."""
    directory = Path(input_directory)
    if not directory.is_dir():
        raise ValueError(f"Input directory does not exist: {directory}")

    overlap = _read_symmetric_matrix(directory / "s.dat")
    n_basis = overlap.shape[0]

    kinetic = _read_symmetric_matrix(
        directory / "t.dat",
        expected_size=n_basis,
    )
    nuclear_attraction = _read_symmetric_matrix(
        directory / "v.dat",
        expected_size=n_basis,
    )
    electron_repulsion = _read_eri(
        directory / "eri.dat",
        n_basis=n_basis,
    )
    nuclear_repulsion = _read_scalar(directory / "enuc.dat")
    n_electrons = (
        _read_total_nuclear_charge(directory / "geom.dat")
        - molecular_charge
    )

    return AOIntegrals(
        n_electrons=n_electrons,
        nuclear_repulsion=nuclear_repulsion,
        overlap=overlap,
        kinetic=kinetic,
        nuclear_attraction=nuclear_attraction,
        electron_repulsion=electron_repulsion,
    )


def _read_symmetric_matrix(
    path: Path,
    expected_size: int | None = None,
) -> np.ndarray:
    rows = _read_numeric_rows(path, expected_columns=3)
    inferred_size = max(max(int(row[0]), int(row[1])) for row in rows)

    if expected_size is not None and inferred_size > expected_size:
        raise ValueError(
            f"Index in {path} exceeds basis size {expected_size}"
        )

    size = expected_size or inferred_size
    matrix = np.zeros((size, size), dtype=float)

    for raw_i, raw_j, value in rows:
        i = int(raw_i) - 1
        j = int(raw_j) - 1
        _check_index(i, size, path)
        _check_index(j, size, path)
        matrix[i, j] = value
        matrix[j, i] = value

    return matrix


def _read_eri(path: Path, n_basis: int) -> np.ndarray:
    rows = _read_numeric_rows(path, expected_columns=5)
    eri = np.zeros((n_basis,) * 4, dtype=float)

    for raw_i, raw_j, raw_k, raw_l, value in rows:
        i, j, k, l = (
            int(raw_i) - 1,
            int(raw_j) - 1,
            int(raw_k) - 1,
            int(raw_l) - 1,
        )
        for index in (i, j, k, l):
            _check_index(index, n_basis, path)

        permutations = (
            (i, j, k, l),
            (j, i, k, l),
            (i, j, l, k),
            (j, i, l, k),
            (k, l, i, j),
            (l, k, i, j),
            (k, l, j, i),
            (l, k, j, i),
        )
        for permutation in permutations:
            eri[permutation] = value

    return eri


def _read_scalar(path: Path) -> float:
    text = _require_file(path).read_text(encoding="utf-8").strip()
    try:
        return float(text)
    except ValueError as error:
        raise ValueError(f"Expected one number in {path}") from error


def _read_total_nuclear_charge(path: Path) -> int:
    lines = _require_file(path).read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"Geometry file is empty: {path}")

    try:
        n_atoms = int(lines[0].strip())
    except ValueError as error:
        raise ValueError(
            f"First line of {path} must be the atom count"
        ) from error

    atom_lines = [line for line in lines[1:] if line.strip()]
    if len(atom_lines) != n_atoms:
        raise ValueError(
            f"{path} declares {n_atoms} atoms but contains "
            f"{len(atom_lines)} atom rows"
        )

    total_charge = 0
    for line_number, line in enumerate(atom_lines, start=2):
        fields = line.split()
        if not fields:
            continue
        try:
            nuclear_charge = float(fields[0])
        except ValueError as error:
            raise ValueError(
                f"Invalid nuclear charge in {path}:{line_number}"
            ) from error

        rounded_charge = round(nuclear_charge)
        if not np.isclose(nuclear_charge, rounded_charge):
            raise ValueError(
                f"Non-integral nuclear charge in {path}:{line_number}"
            )
        total_charge += int(rounded_charge)

    return total_charge


def _read_numeric_rows(
    path: Path,
    expected_columns: int,
) -> list[tuple[float, ...]]:
    rows: list[tuple[float, ...]] = []
    for line_number, line in enumerate(
        _require_file(path).read_text(encoding="utf-8").splitlines(),
        start=1,
    ):
        if not line.strip():
            continue

        fields = line.split()
        if len(fields) != expected_columns:
            raise ValueError(
                f"Expected {expected_columns} columns in "
                f"{path}:{line_number}; got {len(fields)}"
            )
        try:
            rows.append(tuple(float(field) for field in fields))
        except ValueError as error:
            raise ValueError(
                f"Invalid number in {path}:{line_number}"
            ) from error

    if not rows:
        raise ValueError(f"No data found in {path}")
    return rows


def _require_file(path: Path) -> Path:
    if not path.is_file():
        raise ValueError(f"Required input file does not exist: {path}")
    return path


def _check_index(index: int, size: int, path: Path) -> None:
    if not 0 <= index < size:
        raise ValueError(
            f"AO index {index + 1} in {path} is outside 1..{size}"
        )
