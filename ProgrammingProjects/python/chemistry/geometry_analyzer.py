import numpy as np

from .molecule import Molecule


class GeometryAnalyzer:
    def __init__(self, molecule: Molecule) -> None:
        self._molecule = molecule

    def distance(self, i: int, j: int) -> float:
        displacement = (
            self._molecule.atom(i).position
            - self._molecule.atom(j).position
        )
        return float(np.linalg.norm(displacement))

    def distances_matrix(self) -> np.ndarray:
        count = self._molecule.atom_count
        return np.array(
            [
                [
                    self.distance(i, j)
                    for j in range(count)
                ]
                for i in range(count)
            ]
        )

    def angle_degrees(
        self,
        i: int,
        j: int,
        k: int,
    ) -> float:
        ji = (
            self._molecule.atom(i).position
            - self._molecule.atom(j).position
        )
        jk = (
            self._molecule.atom(k).position
            - self._molecule.atom(j).position
        )

        ji_norm = np.linalg.norm(ji)
        jk_norm = np.linalg.norm(jk)

        if ji_norm == 0.0 or jk_norm == 0.0:
            raise ValueError(
                "Cannot calculate an angle "
                "with coincident atoms"
            )

        cosine = np.clip(
            np.dot(ji, jk) / (ji_norm * jk_norm),
            -1.0,
            1.0,
        )
        return float(np.degrees(np.arccos(cosine)))

    def out_of_plane_angle_degrees(
        self,
        i: int,
        j: int,
        k: int,
        l: int,
    ) -> float:
        ki = (
            self._molecule.atom(i).position
            - self._molecule.atom(k).position
        )
        kj = (
            self._molecule.atom(j).position
            - self._molecule.atom(k).position
        )
        kl = (
            self._molecule.atom(l).position
            - self._molecule.atom(k).position
        )

        plane_normal = np.cross(kj, kl)
        ki_norm = np.linalg.norm(ki)
        normal_norm = np.linalg.norm(plane_normal)

        if ki_norm == 0.0:
            raise ValueError(
                "Out-of-plane atom coincides "
                "with the central atom"
            )

        if normal_norm < 1e-12:
            return 0.0

        sine = np.clip(
            np.dot(ki, plane_normal)
            / (ki_norm * normal_norm),
            -1.0,
            1.0,
        )
        return float(np.degrees(np.arcsin(sine)))

    def torsional_angle_degrees(
        self,
        i: int,
        j: int,
        k: int,
        l: int,
    ) -> float:
        ij = (
            self._molecule.atom(j).position
            - self._molecule.atom(i).position
        )
        jk = (
            self._molecule.atom(k).position
            - self._molecule.atom(j).position
        )
        kl = (
            self._molecule.atom(l).position
            - self._molecule.atom(k).position
        )

        ijk_normal = np.cross(ij, jk)
        jkl_normal = np.cross(jk, kl)

        ijk_norm = np.linalg.norm(ijk_normal)
        jkl_norm = np.linalg.norm(jkl_normal)

        if ijk_norm < 1e-12 or jkl_norm < 1e-12:
            raise ValueError(
                "Torsional angle is undefined "
                "for collinear atoms"
            )

        cosine = np.clip(
            np.dot(ijk_normal, jkl_normal)
            / (ijk_norm * jkl_norm),
            -1.0,
            1.0,
        )
        angle = float(np.degrees(np.arccos(cosine)))

        sign = np.dot(
            jk,
            np.cross(ijk_normal, jkl_normal),
        )
        if sign < 0.0:
            angle *= -1.0

        return angle
