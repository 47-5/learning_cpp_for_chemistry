from typing import List
import numpy as np

from atomic_mass import ATOMIC_MASSES
from physical_constants import BOHR_TO_ANGSTROM, BOHR_TO_CM, AMU_TO_GRAM, PLANCK_CONSTANT, SPEED_OF_LIGHT


class Atom:
    def __init__(self, atom_number, x, y, z):
        self.atom_number = int(atom_number)
        if not 1 <= self.atom_number < len(ATOMIC_MASSES):
            raise ValueError(
                f"Invalid atomic number: {self.atom_number}"
            )

        self.x, self.y, self.z = float(x), float(y), float(z)
        self.mass = ATOMIC_MASSES[self.atom_number]


class Molecule:
    def __init__(self, atoms: List[Atom]):
        self.atoms = atoms

    @staticmethod
    def from_file(file_path):
        f = open(file_path).readlines()
        N = int(f[0])
        atoms = []
        for l in f[1:N+1]:
            atom_number, x, y, z = l.strip().split()
            atoms.append(Atom(atom_number, x, y, z))
        if len(atoms) != N:
            raise ValueError(
                f"Expected {N} atoms, but found {len(atoms)}"
            )
        return Molecule(atoms=atoms)

    @property
    def N(self):
        return len(self.atoms)

    @property
    def positions(self):
        return np.array([[i.x, i.y, i.z] for i in self.atoms])

    @property
    def atomic_numbers(self):
        return np.array([i.atom_number for i in self.atoms])

    @property
    def masses(self):
        return np.array([i.mass for i in self.atoms])

    def get_distance_by_ij_index(self, i, j):
        assert 0 <= i < self.N, print('atom index out of N_atoms')
        assert 0 <= j < self.N, print('atom index out of N_atoms')

        distance = np.linalg.norm(self.positions[i] - self.positions[j])
        return distance

    def get_distances_matrix(self):
        D = np.array([[self.get_distance_by_ij_index(i, j) for j in range(self.N)] for i in range(self.N)])
        return D

    def print_distances(self):
        print('Interatomic distances:')
        D = self.get_distances_matrix()
        for i in range(self.N):
            for j in range(self.N):
                if i > j:
                    print(f"{i} {j} {D[i,j]}")

    def get_angele_by_ijk_index(self, i, j, k):
        assert 0 <= i < self.N, print('atom index out of N_atoms')
        assert 0 <= j < self.N, print('atom index out of N_atoms')
        assert 0 <= k < self.N, print('atom index out of N_atoms')

        u = self.positions[i] - self.positions[j]
        v = self.positions[k] - self.positions[j]

        norm_u = np.linalg.norm(u)
        norm_v = np.linalg.norm(v)

        if norm_u == 0.0 or norm_v == 0.0:
            raise ValueError("Cannot calculate an angle with coincident atoms")
        cos_theta = np.dot(u, v) / (norm_u * norm_v)
        cos_theta = np.clip(cos_theta, -1.0, 1.0)

        return np.degrees(np.arccos(cos_theta))

    def print_angles(self):
        print('Bond angles:')
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    if i < j < k and self.get_distance_by_ij_index(i, j) <= 4.0 and self.get_distance_by_ij_index(j, k) <= 4.0:
                        print(f"{i} {j} {k} {self.get_angele_by_ijk_index(i, j, k):10.6f}")

    def get_out_of_plane_angle_by_ijkl_index(self, i, j, k, l):
        assert 0 <= i < self.N, print('atom index out of N_atoms')
        assert 0 <= j < self.N, print('atom index out of N_atoms')
        assert 0 <= k < self.N, print('atom index out of N_atoms')
        assert 0 <= l < self.N, print('atom index out of N_atoms')

        v_ki = self.positions[i] - self.positions[k]
        v_kj = self.positions[j] - self.positions[k]
        v_kl = self.positions[l] - self.positions[k]

        n_jkl = np.cross(v_kj, v_kl)

        norm_v_ki = np.linalg.norm(v_ki)
        norm_n_jkl = np.linalg.norm(n_jkl)

        if norm_v_ki == 0.0:
            raise ValueError(
                "Out-of-plane atom coincides with the central atom"
            )

        if norm_n_jkl < 1e-12:
            return 0.0

        sin_phi = np.dot(v_ki, n_jkl) / (norm_v_ki * norm_n_jkl)

        sin_phi = np.clip(sin_phi, -1.0, 1.0)

        return np.degrees(np.arcsin(sin_phi))

    def print_out_of_plane_angles(self):
        print("Out-of-plane angles:")

        D = self.get_distances_matrix()

        for i in range(self.N):
            for k in range(self.N):
                for j in range(self.N):
                    for l in range(self.N):
                        if len({i, j, k, l}) < 4:
                            continue

                        if j <= l:
                            continue

                        if (
                            D[i, k] < 4.0
                            and D[j, k] < 4.0
                            and D[k, l] < 4.0
                        ):
                            angle = (
                                self.get_out_of_plane_angle_by_ijkl_index(
                                    i, j, k, l
                                )
                            )

                            print(
                                f"{i}-{j}-{k}-{l} "
                                f"{angle:10.6f}"
                            )

    def get_torsional_angle_by_ijkl_index(self, i, j, k, l):
        v_ij = self.positions[j] - self.positions[i]
        v_jk = self.positions[k] - self.positions[j]
        v_kl = self.positions[l] - self.positions[k]

        n_ijk = np.cross(v_ij, v_jk)
        n_jkl = np.cross(v_jk, v_kl)

        norm_n_ijk = np.linalg.norm(n_ijk)
        norm_n_jkl = np.linalg.norm(n_jkl)

        if norm_n_ijk == 0.0 or norm_n_jkl == 0.0:
            raise ValueError(
                "Torsional angle is undefined for collinear atoms"
            )

        cos_tau = np.dot(n_ijk, n_jkl) / (norm_n_ijk * norm_n_jkl)
        cos_tau = np.clip(cos_tau, -1.0, 1.0)
        tau = np.degrees(np.arccos(cos_tau))

        sign_tau = np.dot(v_jk, np.cross(n_ijk, n_jkl))
        if sign_tau < 0:
            tau *= -1

        return tau

    def print_torsional_angles(self):
        D = self.get_distances_matrix()
        for l in range(self.N):
            for k in range(l+1, self.N):
                for j in range(k+1, self.N):
                    for i in range(j+1, self.N):
                        if D[i, j] <= 4.0 and D[j, k] <= 4.0 and D[k, l] <= 4.0:
                            print(f"{i}-{j}-{k}-{l} {self.get_torsional_angle_by_ijkl_index(i, j, k, l)}")

    def get_center_of_mass(self):
        positions = self.positions
        masses = self.masses.reshape(-1, 1)
        M = np.sum(masses)
        X_M = np.sum(masses * positions, axis=0)
        center_of_mass = X_M / M
        return center_of_mass

    def get_moment_of_inertia_tensor(self):
        centered_positions = self.positions - self.get_center_of_mass()

        I = np.zeros((3, 3))

        for i in range(self.N):
            x, y, z = centered_positions[i]
            mass = self.masses[i]

            I[0, 0] += mass * (y ** 2 + z ** 2)
            I[1, 1] += mass * (x ** 2 + z ** 2)
            I[2, 2] += mass * (x ** 2 + y ** 2)

            I[0, 1] -= mass * x * y
            I[0, 2] -= mass * x * z
            I[1, 2] -= mass * y * z

        I[1, 0] = I[0, 1]
        I[2, 0] = I[0, 2]
        I[2, 1] = I[1, 2]

        return I

    def get_principal_moments_of_inertia(self):
        I = self.get_moment_of_inertia_tensor()

        eigenvalues, eigenvectors = np.linalg.eigh(I)

        return eigenvalues

    def get_rotor_type(self):
        I_A, I_B, I_C = self.get_principal_moments_of_inertia()

        tolerance = 1e-4

        if self.N == 2:
            return "diatomic"

        if I_A < tolerance:
            return "linear polyatomic"

        AB_equal = abs(I_A - I_B) < tolerance
        BC_equal = abs(I_B - I_C) < tolerance

        if AB_equal and BC_equal:
            return "spherical top"

        if AB_equal and not BC_equal:
            return "oblate symmetric top"

        if not AB_equal and BC_equal:
            return "prolate symmetric top"

        return "asymmetric top"

    def get_principal_moments_in_different_units(self):
        moments_amu_bohr2 = (
            self.get_principal_moments_of_inertia()
        )

        moments_amu_angstrom2 = (
            moments_amu_bohr2
            * BOHR_TO_ANGSTROM ** 2
        )

        moments_g_cm2 = (
            moments_amu_bohr2
            * AMU_TO_GRAM
            * BOHR_TO_CM ** 2
        )

        return (
            moments_amu_bohr2,
            moments_amu_angstrom2,
            moments_g_cm2,
        )

    def get_rotational_constants(self):
        _, _, moments_g_cm2 = (
            self.get_principal_moments_in_different_units()
        )

        constants_hz = PLANCK_CONSTANT / (
            8 * np.pi ** 2 * moments_g_cm2
        )

        constants_mhz = constants_hz / 1e6

        constants_cm_inverse = (
            constants_hz / SPEED_OF_LIGHT
        )

        return constants_mhz, constants_cm_inverse
