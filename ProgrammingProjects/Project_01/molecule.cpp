//
// Created by LiuRuichen on 2026/8/4.
//

#include "molecule.h"
#include "atomic_masses.h"
#include "physical_constants.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <set>




Molecule::Molecule(std::vector<Atom> atoms) : _atoms(std::move(atoms)) {}

Molecule Molecule::from_file(const std::filesystem::path& path) {
    // 创建一个名为file的std::ifstream对象，用path作为构造函数的传入参数
    std::ifstream file{path};
    // 判断传入的路径是否能打开
    if (!file) {
        throw std::runtime_error("Cannot open file " + path.string());
    }

    // 读取原子数目
    std::size_t n_atom;
    if (!(file >> n_atom)) {
        throw std::runtime_error(
            "Failed to read atom count"
        );
    }

    // 创建保存Atom对象的vector
    std::vector<Atom> atoms;
    atoms.reserve(n_atom); // 这是提前让这个数组有足够保存所有atom对象的内存，不写也不错，只是会在后续不断地临时申请更多内存，降低效率

    // 依次读取原子信息
    for (std::size_t i = 0; i < n_atom; i++) {
        int atomic_number;
        double x, y, z;
        file >> atomic_number >> x >> y >> z;
        if (!file) {
            throw std::runtime_error("Failed to read atom at index " + std::to_string(i));
        }

        if (atomic_number <= 0 || static_cast<std::size_t>(atomic_number) >= ATOMIC_MASSES.size()) {
            throw std::runtime_error(
                "Invalid atomic number at index " +
                std::to_string(i)
            );
}

        Atom atom{atomic_number, ATOMIC_MASSES[atomic_number], Eigen::Vector3d{x, y, z}};
        atoms.push_back(atom);
    }

    return Molecule{std::move(atoms)};
}

std::size_t Molecule::atomCount() const
{
    return _atoms.size();
}

const Atom& Molecule::atom(std::size_t index) const {
    return _atoms[index];
}

double Molecule::get_distance_by_ij_index(std::size_t i, std::size_t j) const {
    if (i >= _atoms.size() || j >= _atoms.size()) {
        throw std::out_of_range("Atom index out of range");
    }

    const double d = (_atoms[i].position - _atoms[j].position).norm();
    return d;
}

void Molecule::print_distances() const {
    std::cout << "Interatomic distances (bohr):\n";
    for (std::size_t i = 0; i < _atoms.size(); i++) {
        for (std::size_t j = 0; j < i; j++) {
            std::cout << i <<" " << j  << " "<< get_distance_by_ij_index(i, j) << "\n";
        }
    }
}

double Molecule::get_angele_by_ijk_index(std::size_t i, std::size_t j, std::size_t k) const {
    if (i >= _atoms.size() || j >= _atoms.size() || k >= _atoms.size()) {
        throw std::out_of_range("Atom index out of range");
    }

    const Eigen::Vector3d v_ij = _atoms[i].position - _atoms[j].position;
    const Eigen::Vector3d v_ik = _atoms[k].position - _atoms[j].position;

    const double norm_ij = v_ij.norm();
    const double norm_ik = v_ik.norm();

    if (norm_ij == 0 || norm_ik == 0) {
        throw std::runtime_error("Cannot calculate angle with coincident atoms");
    }

    double cos_theta = v_ij.dot(v_ik) / (norm_ij * norm_ik);
    cos_theta = std::clamp(cos_theta, -1.0, 1.0);

    const double theta = std::acos(cos_theta) * 180.0 / PI;
    return theta;
}

void Molecule::print_angles() const {
    std::cout << "Bond angles:\n";
    for (std::size_t i = 0; i < _atoms.size(); i++) {
        for (std::size_t j = i + 1; j < _atoms.size(); j++) {
            for (std::size_t k = j + 1; k < _atoms.size(); k++) {
                if (get_distance_by_ij_index(i, j) <= 4.0 && get_distance_by_ij_index(j, k) <= 4.0) {
                    std::cout << i << " " << j << " " << k << " ";
                    std::cout << std::fixed << std::setprecision(6) << get_angele_by_ijk_index(i, j, k) << "\n";
                }
            }
        }
    }
}

double Molecule::get_out_of_plane_angle_by_ijkl_index(std::size_t i, std::size_t j, std::size_t k,
                                                      std::size_t l) const {
    if (i >= _atoms.size() || j >= _atoms.size() || k >= _atoms.size() || l >= _atoms.size()) {
        throw std::out_of_range("Atom index out of range");
    }

    Eigen::Vector3d v_ik{_atoms[i].position - _atoms[k].position};
    Eigen::Vector3d v_jk{_atoms[j].position - _atoms[k].position};
    Eigen::Vector3d v_lk{_atoms[l].position - _atoms[k].position};

    Eigen::Vector3d n_jkl = v_jk.cross(v_lk); // 求jkl三点确定的面的法向量

    const double norm_ik = v_ik.norm();
    const double norm_jkl = n_jkl.norm();

    if (norm_ik == 0) {
        throw std::runtime_error("Out-of-plane atom coincides with the central atom");
    }

    if (norm_jkl < 1e-12) {
        return 0.0;
    }

    double sin_theta = v_ik.dot(n_jkl) / (norm_ik * norm_jkl);
    sin_theta = std::clamp(sin_theta, -1.0, 1.0);

    const double theta = std::asin(sin_theta) * 180.0 / PI;
    return theta;
}

void Molecule::print_out_of_plane_angles() const {
    std::cout << "Out-of-plane angles:\n";
    for (std::size_t i = 0; i < _atoms.size(); i++) {
        for (std::size_t j = 0; j < _atoms.size(); j++) {
            for (std::size_t k = 0; k < _atoms.size(); k++) {
                for (std::size_t l = 0; l < _atoms.size(); l++) {
                    if (std::set<std::size_t>{i, j, k, l}.size() < 4) {
                        continue;
                    }
                    if (j <= l) {continue;}
                    // 作业规定：三个原子都必须位于中心原子 k
                    // 的 4.0 bohr 范围内
                    if (
                        get_distance_by_ij_index(i, k) < 4.0 &&
                        get_distance_by_ij_index(j, k) < 4.0 &&
                        get_distance_by_ij_index(k, l) < 4.0
                    ) {
                        const double angle =
                            get_out_of_plane_angle_by_ijkl_index(
                                i, j, k, l
                            );
                        std::cout << i << " " << j << " " << k << " " << l << " " << angle << "\n";

                    }
                }
            }
        }
    }
}

double Molecule::get_torsional_angle_by_ijkl_index(
    std::size_t i, std::size_t j, std::size_t k, std::size_t l
) const {
    if (i >= _atoms.size() || j >= _atoms.size() ||
        k >= _atoms.size() || l >= _atoms.size()) {
        throw std::out_of_range("Atom index out of range");
    }

    const Eigen::Vector3d v_ij = _atoms[j].position - _atoms[i].position;
    const Eigen::Vector3d v_jk = _atoms[k].position - _atoms[j].position;
    const Eigen::Vector3d v_kl = _atoms[l].position - _atoms[k].position;

    const Eigen::Vector3d n_ijk = v_ij.cross(v_jk);
    const Eigen::Vector3d n_jkl = v_jk.cross(v_kl);

    const double norm_n_ijk = n_ijk.norm();
    const double norm_n_jkl = n_jkl.norm();

    if (norm_n_ijk < 1e-12 || norm_n_jkl < 1e-12) {
        throw std::runtime_error(
            "Torsional angle is undefined for collinear atoms"
        );
    }

    double cos_tau = n_ijk.dot(n_jkl) / (norm_n_ijk * norm_n_jkl);
    cos_tau = std::clamp(cos_tau, -1.0, 1.0);

    double tau = std::acos(cos_tau) * 180.0 / PI;

    const double sign_tau = v_jk.dot(n_ijk.cross(n_jkl));
    if (sign_tau < 0.0) {
        tau *= -1.0;
    }

    return tau;
}

void Molecule::print_torsional_angles() const {
    std::cout << "Torsional angles:\n";

    for (std::size_t i = 0; i < _atoms.size(); ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            for (std::size_t k = 0; k < j; ++k) {
                for (std::size_t l = 0; l < k; ++l) {
                    if (
                        get_distance_by_ij_index(i, j) < 4.0 &&
                        get_distance_by_ij_index(j, k) < 4.0 &&
                        get_distance_by_ij_index(k, l) < 4.0
                    ) {
                        const double angle = get_torsional_angle_by_ijkl_index(i, j, k, l);

                        std::cout
                            << std::setw(2) << i << "-"
                            << std::setw(2) << j << "-"
                            << std::setw(2) << k << "-"
                            << std::setw(2) << l << " "
                            << std::setw(10) << std::fixed << std::setprecision(6)
                            << angle << '\n';
                    }
                }
            }
        }
    }
}

Eigen::Vector3d Molecule::get_center_of_mass() const {
    double X, Y, Z;
    double sum_Mx{0.0}, sum_My{0.0}, sum_Mz{0.0}, sum_M{0.0};
    for (std::size_t i = 0; i < _atoms.size(); ++i) {
        sum_Mx += _atoms[i].atomic_mass * _atoms[i].position[0];
        sum_My += _atoms[i].atomic_mass * _atoms[i].position[1];
        sum_Mz += _atoms[i].atomic_mass * _atoms[i].position[2];
        sum_M += _atoms[i].atomic_mass;
    }
    X = sum_Mx / sum_M;
    Y = sum_My / sum_M;
    Z = sum_Mz / sum_M;
    return Eigen::Vector3d{X, Y, Z};
}

Eigen::Matrix3d Molecule::get_moment_of_inertia_tensor() const {
    Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
    Eigen::Vector3d center_of_mass = get_center_of_mass();
    for (const Atom& atom : _atoms) {
        Eigen::Vector3d centered_positions = atom.position - center_of_mass;
        I(0, 0) += atom.atomic_mass * (centered_positions[1] * centered_positions[1] + centered_positions[2] * centered_positions[2]);
        I(1, 1) += atom.atomic_mass * (centered_positions[0] * centered_positions[0] + centered_positions[2] * centered_positions[2]);
        I(2, 2) += atom.atomic_mass * (centered_positions[0] * centered_positions[0] + centered_positions[1] * centered_positions[1]);

        I(0, 1) -= atom.atomic_mass * centered_positions[0] * centered_positions[1];
        I(0, 2) -= atom.atomic_mass * centered_positions[0] * centered_positions[2];
        I(1, 2) -= atom.atomic_mass * centered_positions[1] * centered_positions[2];
    }
    I(1, 0) = I(0, 1);
    I(2, 0) = I(0, 2);
    I(2, 1) = I(1, 2);
    return I;
}

Eigen::Vector3d Molecule::get_principal_moments_of_inertia() const {
    const Eigen::Matrix3d inertia = get_moment_of_inertia_tensor();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(inertia);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error(
            "Failed to diagonalize inertia tensor"
        );
    }

    return solver.eigenvalues();
}

std::string Molecule::get_rotor_type() const {
    double tolerance = 1e-4;

    if (_atoms.size() == 2) {
        return "diatomic";
    }

    const Eigen::Vector3d I_ABC = get_principal_moments_of_inertia();

    if (std::abs(I_ABC[0]) < tolerance) {
        return "linear polyatomic";
    }

    const bool AB_equal = (std::abs(I_ABC[0] - I_ABC[1]) < tolerance);
    const bool BC_equal = (std::abs(I_ABC[1] - I_ABC[2]) < tolerance);

    if (AB_equal && BC_equal) {
        return "spherical top";
    }

    if (AB_equal && !BC_equal) {
        return "oblate symmetric top";
    }

    if (!AB_equal && BC_equal) {
        return "prolate symmetric top";
    }
    return "asymmetric top";
}
Eigen::Vector3d Molecule::get_principal_moments_amu_angstrom2() const {
    const Eigen::Vector3d moments = get_principal_moments_of_inertia();
    return moments * (BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM);
}

Eigen::Vector3d Molecule::get_principal_moments_g_cm2() const {
    const Eigen::Vector3d moments = get_principal_moments_of_inertia();
    return moments * AMU_TO_GRAM * (BOHR_TO_CM * BOHR_TO_CM);
}

Eigen::Vector3d Molecule::get_rotational_constants_mhz() const {
    const Eigen::Vector3d moments = get_principal_moments_g_cm2();
    Eigen::Vector3d constants;

    for (Eigen::Index i = 0; i < moments.size(); ++i) {
        if (moments[i] <= 0.0) {
            throw std::runtime_error(
                "Cannot calculate rotational constant from zero moment of inertia"
            );
        }

        const double constant_hz = PLANCK_CONSTANT_ERG_S /
            (8.0 * PI * PI * moments[i]);
        constants[i] = constant_hz / 1e6;
    }

    return constants;
}

Eigen::Vector3d Molecule::get_rotational_constants_cm_inverse() const {
    const Eigen::Vector3d constants_mhz = get_rotational_constants_mhz();
    return constants_mhz * 1e6 / SPEED_OF_LIGHT_CM_S;
}
