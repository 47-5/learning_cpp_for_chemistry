//
// Created by LiuRuichen on 2026/8/4.
//

#ifndef LEARNING_CPP_FOR_CHEMISTRY_MOLECULE_H
#define LEARNING_CPP_FOR_CHEMISTRY_MOLECULE_H

#include <cstddef>
#include <filesystem>
#include <vector>
#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

struct Atom {
    int atomic_number;
    double atomic_mass;
    Eigen::Vector3d position;
};

class Molecule {
public:
    explicit Molecule(std::vector<Atom> atoms);

    static Molecule from_file(const std::filesystem::path& path);

    [[nodiscard]] std::size_t atomCount() const;
    [[nodiscard]] const Atom& atom(std::size_t index) const;

    [[nodiscard]] double get_distance_by_ij_index(std::size_t i, std::size_t j) const;
    void print_distances() const;

    [[nodiscard]] double get_angele_by_ijk_index(std::size_t i, std::size_t j, std::size_t k) const;
    void print_angles() const;

    [[nodiscard]] double get_out_of_plane_angle_by_ijkl_index(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    void print_out_of_plane_angles() const;

    [[nodiscard]] double get_torsional_angle_by_ijkl_index(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    void print_torsional_angles() const;

    [[nodiscard]] Eigen::Vector3d get_center_of_mass() const;
    [[nodiscard]] Eigen::Matrix3d get_moment_of_inertia_tensor() const;
    [[nodiscard]] Eigen::Vector3d get_principal_moments_of_inertia() const;
    [[nodiscard]] std::string get_rotor_type() const;

    [[nodiscard]] Eigen::Vector3d get_principal_moments_amu_angstrom2() const;
    [[nodiscard]] Eigen::Vector3d get_principal_moments_g_cm2() const;
    [[nodiscard]] Eigen::Vector3d get_rotational_constants_mhz() const;
    [[nodiscard]] Eigen::Vector3d get_rotational_constants_cm_inverse() const;

private:
    std::vector<Atom> _atoms;
};

#endif // LEARNING_CPP_FOR_CHEMISTRY_MOLECULE_H
