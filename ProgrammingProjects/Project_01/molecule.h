//
// Created by LiuRuichen on 2026/8/4.
//

#ifndef LEARNING_CPP_FOR_CHEMISTRY_MOLECULE_H
#define LEARNING_CPP_FOR_CHEMISTRY_MOLECULE_H

#include <cstddef>
#include <vector>

#include <Eigen/Core>

struct Atom {
    int atomic_number;
    double atomic_mass;
    Eigen::Vector3d position;
};

class Molecule {
public:
    explicit Molecule(std::vector<Atom> atoms);

    [[nodiscard]] std::size_t atomCount() const;
    [[nodiscard]] double distanceBetweenAtoms(std::size_t atom1_index, std::size_t atom2_index) const;

private:
    std::vector<Atom> _atoms;
};

#endif // LEARNING_CPP_FOR_CHEMISTRY_MOLECULE_H
