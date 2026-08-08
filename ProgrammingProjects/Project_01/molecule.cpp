//
// Created by LiuRuichen on 2026/8/4.
//

#include "molecule.h"

#include <stdexcept>
#include <utility>

Molecule::Molecule(std::vector<Atom> atoms) : _atoms(std::move(atoms)) {}

std::size_t Molecule::atomCount() const
{
    return _atoms.size();
}

double Molecule::distanceBetweenAtoms(std::size_t atom1_index, std::size_t atom2_index) const
{
    if (atom1_index >= _atoms.size() || atom2_index >= _atoms.size()) {
        throw std::out_of_range("Atom index out of range");
    }

    const double d = (_atoms[atom1_index].position - _atoms[atom2_index].position).norm();
    return d;
}
