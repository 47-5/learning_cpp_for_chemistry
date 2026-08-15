#include "minihf/integrals.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace minihf {
namespace {

struct MatrixEntry {
    Eigen::Index row;
    Eigen::Index column;
    double value;
};

std::ifstream open_input(const std::filesystem::path& path)
{
    std::ifstream input{path};
    if (!input) {
        throw std::runtime_error(
            "Cannot open required input file: " + path.string());
    }
    return input;
}

Eigen::MatrixXd read_symmetric_matrix(
    const std::filesystem::path& path,
    Eigen::Index expected_size = 0)
{
    std::ifstream input = open_input(path);
    std::vector<MatrixEntry> entries;
    Eigen::Index maximum_index = 0;

    long long file_row = 0;
    long long file_column = 0;
    double value = 0.0;

    while (input >> file_row >> file_column >> value) {
        if (file_row <= 0 || file_column <= 0) {
            throw std::runtime_error(
                "AO indices must start at 1 in " + path.string());
        }

        const Eigen::Index row =
            static_cast<Eigen::Index>(file_row - 1);
        const Eigen::Index column =
            static_cast<Eigen::Index>(file_column - 1);

        maximum_index = std::max(
            maximum_index,
            std::max(row, column) + 1);
        entries.push_back(MatrixEntry{row, column, value});
    }

    if (!input.eof()) {
        throw std::runtime_error(
            "Invalid matrix data in " + path.string());
    }
    if (entries.empty()) {
        throw std::runtime_error(
            "No matrix data found in " + path.string());
    }

    const Eigen::Index size =
        expected_size == 0 ? maximum_index : expected_size;
    if (maximum_index > size) {
        throw std::runtime_error(
            "AO index exceeds the expected basis size in " +
            path.string());
    }

    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(size, size);
    for (const MatrixEntry& entry : entries) {
        matrix(entry.row, entry.column) = entry.value;
        matrix(entry.column, entry.row) = entry.value;
    }
    return matrix;
}

double read_scalar(const std::filesystem::path& path)
{
    std::ifstream input = open_input(path);
    double value = 0.0;
    if (!(input >> value)) {
        throw std::runtime_error(
            "Expected a number in " + path.string());
    }
    return value;
}

int read_electron_count(
    const std::filesystem::path& geometry_path,
    int molecular_charge)
{
    std::ifstream input = open_input(geometry_path);

    int atom_count = 0;
    if (!(input >> atom_count) || atom_count <= 0) {
        throw std::runtime_error(
            "Invalid atom count in " + geometry_path.string());
    }

    int total_nuclear_charge = 0;
    for (int atom = 0; atom < atom_count; ++atom) {
        double nuclear_charge = 0.0;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

        if (!(input >> nuclear_charge >> x >> y >> z)) {
            throw std::runtime_error(
                "Not enough atom rows in " + geometry_path.string());
        }

        const long long rounded_charge = std::llround(nuclear_charge);
        if (std::abs(
                nuclear_charge - static_cast<double>(rounded_charge))
            > 1.0e-12) {
            throw std::runtime_error(
                "Non-integral nuclear charge in " +
                geometry_path.string());
        }
        total_nuclear_charge += static_cast<int>(rounded_charge);
    }

    return total_nuclear_charge - molecular_charge;
}

std::size_t dense_eri_size(Eigen::Index basis_size)
{
    const std::size_t size = static_cast<std::size_t>(basis_size);
    return size * size * size * size;
}

std::size_t flat_eri_index(
    Eigen::Index mu,
    Eigen::Index nu,
    Eigen::Index lambda,
    Eigen::Index sigma,
    Eigen::Index basis_size)
{
    const std::size_t size = static_cast<std::size_t>(basis_size);
    return (((static_cast<std::size_t>(mu) * size
              + static_cast<std::size_t>(nu))
             * size
             + static_cast<std::size_t>(lambda))
            * size
            + static_cast<std::size_t>(sigma));
}

std::vector<double> read_eri(
    const std::filesystem::path& path,
    Eigen::Index basis_size)
{
    std::ifstream input = open_input(path);
    std::vector<double> eri(dense_eri_size(basis_size), 0.0);

    long long file_i = 0;
    long long file_j = 0;
    long long file_k = 0;
    long long file_l = 0;
    double value = 0.0;

    while (input >> file_i >> file_j >> file_k >> file_l >> value) {
        const std::array<long long, 4> file_indices{
            file_i,
            file_j,
            file_k,
            file_l};
        for (const long long index : file_indices) {
            if (index <= 0 || index > basis_size) {
                throw std::runtime_error(
                    "ERI index is outside the AO basis in " +
                    path.string());
            }
        }

        const Eigen::Index i = static_cast<Eigen::Index>(file_i - 1);
        const Eigen::Index j = static_cast<Eigen::Index>(file_j - 1);
        const Eigen::Index k = static_cast<Eigen::Index>(file_k - 1);
        const Eigen::Index l = static_cast<Eigen::Index>(file_l - 1);

        const std::array<std::array<Eigen::Index, 4>, 8> permutations{{
            {{i, j, k, l}},
            {{j, i, k, l}},
            {{i, j, l, k}},
            {{j, i, l, k}},
            {{k, l, i, j}},
            {{l, k, i, j}},
            {{k, l, j, i}},
            {{l, k, j, i}},
        }};

        for (const auto& indices : permutations) {
            eri[flat_eri_index(
                indices[0],
                indices[1],
                indices[2],
                indices[3],
                basis_size)] = value;
        }
    }

    if (!input.eof()) {
        throw std::runtime_error(
            "Invalid ERI data in " + path.string());
    }
    return eri;
}

} // namespace

AOIntegrals::AOIntegrals(
    int electron_count,
    double nuclear_repulsion,
    Eigen::MatrixXd overlap,
    Eigen::MatrixXd kinetic,
    Eigen::MatrixXd nuclear_attraction,
    std::vector<double> electron_repulsion)
    : electron_count_{electron_count},
      nuclear_repulsion_{nuclear_repulsion},
      overlap_{std::move(overlap)},
      kinetic_{std::move(kinetic)},
      nuclear_attraction_{std::move(nuclear_attraction)},
      electron_repulsion_{std::move(electron_repulsion)}
{
    validate();
}

AOIntegrals AOIntegrals::from_directory(
    const std::filesystem::path& input_directory,
    int molecular_charge)
{
    if (!std::filesystem::is_directory(input_directory)) {
        throw std::invalid_argument(
            "Input directory does not exist: " +
            input_directory.string());
    }

    Eigen::MatrixXd overlap =
        read_symmetric_matrix(input_directory / "s.dat");
    const Eigen::Index basis_size = overlap.rows();

    return AOIntegrals{
        read_electron_count(
            input_directory / "geom.dat",
            molecular_charge),
        read_scalar(input_directory / "enuc.dat"),
        std::move(overlap),
        read_symmetric_matrix(
            input_directory / "t.dat",
            basis_size),
        read_symmetric_matrix(
            input_directory / "v.dat",
            basis_size),
        read_eri(input_directory / "eri.dat", basis_size)};
}

int AOIntegrals::electron_count() const noexcept
{
    return electron_count_;
}

Eigen::Index AOIntegrals::basis_size() const noexcept
{
    return overlap_.rows();
}

double AOIntegrals::nuclear_repulsion() const noexcept
{
    return nuclear_repulsion_;
}

const Eigen::MatrixXd& AOIntegrals::overlap() const noexcept
{
    return overlap_;
}

const Eigen::MatrixXd& AOIntegrals::kinetic() const noexcept
{
    return kinetic_;
}

const Eigen::MatrixXd& AOIntegrals::nuclear_attraction() const noexcept
{
    return nuclear_attraction_;
}

Eigen::MatrixXd AOIntegrals::core_hamiltonian() const
{
    return kinetic_ + nuclear_attraction_;
}

double AOIntegrals::eri(
    Eigen::Index mu,
    Eigen::Index nu,
    Eigen::Index lambda,
    Eigen::Index sigma) const
{
    return electron_repulsion_.at(
        eri_index(mu, nu, lambda, sigma));
}

std::size_t AOIntegrals::eri_index(
    Eigen::Index mu,
    Eigen::Index nu,
    Eigen::Index lambda,
    Eigen::Index sigma) const
{
    const Eigen::Index size = basis_size();
    for (const Eigen::Index index : {mu, nu, lambda, sigma}) {
        if (index < 0 || index >= size) {
            throw std::out_of_range("AO ERI index is outside the basis");
        }
    }
    return flat_eri_index(mu, nu, lambda, sigma, size);
}

void AOIntegrals::validate() const
{
    if (electron_count_ <= 0) {
        throw std::invalid_argument(
            "Electron count must be positive");
    }
    if (overlap_.rows() == 0 || overlap_.rows() != overlap_.cols()) {
        throw std::invalid_argument(
            "Overlap matrix must be non-empty and square");
    }

    const Eigen::Index size = overlap_.rows();
    const auto expected_shape = [size](const Eigen::MatrixXd& matrix) {
        return matrix.rows() == size && matrix.cols() == size;
    };

    if (!expected_shape(kinetic_)
        || !expected_shape(nuclear_attraction_)) {
        throw std::invalid_argument(
            "All one-electron matrices must have the same shape");
    }

    if (!overlap_.isApprox(overlap_.transpose())
        || !kinetic_.isApprox(kinetic_.transpose())
        || !nuclear_attraction_.isApprox(
            nuclear_attraction_.transpose())) {
        throw std::invalid_argument(
            "One-electron matrices must be symmetric");
    }

    if (electron_repulsion_.size() != dense_eri_size(size)) {
        throw std::invalid_argument(
            "Dense ERI storage has the wrong size");
    }
}

} // namespace minihf
