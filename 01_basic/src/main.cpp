#include <iostream>
#include <numeric>
#include <vector>

int main() {
    const std::vector<double> orbital_energies{-0.5782, -0.2247, 0.1351, 0.3114};
    const double sum =
        std::accumulate(orbital_energies.begin(), orbital_energies.end(), 0.0);
    const double mean = sum / static_cast<double>(orbital_energies.size());

    std::cout << "Learning C++ for computational chemistry\n";
    std::cout << "Number of toy orbital energies: " << orbital_energies.size() << '\n';
    std::cout << "Mean energy: " << mean << " hartree\n";

    return 0;
}
