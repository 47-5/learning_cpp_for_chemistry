# Project 3 C++ scaffold

This is the C++17/Eigen version of the closed-shell RHF teaching project.
It mirrors the Python responsibility boundaries without copying Python
syntax into C++.

```text
Project_03/
├── include/minihf/
│   ├── integrals.h        # AOIntegrals and dense ERI access
│   └── rhf_calculator.h   # SCF data objects and RHFCalculator
├── src/
│   ├── integrals.cpp      # Project 3 input readers
│   └── rhf_calculator.cpp # Ordered numerical TODOs
├── tests/
│   └── scaffold_test.cpp  # Input and construction smoke test
├── main.cpp
└── CMakeLists.txt
```

The teaching implementation stores ERIs in a dense one-dimensional vector
and exposes them through

```cpp
integrals.eri(mu, nu, lambda, sigma)
```

which returns \((\mu\nu|\lambda\sigma)\). This keeps the first C++ version
close to the Python four-dimensional array. Eightfold compression can be a
later optimization.

Implement `RHFCalculator` in this order:

1. `build_orthogonalizer()`
2. `solve_roothaan()`
3. `build_density()`
4. `build_fock()`
5. `electronic_energy()`
6. `run()`

From the repository root:

```bash
cmake -S . -B build -G Ninja
cmake --build build --target project03 project03_scaffold_test

./build/ProgrammingProjects/cpp/Project_03/project03 \
  ProgrammingProjects/resources/Project_03/input/h2o/STO-3G

ctest --test-dir build -R project03_scaffold --output-on-failure
```
