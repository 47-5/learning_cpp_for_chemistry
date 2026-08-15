# Programming projects

The programming projects are organized by language while sharing the same
input data and theory notes.

```text
ProgrammingProjects/
├── cpp/
│   ├── common/             # Reusable chemistry_core C++ library
│   ├── Project_01/         # C++ Project 1 entry point
│   └── Project_02/         # C++ Project 2 entry point
├── python/
│   ├── chemistry/          # Reusable Python package
│   ├── Project_01/         # Python Project 1 entry point
│   ├── Project_02/         # Python Project 2 entry point
│   └── Project_03/         # Python RHF/SCF scaffold
└── resources/
    ├── Project_01/         # Project 1 inputs and notes
    ├── Project_02/         # Project 2 inputs and notes
    └── Project_03/         # Project 3 AO integrals and theory
```

## Responsibility boundaries

`Molecule` owns only the basic molecular data: atoms, atomic masses,
Cartesian positions, and geometry-file loading.

Analysis algorithms are separate objects initialized from a molecule:

```cpp
const chemistry::Molecule molecule =
    chemistry::Molecule::from_file(path);

const chemistry::GeometryAnalyzer geometry{molecule};
const chemistry::RotationalAnalyzer rotation{molecule};
```

The analyzers hold a reference to the molecule, so the molecule must outlive
them. Project entry points are responsible for command-line parsing and
printing; analyzers return numerical results without writing to standard
output.

## C++ targets

Both executables depend on the same `chemistry_core` library:

```text
chemistry_core
     ↑
 ┌───┴────┐
project01 project02
```

From the repository root:

```bash
cmake -S . -B build -G Ninja
cmake --build build --parallel

./build/ProgrammingProjects/cpp/Project_01/project01 \
    ProgrammingProjects/resources/Project_01/input/allene.dat

./build/ProgrammingProjects/cpp/Project_02/project02 \
    ProgrammingProjects/resources/Project_02/input/h2o_geom.txt \
    ProgrammingProjects/resources/Project_02/input/h2o_hessian.txt
```

See `python/README.md` for the Python commands.
