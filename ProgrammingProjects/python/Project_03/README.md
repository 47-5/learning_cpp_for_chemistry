# Project 3 Python RHF implementation

This directory is the Python entry point for the closed-shell RHF/SCF
exercise. Reusable data and method classes live in
`chemistry/electronic_structure/`.

## Current structure

```text
ProgrammingProjects/python/
├── chemistry/electronic_structure/
│   ├── integrals.py   # AOIntegrals and Project 3 file readers
│   ├── models.py      # SCFOptions, SCFIteration, RHFResult
│   └── rhf.py         # RHFCalculator and SCF implementation
└── Project_03/
    ├── main.py        # H2O/STO-3G entry point
    └── test_scaffold.py
```

`AOIntegrals` uses a dense four-dimensional ERI array in chemist notation:

```python
eri[mu, nu, lam, sig] == (mu nu | lam sig)
```

This is deliberately simple for the small teaching inputs. ERI compression
can be added later without changing the public `RHFCalculator` interface.

## Run the calculation

From `ProgrammingProjects/python`:

```bash
python -m Project_03.main
python -m unittest Project_03.test_scaffold -v
```

The entry point loads H2O/STO-3G, runs the RHF calculation, prints every
SCF iteration, and reports the converged electronic and total energies.

## Numerical structure

`RHFCalculator` keeps the numerical steps separate in this order:

1. `_build_orthogonalizer()`
2. `_solve_roothaan()`
3. `_build_density()`
4. `_build_fock()`
5. `_electronic_energy()`
6. `run()`

Printing remains in `main.py`; the calculator returns arrays and result
objects without writing to standard output.
