# Python implementations

Run the project entry points from this directory so that Python can find
the shared `chemistry` package. Create an isolated environment and install
its dependency first:

```bash
# One-time WSL setup:
sudo apt install python3-venv python3-pip

cd ProgrammingProjects/python
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Then run:

```bash
cd ProgrammingProjects/python
python3 -m Project_01.main
python3 -m Project_02.main
python3 -m Project_03.main
```

The entry points use input-path constants near the top of each `main.py`,
so they can also be run directly from PyCharm without program arguments.
Change `INPUT_FILE` in Project 1, `GEOMETRY_FILE` and `HESSIAN_FILE`
in Project 2, or `INPUT_DIRECTORY` in Project 3 to select another input.

The projects import the same domain model and analyzers from
`chemistry/`. Input data and theory notes are shared with the C++
implementations through `../resources/`.

Project 3 contains a working teaching implementation of closed-shell RHF.
Its entry point prints the full SCF history and the converged energy, while
the reusable numerical code lives in `chemistry/electronic_structure/`.
