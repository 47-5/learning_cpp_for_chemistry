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
python3 -m Project_02.main \
    ../resources/Project_02/input/h2o_geom.txt \
    ../resources/Project_02/input/h2o_hessian.txt
```

Both projects import the same domain model and analyzers from
`chemistry/`. Input data and theory notes are shared with the C++
implementations through `../resources/`.
