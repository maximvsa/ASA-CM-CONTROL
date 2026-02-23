# ASA-CM-CONTROL

Initial scaffold for an IDAES steady-state model of the ASA continuous process.

## Python / IDAES compatibility

- Use Python 3.10 for this repository.
- The provided environment file is pinned to Python 3.10 and includes `idaes-pse`.
- In VS Code, select the interpreter for your 3.10 env (for example `.../miniforge3/envs/asa-cm-control-env/python.exe`).
- Run `scripts/setup_windows.bat` from an Anaconda/Miniforge-enabled shell where `conda` is available.

## Repository layout

```text
ASA-CM-CONTROL/
├─ asa_cm_control_flowsheet_graph.yaml
├─ environment.yaml
├─ scripts/
│  ├─ setup_windows.bat
│  └─ run_steady_state.py
├─ src/
│  └─ asa_cm_control/
│     ├─ __init__.py
│     ├─ data/
│     │  └─ README.md
│     └─ models/
│        ├─ __init__.py
│        └─ steady_state.py
└─ tests/
	└─ test_steady_state_smoke.py
```

## Run first steady-state build

From the repository root (after selecting the 3.10 interpreter in VS Code):

```bash
python scripts/run_steady_state.py
```

Expected output includes:
- `Steady-state model solved successfully.`
- `feed_rate_mol_per_s`, `conversion`, and `product_rate_mol_per_s`

## Run smoke test

```bash
pytest -q
```