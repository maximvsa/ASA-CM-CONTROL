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

## Visualize flowsheet + inspect DOF

Generate a connectivity diagram and a DOF report from the current model scaffold:

```bash
python scripts/visualize_and_report_dof.py
```

Artifacts are written to:
- `artifacts/flowsheet/steady_state_flowsheet.mmd`
- `artifacts/flowsheet/steady_state_dof_report.txt`

Optional: include current stream flow values (if available from current model state):

```bash
python scripts/visualize_and_report_dof.py --with-values
```

### Launch IDAES built-in visualizer

If `idaes_ui` is installed:

```bash
python scripts/visualize_idaes_ui.py
```

Custom title/save path:

```bash
python scripts/visualize_idaes_ui.py --title "ASA CM Demo" --save asa_cm_demo_vis.json --save-dir artifacts/flowsheet
```