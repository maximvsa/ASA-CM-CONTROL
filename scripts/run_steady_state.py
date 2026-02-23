from __future__ import annotations

import sys
from pathlib import Path

import pyomo.environ as pyo

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from asa_cm_control.models.steady_state import build_model, solve_model, summarize_model


def main() -> int:
    model = build_model()
    results = solve_model(model)

    termination = results.solver.termination_condition
    status = results.solver.status

    if termination != pyo.TerminationCondition.optimal:
        print(f"Solve did not reach optimal termination. status={status}, termination={termination}")
        return 1

    summary = summarize_model(model)
    print("Steady-state model solved successfully.")
    for key, value in summary.items():
        print(f"- {key}: {value:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
