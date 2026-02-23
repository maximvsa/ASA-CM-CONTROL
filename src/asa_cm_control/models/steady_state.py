from __future__ import annotations

import pyomo.environ as pyo
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver


def build_model(feed_rate_mol_per_s: float = 1.0, conversion: float = 0.6) -> pyo.ConcreteModel:
    """Build a minimal IDAES steady-state model for initial scaffolding.

    This is intentionally simple so the project has a reliable, runnable baseline.
    """
    model = pyo.ConcreteModel(name="asa_cm_control_steady_state")
    model.fs = FlowsheetBlock(dynamic=False)

    model.fs.feed_rate = pyo.Var(
        initialize=feed_rate_mol_per_s,
        bounds=(0, None),
        units=pyo.units.mol / pyo.units.s,
        doc="Total inlet molar flow",
    )
    model.fs.conversion = pyo.Var(
        initialize=conversion,
        bounds=(0, 1),
        units=pyo.units.dimensionless,
        doc="Lumped ASA formation conversion",
    )
    model.fs.product_rate = pyo.Var(
        initialize=feed_rate_mol_per_s * conversion,
        bounds=(0, None),
        units=pyo.units.mol / pyo.units.s,
        doc="Lumped product molar flow",
    )

    model.fs.product_rate_constraint = pyo.Constraint(
        expr=model.fs.product_rate == model.fs.feed_rate * model.fs.conversion
    )

    model.fs.feed_rate.fix(feed_rate_mol_per_s)
    model.fs.conversion.fix(conversion)

    return model


def solve_model(model: pyo.ConcreteModel, solver_name: str = "ipopt"):
    """Solve the steady-state model with an available NLP solver."""
    solver = get_solver(solver=solver_name)
    results = solver.solve(model, tee=False)
    return results


def summarize_model(model: pyo.ConcreteModel) -> dict[str, float]:
    """Return key solved values as plain floats."""
    return {
        "feed_rate_mol_per_s": pyo.value(model.fs.feed_rate),
        "conversion": pyo.value(model.fs.conversion),
        "product_rate_mol_per_s": pyo.value(model.fs.product_rate),
    }
