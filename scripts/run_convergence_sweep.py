import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pyomo.environ import value, check_optimal_termination
from asa_cm_control.steady_state_demo import build_model, initialize_model_sequential, solve_model


def normalize(comp):
    total = sum(comp.values())
    return {key: val / total for key, val in comp.items()}


candidates = [
    {
        "name": "A_baseline",
        "s410": {
            "flow_mol": 6.0,
            "temperature": 308.15,
            "pressure": 101325.0,
            "comp": {"SA": 0.09, "AA": 0.30, "ASA": 0.40, "AcOH": 0.20, "H2SO4": 0.01},
        },
        "s442": {
            "flow_mol": 2.0,
            "temperature": 305.15,
            "pressure": 101325.0,
            "comp": {"SA": 0.005, "AA": 0.75, "ASA": 0.004, "AcOH": 0.24, "H2SO4": 0.001},
        },
    },
    {
        "name": "B_more_AA_reflux",
        "s410": {
            "flow_mol": 5.5,
            "temperature": 309.15,
            "pressure": 101325.0,
            "comp": {"SA": 0.07, "AA": 0.34, "ASA": 0.38, "AcOH": 0.20, "H2SO4": 0.01},
        },
        "s442": {
            "flow_mol": 2.5,
            "temperature": 304.15,
            "pressure": 101325.0,
            "comp": {"SA": 0.003, "AA": 0.80, "ASA": 0.003, "AcOH": 0.193, "H2SO4": 0.001},
        },
    },
    {
        "name": "C_lower_ASA_recycle",
        "s410": {
            "flow_mol": 4.5,
            "temperature": 307.15,
            "pressure": 101325.0,
            "comp": {"SA": 0.12, "AA": 0.32, "ASA": 0.28, "AcOH": 0.27, "H2SO4": 0.01},
        },
        "s442": {
            "flow_mol": 1.6,
            "temperature": 306.15,
            "pressure": 101325.0,
            "comp": {"SA": 0.008, "AA": 0.68, "ASA": 0.008, "AcOH": 0.303, "H2SO4": 0.001},
        },
    },
    {
        "name": "D_higher_recycle_flow",
        "s410": {
            "flow_mol": 8.0,
            "temperature": 310.15,
            "pressure": 101325.0,
            "comp": {"SA": 0.08, "AA": 0.29, "ASA": 0.42, "AcOH": 0.20, "H2SO4": 0.01},
        },
        "s442": {
            "flow_mol": 3.0,
            "temperature": 304.65,
            "pressure": 101325.0,
            "comp": {"SA": 0.004, "AA": 0.72, "ASA": 0.005, "AcOH": 0.270, "H2SO4": 0.001},
        },
    },
    {
        "name": "E_conservative",
        "s410": {
            "flow_mol": 3.8,
            "temperature": 306.65,
            "pressure": 101325.0,
            "comp": {"SA": 0.14, "AA": 0.33, "ASA": 0.22, "AcOH": 0.30, "H2SO4": 0.01},
        },
        "s442": {
            "flow_mol": 1.2,
            "temperature": 305.65,
            "pressure": 101325.0,
            "comp": {"SA": 0.010, "AA": 0.63, "ASA": 0.010, "AcOH": 0.349, "H2SO4": 0.001},
        },
    },
]

solver_options = {
    "max_iter": 5000,
    "tol": 1e-6,
    "acceptable_tol": 1e-5,
    "print_level": 0,
}

rows = []
for candidate in candidates:
    t0 = 0.0
    s410_comp = normalize(candidate["s410"]["comp"])
    s442_comp = normalize(candidate["s442"]["comp"])

    guess_410 = {
        "flow_mol": {t0: candidate["s410"]["flow_mol"]},
        "temperature": {t0: candidate["s410"]["temperature"]},
        "pressure": {t0: candidate["s410"]["pressure"]},
        "mole_frac_comp": {(t0, key): val for key, val in s410_comp.items()},
    }
    guess_442 = {
        "flow_mol": {t0: candidate["s442"]["flow_mol"]},
        "temperature": {t0: candidate["s442"]["temperature"]},
        "pressure": {t0: candidate["s442"]["pressure"]},
        "mole_frac_comp": {(t0, key): val for key, val in s442_comp.items()},
    }

    print(f"\n=== Candidate {candidate['name']} ===")
    try:
        model = build_model()
        initialize_model_sequential(
            model,
            run_unit_initialization=True,
            iter_lim=5,
            display_report=False,
            tear_guess_s410mix_410split=guess_410,
            tear_guess_s442_440_reflux=guess_442,
        )
        results = solve_model(model, solver_options=solver_options)
        termination = str(results.solver.termination_condition)
        optimal = bool(check_optimal_termination(results))
        purity = float(value(model.fs.final_product_purity))
        margin = float(value(model.fs.net_operating_margin_per_s))
        rows.append((candidate["name"], termination, optimal, purity, margin))
        print(
            f"termination={termination}; optimal={optimal}; "
            f"purity={purity:.6f}; margin={margin:.6f}"
        )
    except Exception as error:
        rows.append((candidate["name"], "error", False, float("nan"), float("nan")))
        print(f"error={error}")

print("\n=== Sweep Summary ===")
print("name,termination,optimal,purity,margin")
for row in rows:
    print(f"{row[0]},{row[1]},{row[2]},{row[3]},{row[4]}")
