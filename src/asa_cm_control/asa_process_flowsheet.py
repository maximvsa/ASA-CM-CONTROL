# To-Do:
# - Add crystallizer unit model and test with same property package
# - 


"""First-pass ASA process flowsheet model assembly and solve workflow.

This module defines reusable functions to build, initialize, solve, and report
the ASA CSTR flowsheet using the custom thermophysical and reaction packages.

Intended entry point:
    run_asa_process.py (repository root)

Run from repository root:
    python run_asa_process.py

The module is intentionally importable and does not provide a direct
``if __name__ == "__main__"`` execution block.
"""

from pyomo.environ import ConcreteModel, SolverFactory
from idaes.core import FlowsheetBlock
from asa_cm_control.props.asa_thermo_property_package import ASAThermoParameterBlock
from asa_cm_control.props.asa_reaction_property_package import ASAReactionParameterBlock
from idaes.models.unit_models import CSTR
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog

def build_flowsheet():
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    
    model.fs.thermo_params = ASAThermoParameterBlock()
    model.fs.reaction_params = ASAReactionParameterBlock(
        property_package=model.fs.thermo_params,
    )
    
    model.fs.cstr = CSTR(
        property_package=model.fs.thermo_params,
        reaction_package=model.fs.reaction_params,
        has_equilibrium_reactions=True
    )
    
    return model


def set_operating_conditions(model):
    model.fs.cstr.inlet.flow_mol.fix(100)
    model.fs.cstr.inlet.temperature.fix(325.0)
    model.fs.cstr.inlet.pressure.fix(101325)
    
    epsilon = 1e-8
    major = 1 - 5 * epsilon
    model.fs.cstr.inlet.mole_frac_comp[0, "salicylic_acid"].fix(0.19 * major)
    model.fs.cstr.inlet.mole_frac_comp[0, "acetic_anhydride"].fix(0.78 * major)
    model.fs.cstr.inlet.mole_frac_comp[0, "aspirin"].fix(epsilon)
    model.fs.cstr.inlet.mole_frac_comp[0, "acetic_acid"].fix(epsilon)
    model.fs.cstr.inlet.mole_frac_comp[0, "water"].fix(0.02 * major)
    
    model.fs.cstr.inlet.mole_frac_comp[0, "H2SO4"].fix(0.01 * major)
    model.fs.cstr.inlet.mole_frac_comp[0, "H_plus"].fix(epsilon)
    model.fs.cstr.inlet.mole_frac_comp[0, "HSO4_minus"].fix(epsilon)
    model.fs.cstr.inlet.mole_frac_comp[0, "SO4_2minus"].fix(epsilon)
    
    model.fs.cstr.volume.fix(0.01)


def initialize_model(model):
    print("Degrees of Freedom =", degrees_of_freedom(model.fs))
    model.fs.cstr.initialize(outlvl=idaeslog.INFO)


def solve_model(model, tee=True):
    solver = SolverFactory("ipopt")
    return solver.solve(model, tee=tee)


def report_results(model):
    model.fs.cstr.report()


def main():
    """Execute the end-to-end flowsheet workflow.

    This function is invoked by the root launcher script ``run_asa_process.py``.
    """
    model = build_flowsheet()
    set_operating_conditions(model)
    initialize_model(model)
    solve_model(model, tee=True)
    report_results(model)