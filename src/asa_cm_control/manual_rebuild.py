from pathlib import Path
import sys

SRC_DIR = Path(__file__).resolve().parents[1]
if str(SRC_DIR) not in sys.path:
	sys.path.insert(0, str(SRC_DIR))

from pyomo.environ import ConcreteModel, SolverFactory
from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.base.generic_reaction import GenericReactionParameterBlock
from idaes.models.unit_models import CSTR
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog

from asa_cm_control.props.thermo_config import configuration as thermo_config
from asa_cm_control.props.reaction_config import configuration as reaction_config

def build_model():
	m = ConcreteModel()
	m.fs = FlowsheetBlock(dynamic=False)

	m.fs.thermo_params = GenericParameterBlock(**thermo_config)
	m.fs.reaction_params = GenericReactionParameterBlock(
		property_package=m.fs.thermo_params,
		**reaction_config
	)

	m.fs.cstr = CSTR(
		property_package=m.fs.thermo_params,
		reaction_package=m.fs.reaction_params,
	)

	return m


def set_operating_conditions(m):
	m.fs.cstr.inlet.flow_mol.fix(1.0)
	m.fs.cstr.inlet.temperature.fix(350)
	m.fs.cstr.inlet.pressure.fix(101325)

	epsilon = 1e-8
	major = 1 - 3 * epsilon
	m.fs.cstr.inlet.mole_frac_comp[0, "salicylic_acid"].fix(0.495 * major)
	m.fs.cstr.inlet.mole_frac_comp[0, "acetic_anhydride"].fix(0.495 * major)
	m.fs.cstr.inlet.mole_frac_comp[0, "sulfuric_acid"].fix(0.01 * major)
	m.fs.cstr.inlet.mole_frac_comp[0, "aspirin"].fix(epsilon)
	m.fs.cstr.inlet.mole_frac_comp[0, "acetic_acid"].fix(epsilon)
	m.fs.cstr.inlet.mole_frac_comp[0, "water"].fix(epsilon)

	m.fs.cstr.volume.fix(10.0)


def initialize_model(m):
	print("Degrees of Freedom =", degrees_of_freedom(m.fs.cstr))
	m.fs.cstr.initialize(outlvl=idaeslog.INFO)


def solve_model(m):

	solver = SolverFactory("ipopt")
	return solver.solve(m, tee=True)


def main():
	m = build_model()
	set_operating_conditions(m)
	initialize_model(m)
	status = solve_model(m)
	m.fs.cstr.report()
	return m, status


if __name__ == "__main__":
	main()