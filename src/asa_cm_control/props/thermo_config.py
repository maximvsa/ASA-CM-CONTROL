from pyomo.environ import units as pyunits
from idaes.core import Component, LiquidPhase
from idaes.models.properties.modular_properties.base.generic_property import StateIndex
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.pure.ConstantProperties import Constant

configuration = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "components": {
        "salicylic_acid": {
            "type": Component,
            "cp_mol_liq_comp": Constant.cp_mol_liq_comp,
            "enth_mol_liq_comp": Constant.enth_mol_liq_comp,
            "dens_mol_liq_comp": Constant.dens_mol_liq_comp,
            "parameter_data": {
                "mw": (0.13812, pyunits.kg / pyunits.mol),
                "cp_mol_liq_comp_coeff": (185.0, pyunits.J / pyunits.mol / pyunits.K),
                "dens_mol_liq_comp_coeff": (6500.0, pyunits.mol / pyunits.m**3),
            },
        },
        "acetic_anhydride": {
            "type": Component,
            "cp_mol_liq_comp": Constant.cp_mol_liq_comp,
            "enth_mol_liq_comp": Constant.enth_mol_liq_comp,
            "dens_mol_liq_comp": Constant.dens_mol_liq_comp,
            "parameter_data": {
                "mw": (0.10209, pyunits.kg / pyunits.mol),
                "cp_mol_liq_comp_coeff": (155.0, pyunits.J / pyunits.mol / pyunits.K),
                "dens_mol_liq_comp_coeff": (10600.0, pyunits.mol / pyunits.m**3),
            },
        },
        "aspirin": {
            "type": Component,
            "cp_mol_liq_comp": Constant.cp_mol_liq_comp,
            "enth_mol_liq_comp": Constant.enth_mol_liq_comp,
            "dens_mol_liq_comp": Constant.dens_mol_liq_comp,
            "parameter_data": {
                "mw": (0.18016, pyunits.kg / pyunits.mol),
                "cp_mol_liq_comp_coeff": (220.0, pyunits.J / pyunits.mol / pyunits.K),
                "dens_mol_liq_comp_coeff": (6700.0, pyunits.mol / pyunits.m**3),
            },
        },
        "acetic_acid": {
            "type": Component,
            "cp_mol_liq_comp": Constant.cp_mol_liq_comp,
            "enth_mol_liq_comp": Constant.enth_mol_liq_comp,
            "dens_mol_liq_comp": Constant.dens_mol_liq_comp,
            "parameter_data": {
                "mw": (0.06005, pyunits.kg / pyunits.mol),
                "cp_mol_liq_comp_coeff": (125.0, pyunits.J / pyunits.mol / pyunits.K),
                "dens_mol_liq_comp_coeff": (17500.0, pyunits.mol / pyunits.m**3),
            },
        },
        "sulfuric_acid": {
            "type": Component,
            "cp_mol_liq_comp": Constant.cp_mol_liq_comp,
            "enth_mol_liq_comp": Constant.enth_mol_liq_comp,
            "dens_mol_liq_comp": Constant.dens_mol_liq_comp,
            "parameter_data": {
                "mw": (0.098079, pyunits.kg / pyunits.mol),
                "cp_mol_liq_comp_coeff": (140.0, pyunits.J / pyunits.mol / pyunits.K),
                "dens_mol_liq_comp_coeff": (18800.0, pyunits.mol / pyunits.m**3),
            },
        },
        "water": {
            "type": Component,
            "cp_mol_liq_comp": Constant.cp_mol_liq_comp,
            "enth_mol_liq_comp": Constant.enth_mol_liq_comp,
            "dens_mol_liq_comp": Constant.dens_mol_liq_comp,
            "parameter_data": {
                "mw": (0.01801528, pyunits.kg / pyunits.mol),
                "cp_mol_liq_comp_coeff": (75.3, pyunits.J / pyunits.mol / pyunits.K),
                "dens_mol_liq_comp_coeff": (55500.0, pyunits.mol / pyunits.m**3),
            },
        },
    },
    "phases": {
        "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},
    },
    "state_definition": FTPx,
    "include_enthalpy_of_formation": False,
    "state_bounds": {
        "flow_mol": (0, 10, 100, pyunits.mol / pyunits.s),
        "temperature": (250, 298.15, 500, pyunits.K),
        "pressure": (1e4, 101325, 1e6, pyunits.Pa),
    },
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    "state_components": StateIndex.true,
}