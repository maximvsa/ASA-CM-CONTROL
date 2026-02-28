from pyomo.environ import units as pyunits
from idaes.models.properties.modular_properties.base.generic_reaction import ConcentrationForm
from idaes.models.properties.modular_properties.reactions.rate_constant import arrhenius
from idaes.models.properties.modular_properties.reactions.rate_forms import power_law_rate
from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn

configuration = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "rate_reactions": {
        "R_sa_to_asa": {
            "stoichiometry": {
                ("Liq", "salicylic_acid"): -1,
                ("Liq", "acetic_anhydride"): -1,
                ("Liq", "aspirin"): 1,
                ("Liq", "acetic_acid"): 1,
                ("Liq", "sulfuric_acid"): 0,
                ("Liq", "water"): 0,
            },
            "heat_of_reaction": constant_dh_rxn,
            "rate_constant": arrhenius,
            "rate_form": power_law_rate,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (-4.5e4, pyunits.J / pyunits.mol),
                "arrhenius_const": (5.0e4, pyunits.mol / pyunits.m**3 / pyunits.s),
                "energy_activation": (3.5e4, pyunits.J / pyunits.mol),
            },
        }
    },
}
