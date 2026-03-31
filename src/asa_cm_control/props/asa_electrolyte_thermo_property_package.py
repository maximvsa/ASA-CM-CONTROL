from idaes.core import (
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    Component,
    LiquidPhase,
    VaporPhase,
    SolidPhase,
    MaterialFlowBasis,
    MaterialBalanceType,
    EnergyBalanceType,
)
from pyomo.environ import (Var,
    Expression,
    Constraint,
    NonNegativeReals,
    units as pyunits,
    Reals,
    Set,
    exp,
    log,
)
from idaes.core.util.initialization import fix_state_vars, revert_state_vars
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog

@declare_process_block_class("ASAElectroThermoParameterBlock")
class ASAElectroThermoParameterData(PhysicalParameterBlock):
    def build(self):
        super().build()
        self._state_block_class = ASAElectroThermoStateBlock
        
        self.salicylic_acid = Component()
        self.acetic_anhydride = Component()
        self.aspirin = Component()
        self.acetic_acid = Component()
        self.water = Component()
        
        self.H2SO4 = Component()
        self.H_plus = Component()
        self.HSO4_minus = Component()
        self.SO4_2minus = Component()
        
        self.liquid = LiquidPhase()
        
        self.pressure_ref = Var(
            initialize=101325.0,
            units=pyunits.Pa,
            doc="Reference pressure",
        )
        
        self.pressure_ref.fix()
        
        self.temperature_ref = Var(
            initialize=298.15,
            units=pyunits.K,
            doc="Reference temperature",
        )
        
        self.temperature_ref.fix()
        
        self.mw_comp = Var(
            self.component_list,
            initialize={
                "salicylic_acid": 0.13812,
                "acetic_anhydride": 0.10209,
                "aspirin": 0.18016,
                "acetic_acid": 0.06005,
                "water": 0.01801528,
                "H2SO4": 0,
                "H_plus": 0,
                "HSO4_minus": 0,
                "SO4_2minus": 0,
            },
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weight of each component in kg/mol",
        )
        
        self.mw_comp.fix()
        
        self.cp_mol_liq_comp = Var(
            self.component_list,
            initialize={
                "salicylic_acid": 210.0,
                "acetic_anhydride": 168.2,
                "aspirin": 270.0,
                "acetic_acid": 123.1,
                "water": 75.37,
                "H2SO4": 0,
                "H_plus": 0,
                "HSO4_minus": 0,
                "SO4_2minus": 0,
            },
            units=pyunits.J / pyunits.mol / pyunits.K,
            doc="Constant molar heat capacity for liquid phase by component"
        )
        
        self.cp_mol_liq_comp.fix()
        
        self.density_liq_comp = Var(
            self.component_list,
            initialize={
                "salicylic_acid": 1440.0,
                "acetic_anhydride": 1060.0,
                "aspirin": 1330.0,
                "acetic_acid": 1053.0,
                "water": 1000.0,
                "H2SO4": 0,
                "H_plus": 0,
                "HSO4_minus": 0,
                "SO4_2minus": 0,
            },
            units=pyunits.kg / pyunits.m**3,
            doc="Density of pure liquid component at reference conditions (Some are estimated based on literature values at 25C, may need to be updated with more accurate values or temperature dependence)"
        )
        
        self.density_liq_comp.fix()
        
        self.dh_form_liq_comp = Var(
            self.component_list,
            initialize={
                "salicylic_acid": -585260,
                "acetic_anhydride": -625000,
                "aspirin": -758200,
                "acetic_acid": -484000,
                "water": -286000,
                "H2SO4": 0,
                "H_plus": 0,
                "HSO4_minus": 0,
                "SO4_2minus": 0,
            },
            units=pyunits.J / pyunits.mol,
            doc="Standard enthalpy of formation at reference conditions of 298.15K for liquid phase by component"
        )
        
        self.dh_form_liq_comp.fix()
        
        self.charge_comp = Var(
            self.component_list,
            initialize={
                "salicylic_acid": 0,
                "acetic_anhydride": 0,
                "aspirin": 0,
                "acetic_acid": 0,
                "water": 0,
                "H2SO4": 0,
                "H_plus": 1,
                "HSO4_minus": -1,
                "SO4_2minus": -2,
            },
            units=pyunits.dimensionless,
            doc=""
        )
        
        self.charge_comp.fix()
        
        self.Ka1 = Var(
            self.component_list,
            initialize={
                "salicylic_acid": 0,
                "acetic_anhydride": 0,
                "aspirin": 0,
                "acetic_acid": 0,
                "water": 0,
                "H2SO4": 0,
                "H_plus": 0,
                "HSO4_minus": 0,
                "SO4_2minus": 0,
            },
            units=pyunits.dimensionless,
            doc=""
        )
        
        self.Ka1.fix()
        
        self.Ka2 = Var(
            self.component_list,
            initialize={
                "salicylic_acid": 0,
                "acetic_anhydride": 0,
                "aspirin": 0,
                "acetic_acid": 0,
                "water": 0,
                "H2SO4": 0,
                "H_plus": 0,
                "HSO4_minus": 0,
                "SO4_2minus": 0,
            },
            units=pyunits.dimensionless,
            doc=""
        )
        
        self.Ka2.fix()