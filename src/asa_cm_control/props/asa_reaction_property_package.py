# To-Do:
# - Fill out the doc string kwarg pass-ins for everything
# - Check if there is any other reactions that could be added

from idaes.core import (
    declare_process_block_class,
    ReactionParameterBlock,
    ReactionBlockBase,
    ReactionBlockDataBase,
    MaterialFlowBasis,
)
from pyomo.environ import (
    Var,
    Param,
    Set,
    Expression,
    Constraint,
    NonNegativeReals,
    units as pyunits,
    exp,
    Reals,
    PositiveReals,
    value,
)
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog

@declare_process_block_class("ASAReactionParameterBlock")
class ASAReactionParameterData(ReactionParameterBlock):
    
    def build(self):
        
        super().build()
        self._reaction_block_class = ASAReactionBlock  # pyright: ignore[reportUndefinedVariable]
        
        property_package = self.config.property_package
        if property_package is None:
            raise ConfigurationError("ASAReactionParameterBlock requires property_package.")
        
        self.rate_reaction_idx = Set(
            initialize=[
                "r1_aspirin_synthesis",
                "r2_acetic_anhydride_hydrolysis",
                "r3_aspirin_hydrolysis",
            ]
        )
        
        self.rate_reaction_stoichiometry = {}
        for rxn in self.rate_reaction_idx:
            for phase in property_package.phase_list:
                for comp in property_package.component_list:
                    self.rate_reaction_stoichiometry[(rxn, phase, comp)] = 0
        
        self.rate_reaction_stoichiometry[("r1_aspirin_synthesis", "liquid", "salicylic_acid")] = -1
        self.rate_reaction_stoichiometry[("r1_aspirin_synthesis", "liquid", "acetic_anhydride")] = -1
        self.rate_reaction_stoichiometry[("r1_aspirin_synthesis", "liquid", "aspirin")] = 1
        self.rate_reaction_stoichiometry[("r1_aspirin_synthesis", "liquid", "acetic_acid")] = 1
        
        self.rate_reaction_stoichiometry[("r2_acetic_anhydride_hydrolysis", "liquid", "acetic_anhydride")] = -1
        self.rate_reaction_stoichiometry[("r2_acetic_anhydride_hydrolysis", "liquid", "acetic_acid")] = 2
        self.rate_reaction_stoichiometry[("r2_acetic_anhydride_hydrolysis", "liquid", "water")] = -1
        
        self.rate_reaction_stoichiometry[("r3_aspirin_hydrolysis", "liquid", "salicylic_acid")] = 1
        self.rate_reaction_stoichiometry[("r3_aspirin_hydrolysis", "liquid", "aspirin")] = -1
        self.rate_reaction_stoichiometry[("r3_aspirin_hydrolysis", "liquid", "acetic_acid")] = 1
        self.rate_reaction_stoichiometry[("r3_aspirin_hydrolysis", "liquid", "water")] = -1
        
        self.equilibrium_reaction_idx = Set(
            initialize=[
                "e1_h2so4_dissociation",
                "e2_hso4_dissociation",
            ]
        )
        
        self.equilibrium_reaction_stoichiometry = {}
        for reaction in self.equilibrium_reaction_idx:
            for phase in property_package.phase_list:
                for component in property_package.component_list:
                    self.equilibrium_reaction_stoichiometry[(reaction, phase, component)] = 0
        
        self.equilibrium_reaction_stoichiometry[("e1_h2so4_dissociation", "liquid", "H2SO4")] = -1
        self.equilibrium_reaction_stoichiometry[("e1_h2so4_dissociation", "liquid", "H_plus")] = 1
        self.equilibrium_reaction_stoichiometry[("e1_h2so4_dissociation", "liquid", "HSO4_minus")] = 1
        
        self.equilibrium_reaction_stoichiometry[("e2_hso4_dissociation", "liquid", "HSO4_minus")] = -1
        self.equilibrium_reaction_stoichiometry[("e2_hso4_dissociation", "liquid", "H_plus")] = 1
        self.equilibrium_reaction_stoichiometry[("e2_hso4_dissociation", "liquid", "SO4_2minus")] = 1
        
        self.Keq_e1 = Var(
            initialize=1e3,
            domain=PositiveReals,
            units=pyunits.dimensionless,
        )
        
        self.Keq_e1.fix()
        
        self.Keq_e2 = Var(
            initialize=1.2e-2,
            domain=PositiveReals,
            units=pyunits.dimensionless,
        )
        
        self.Keq_e2.fix()
        
        reaction_rate_units = pyunits.mol / pyunits.m**3 / pyunits.s
        
        self.A0_1 = Var(
            initialize=0.0,
            domain=NonNegativeReals,
            units=reaction_rate_units,
            doc="",
        )
        
        self.A0_1.fix()
        
        self.Ea0_1 = Var(
            initialize=0.0,
            domain=NonNegativeReals,
            units=pyunits.J/pyunits.mol,
            doc="",
        )
        
        self.Ea0_1.fix()
        
        self.Acat_1 = Var(
            initialize=1.1e12,
            domain=NonNegativeReals,
            units=reaction_rate_units,
            doc="",
        )
        
        self.Acat_1.fix()
        
        self.Ea_cat_1 = Var(
            initialize=6.44e4,
            domain=NonNegativeReals,
            units=pyunits.J/pyunits.mol,
            doc="",
        )
        
        self.Ea_cat_1.fix()
        
        self.m_1 = Var(
            initialize=0.0,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="",
        )
        
        self.m_1.fix()
        
        self.alpha_1 = Var(
            initialize=1.0,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="",
        )
        
        self.alpha_1.fix()
        
        self.beta_1 = Var(
            initialize=1.0,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="",
        )
        
        self.beta_1.fix()
        
        
        # PARAMETER SET FOR REACTION 2
        
        self.A0_2 = Var(
            initialize=1.5e9,
            domain=NonNegativeReals,
            units=reaction_rate_units,
            doc="",
        )
        
        self.A0_2.fix()
        
        self.Ea0_2 = Var(
            initialize=5.01e4,
            domain=NonNegativeReals,
            units=pyunits.J/pyunits.mol,
            doc="",
        )
        
        self.Ea0_2.fix()
        
        self.Acat_2 = Var(
            initialize=1.0e13,
            domain=NonNegativeReals,
            units=reaction_rate_units,
            doc="",
        )
        
        self.Acat_2.fix()
        
        self.Ea_cat_2 = Var(
            initialize=5.78e4,
            domain=NonNegativeReals,
            units=pyunits.J/pyunits.mol,
            doc="",
        )
        
        self.Ea_cat_2.fix()
        
        self.m_2 = Var(
            initialize=1.0,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="",
        )
        
        self.m_2.fix()
        
        self.alpha_2 = Var(
            initialize=1.0,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="",
        )
        
        self.alpha_2.fix()
        
        self.beta_2 = Var(
            initialize=1.0,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="",
        )
        
        self.beta_2.fix()
        
        
        # PARAMETER SET FOR REACTION 3
        
        self.A0_3 = Var(
            initialize=1.5e9,
            domain=NonNegativeReals,
            units=reaction_rate_units,
            doc="",
        )
        
        self.A0_3.fix()
        
        self.Ea0_3 = Var(
            initialize=5.01e4,
            domain=NonNegativeReals,
            units=pyunits.J/pyunits.mol,
            doc="",
        )
        
        self.Ea0_3.fix()
        
        self.Acat_3 = Var(
            initialize=1.0e13,
            domain=NonNegativeReals,
            units=reaction_rate_units,
            doc="",
        )
        
        self.Acat_3.fix()
        
        self.Ea_cat_3 = Var(
            initialize=5.78e4,
            domain=NonNegativeReals,
            units=pyunits.J/pyunits.mol,
            doc="",
        )
        
        self.Ea_cat_3.fix()
        
        self.m_3 = Var(
            initialize=1.0,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="",
        )
        
        self.m_3.fix()
        
        self.alpha_3 = Var(
            initialize=1.0,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="",
        )
        
        self.alpha_3.fix()
        
        self.beta_3 = Var(
            initialize=1.0,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="",
        )
        
        self.beta_3.fix()
    
    @classmethod
    def define_metadata(cls, obj):
        
        obj.add_properties(
            {
                'reaction_rate': {'method': '_reaction_rate'},
                'equilibrium_constraint': {'method': '_equilibrium_constraint'},
            }
        )
        
        obj.add_default_units(
            {
                'time': pyunits.s,
                'length': pyunits.m,
                'mass': pyunits.kg,
                'amount': pyunits.mol,
                'temperature': pyunits.K,
            }
        )


# REACTION BLOCK CLASSES

class ASAReactionBlockMethods(ReactionBlockBase):
    
    def initialize(
        self,
        state_vars_fixed=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="reactions")
        init_log.info("Reaction initialization complete (no solve required).")
        
        _ = state_vars_fixed
        _ = solver
        _ = optarg
        
        return None

@declare_process_block_class("ASAReactionBlock", block_class=ASAReactionBlockMethods)
class ASAReactionBlockData(ReactionBlockDataBase):
    
    def build(self):
        super().build()
        if self.config.has_equilibrium:
            self._equilibrium_constraint()
    
    def get_reaction_rate_basis(self):
        return MaterialFlowBasis.molar
    
    def _reaction_rate(self):
        state = self.state_ref
        params = self.params
        eps = 1e-8
        R = 8.314462618 * pyunits.J / pyunits.mol / pyunits.K
        
        a_sa  = state.activity_true_comp["salicylic_acid"] + eps
        a_aa  = state.activity_true_comp["acetic_anhydride"] + eps
        a_h2o = state.activity_true_comp["water"] + eps
        a_asa = state.activity_true_comp["aspirin"] + eps
        a_hplus = state.activity_true_comp["H_plus"] + eps
        
        def reaction_rule(b, reaction):
            
            if reaction == "r1_aspirin_synthesis":
                
                k0 = params.A0_1 * exp(-params.Ea0_1 / (R * state.temperature))
                kcat = params.Acat_1 * exp(-params.Ea_cat_1 / (R * state.temperature))
                
                return (k0 + kcat * a_hplus ** params.m_1) * a_sa ** params.alpha_1 * a_aa ** params.beta_1
            
            if reaction == "r2_acetic_anhydride_hydrolysis":
                
                k0 = params.A0_2 * exp(-params.Ea0_2 / (R * state.temperature))
                kcat = params.Acat_2 * exp(-params.Ea_cat_2 / (R * state.temperature))
                
                return (k0 + kcat * a_hplus ** params.m_2) * a_aa ** params.alpha_2 * a_h2o ** params.beta_2
            
            if reaction == "r3_aspirin_hydrolysis":
                
                k0 = params.A0_3 * exp(-params.Ea0_3 / (R * state.temperature))
                kcat = params.Acat_3 * exp(-params.Ea_cat_3 / (R * state.temperature))
                
                return (k0 + kcat * a_hplus ** params.m_3) * a_asa ** params.alpha_3 * a_h2o ** params.beta_3
            
            raise ConfigurationError(
                f"Unknown reaction '{reaction}' encountered in "
                f"ASAReactionBlockData._reaction_rate for block {self.name}."
            )
        
        self.reaction_rate = Expression(params.rate_reaction_idx, rule=reaction_rule)
    
    def _equilibrium_constraint(self):
        state = self.state_ref
        epsilon = 1e-8
        
        a_h2so4 = state.activity_true_comp["H2SO4"] + epsilon
        a_hplus = state.activity_true_comp["H_plus"] + epsilon
        a_hso4 = state.activity_true_comp["HSO4_minus"] + epsilon
        a_so4 = state.activity_true_comp["SO4_2minus"] + epsilon
        
        def equilibrium_rule(block, reaction):
            if reaction == "e1_h2so4_dissociation":
                return a_hplus * a_hso4 == block.params.Keq_e1 * a_h2so4
            if reaction == "e2_hso4_dissociation":
                return a_hplus * a_so4 == block.params.Keq_e2 * a_hso4
            raise ConfigurationError(f"Unknown equilibrium reaction {reaction}")
        
        self.equilibrium_constraint = Constraint(
            self.params.equilibrium_reaction_idx,
            rule=equilibrium_rule
        )