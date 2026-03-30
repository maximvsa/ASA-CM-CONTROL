# To-Do:
# - Fill out the doc string kwarg pass-ins for everything


"""Custom reaction property package for ASA synthesis kinetics.

This module defines an IDAES reaction package with:
- A reaction parameter block containing stoichiometry and fixed kinetic constants.
- A reaction state block that builds reaction-rate expressions on demand.

Current scope is liquid-phase kinetics for aspirin synthesis and acetic
anhydride hydrolysis.
"""

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


# PARAMETER BLOCK CLASS

@declare_process_block_class("ASAReactionParameterBlock")
class ASAReactionParameterData(ReactionParameterBlock):
    """Reaction parameter block for ASA process chemistry.

    The block stores reaction indices, stoichiometric maps, and fixed kinetic
    parameters used by associated reaction state blocks.
    """

    def build(self):
        """Construct reaction sets, stoichiometry, and kinetic parameters.

        Defines two liquid-phase rate reactions:
        - r1_aspirin_synthesis
        - r2_acetic_anhydride_hydrolysis

        Kinetic constants are stored as fixed Vars for straightforward
        calibration updates and reproducible simulation behavior.
        """
        
        super().build()
        self._reaction_block_class = ASAReactionBlock
        
        property_package = self.config.property_package
        if property_package is None:
            raise ConfigurationError("ASAReactionParameterBlock requires property_package.")
        
        self.rate_reaction_idx = Set(initialize=[
            "r1_aspirin_synthesis",
            "r2_acetic_anhydride_hydrolysis",
            "r3_aspirin_hydrolysis"
            ]
        )
        
        # Empty for now, might be filled in later if any reactions get upgraded to equilibrium
        self.equilibrium_reaction_idx = Set(initialize=[])
        
        self.rate_reaction_stoichiometry = {
            # r1: salicylic_acid + acetic_anhydride -> aspirin + acetic_acid (liquid only)
            ("r1_aspirin_synthesis", "liquid", "salicylic_acid"): -1,
            ("r1_aspirin_synthesis", "liquid", "acetic_anhydride"): -1,
            ("r1_aspirin_synthesis", "liquid", "sulfuric_acid"): 0,
            ("r1_aspirin_synthesis", "liquid", "aspirin"): 1,
            ("r1_aspirin_synthesis", "liquid", "acetic_acid"): 1,
            ("r1_aspirin_synthesis", "liquid", "water"): 0,
            ("r1_aspirin_synthesis", "vapor", "salicylic_acid"): 0,
            ("r1_aspirin_synthesis", "vapor", "acetic_anhydride"): 0,
            ("r1_aspirin_synthesis", "vapor", "sulfuric_acid"): 0,
            ("r1_aspirin_synthesis", "vapor", "aspirin"): 0,
            ("r1_aspirin_synthesis", "vapor", "acetic_acid"): 0,
            ("r1_aspirin_synthesis", "vapor", "water"): 0,
            ("r1_aspirin_synthesis", "solid", "salicylic_acid"): 0,
            ("r1_aspirin_synthesis", "solid", "acetic_anhydride"): 0,
            ("r1_aspirin_synthesis", "solid", "sulfuric_acid"): 0,
            ("r1_aspirin_synthesis", "solid", "aspirin"): 0,
            ("r1_aspirin_synthesis", "solid", "acetic_acid"): 0,
            ("r1_aspirin_synthesis", "solid", "water"): 0,
            
            # r2: acetic_anhydride + water -> 2 acetic_acid (liquid only)
            ("r2_acetic_anhydride_hydrolysis", "liquid", "salicylic_acid"): 0,
            ("r2_acetic_anhydride_hydrolysis", "liquid", "acetic_anhydride"): -1,
            ("r2_acetic_anhydride_hydrolysis", "liquid", "sulfuric_acid"): 0,
            ("r2_acetic_anhydride_hydrolysis", "liquid", "aspirin"): 0,
            ("r2_acetic_anhydride_hydrolysis", "liquid", "acetic_acid"): 2,
            ("r2_acetic_anhydride_hydrolysis", "liquid", "water"): -1,
            ("r2_acetic_anhydride_hydrolysis", "vapor", "salicylic_acid"): 0,
            ("r2_acetic_anhydride_hydrolysis", "vapor", "acetic_anhydride"): 0,
            ("r2_acetic_anhydride_hydrolysis", "vapor", "sulfuric_acid"): 0,
            ("r2_acetic_anhydride_hydrolysis", "vapor", "aspirin"): 0,
            ("r2_acetic_anhydride_hydrolysis", "vapor", "acetic_acid"): 0,
            ("r2_acetic_anhydride_hydrolysis", "vapor", "water"): 0,
            ("r2_acetic_anhydride_hydrolysis", "solid", "salicylic_acid"): 0,
            ("r2_acetic_anhydride_hydrolysis", "solid", "acetic_anhydride"): 0,
            ("r2_acetic_anhydride_hydrolysis", "solid", "sulfuric_acid"): 0,
            ("r2_acetic_anhydride_hydrolysis", "solid", "aspirin"): 0,
            ("r2_acetic_anhydride_hydrolysis", "solid", "acetic_acid"): 0,
            ("r2_acetic_anhydride_hydrolysis", "solid", "water"): 0,
            
            # r3: aspirin + water -> salicylic_acid + acetic_acid (liquid_only)
            ("r3_aspirin_hydrolysis", "liquid", "salicylic_acid"): 1,
            ("r3_aspirin_hydrolysis", "liquid", "acetic_anhydride"): 0,
            ("r3_aspirin_hydrolysis", "liquid", "sulfuric_acid"): 0,
            ("r3_aspirin_hydrolysis", "liquid", "aspirin"): -1,
            ("r3_aspirin_hydrolysis", "liquid", "acetic_acid"): 1,
            ("r3_aspirin_hydrolysis", "liquid", "water"): -1,
            ("r3_aspirin_hydrolysis", "vapor", "salicylic_acid"): 0,
            ("r3_aspirin_hydrolysis", "vapor", "acetic_anhydride"): 0,
            ("r3_aspirin_hydrolysis", "vapor", "sulfuric_acid"): 0,
            ("r3_aspirin_hydrolysis", "vapor", "aspirin"): 0,
            ("r3_aspirin_hydrolysis", "vapor", "acetic_acid"): 0,
            ("r3_aspirin_hydrolysis", "vapor", "water"): 0,
            ("r3_aspirin_hydrolysis", "solid", "salicylic_acid"): 0,
            ("r3_aspirin_hydrolysis", "solid", "acetic_anhydride"): 0,
            ("r3_aspirin_hydrolysis", "solid", "sulfuric_acid"): 0,
            ("r3_aspirin_hydrolysis", "solid", "aspirin"): 0,
            ("r3_aspirin_hydrolysis", "solid", "acetic_acid"): 0,
            ("r3_aspirin_hydrolysis", "solid", "water"): 0,
        }
        
        # EMPTY FOR NOW, WILL FILL IN IF EQUILIBRIUM REACTIONS ARE EVER USED
        self.equilibrium_reaction_stoichiometry = {}
        
        reaction_rate_units = pyunits.mol / pyunits.m**3 / pyunits.s
        
        
        # PARAMETER SET FOR REACTION 1
        
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
        """Declare supported reaction properties and default units.

        Args:
            cls: Reaction parameter block class reference.
            obj: IDAES metadata object to populate.
        """
        
        obj.add_properties(
            {
                'reaction_rate': {'method': '_reaction_rate'},
            }
        )
        
        obj.add_default_units(
            {
                'time': pyunits.s,
                'length': pyunits.m,
                'mass': pyunits.kg,
                'amount': pyunits.mol,
                'temperature': pyunits.K
            }
        )


# REACTION BLOCK CLASSES

class ASAReactionBlockMethods(ReactionBlockBase):
    """Mixin methods attached to the generated ASAReactionBlock class."""
    
    def initialize(
        self,
        state_vars_fixed=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """Initialize reaction block data and log completion.

        Args:
            state_vars_fixed: Placeholder for API compatibility.
            outlvl: IDAES logging level.
            solver: Placeholder for API compatibility.
            optarg: Placeholder for API compatibility.

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="reactions")
        init_log.info("Reaction initialization complete (no solve required).")
        
        _ = state_vars_fixed
        _ = solver
        _ = optarg
        
        return None

@declare_process_block_class("ASAReactionBlock", block_class=ASAReactionBlockMethods)
class ASAReactionBlockData(ReactionBlockDataBase):
    """Reaction state block data class for ASA reaction-rate calculations."""
    
    def build(self):
        """Construct the reaction block data object."""
        
        super().build()
    
    def get_reaction_rate_basis(self):
        """Return the reaction-rate basis.

        Returns:
            MaterialFlowBasis: MaterialFlowBasis.molar.
        """
        return MaterialFlowBasis.molar
    
    def _reaction_rate(self):
        """Build reaction-rate expressions for the configured rate reactions.

        Activity approximation:
            a_i = gamma_i x_i

        Arrhenius terms:
            k_0 = A_0 exp(-E_{a,0} / (R T))
            k_cat = A_cat exp(-E_{a,cat} / (R T))

        Rate forms:
            r_1 = (k_0 + k_cat a_H+^{m_1}) a_SA^{alpha_1} a_AA^{beta_1}
            r_2 = (k_0 + k_cat a_H+^{m_2}) a_AA^{alpha_2} a_H2O^{beta_2}

        where SA is salicylic acid and AA is acetic anhydride.
        """
        
        state = self.state_ref
        params = self.params
        gamma = state.act_coeff_liq_comp
        eps = 1e-12
        R = 8.314462618 * pyunits.J / pyunits.mol / pyunits.K
        
        a_hplus = gamma["sulfuric_acid"] * state.mole_frac_comp["sulfuric_acid"] + eps
        a_sa = gamma["salicylic_acid"] * state.mole_frac_comp["salicylic_acid"] + eps
        a_aa = gamma["acetic_anhydride"] * state.mole_frac_comp["acetic_anhydride"] + eps
        a_h2o = gamma["water"] * state.mole_frac_comp["water"] + eps
        a_asa = gamma["aspirin"] * state.mole_frac_comp["aspirin"] + eps
        
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