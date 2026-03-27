# To-Do:
# - Get the reaction package filled out to the style of the thermo package (differing where it matters)
# - Make sure everything listed in the bullet list below The Reaction Parameter Block section of the documentation gets put in the class
# - Replace placeholder values in global reaction parameters with actual estimates
# - Fill in the doc strings on the global parameters declared in the ReactionParameterBlock class
# - Load the script up with doc strings for each class and function



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

    def build(self):
        
        super().build()
        self._reaction_block_class = ASAReactionBlock
        
        property_package = self.config.property_package
        self._property_package = property_package
        
        self.rate_reaction_idx = Set(initialize=[
            "r1_aspirin_synthesis",
            "r2_acetic_anhydride_hydrolysis",
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
        }
        
        reaction_rate_units = pyunits.mol / pyunits.m**3 / pyunits.s
        
        
        # PARAMETER SET FOR REACTION 1
        
        self.A0_1 = Var(
            initialize=1e5,
            domain=PositiveReals,
            units=reaction_rate_units,
            doc="",
        )
        
        self.A0_1.fix()
        
        self.Ea0_1 = Var(
            initialize=5e4,
            domain=NonNegativeReals,
            units=pyunits.J/pyunits.mol,
            doc="",
        )
        
        self.Ea0_1.fix()
        
        self.Acat_1 = Var(
            initialize=1e5,
            domain=PositiveReals,
            units=reaction_rate_units,
            doc="",
        )
        
        self.Acat_1.fix()
        
        self.Ea_cat_1 = Var(
            initialize=5e4,
            domain=NonNegativeReals,
            units=pyunits.J/pyunits.mol,
            doc="",
        )
        
        self.Ea_cat_1.fix()
        
        self.m_1 = Var(
            initialize=1.0,
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
            initialize=1e5,
            domain=PositiveReals,
            units=reaction_rate_units,
            doc="",
        )
        
        self.A0_2.fix()
        
        self.Ea0_2 = Var(
            initialize=5e4,
            domain=NonNegativeReals,
            units=pyunits.J/pyunits.mol,
            doc="",
        )
        
        self.Ea0_2.fix()
        
        self.Acat_2 = Var(
            initialize=1e5,
            domain=PositiveReals,
            units=reaction_rate_units,
            doc="",
        )
        
        self.Acat_2.fix()
        
        self.Ea_cat_2 = Var(
            initialize=5e4,
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

    @classmethod
    def define_metadata(cls, obj):
        
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

