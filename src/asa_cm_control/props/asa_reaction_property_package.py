# To-Do:
# - Get the reaction package filled out to the style of the thermo package (differing where it matters)
# - Make sure everything listed in the bullet list below The Reaction Parameter Block section of the documentation gets put in the class
# - 

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
        
        if self.config.property_package is None:
            raise ConfigurationError("ASAReactionParameterBlock requires property_package.")
        else:
            property_package = self.config.property_package
        # Use property_package.component_list, property_package.phase_list, property_package.get_metadata(), etc.
        # for validation and setup of stoichiometry/units compatibility.

    @classmethod
    def define_metadata(cls, obj):
        
        obj.add_properties(
            {
                
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