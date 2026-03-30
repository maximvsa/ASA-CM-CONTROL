# To-Do:
# - 


"""Custom liquid-focused thermophysical package for first-pass ASA synthesis modeling.

This module defines a manual IDAES property package with:
- A parameter block containing component/phase definitions and fixed global constants.
- A state block with required interface methods for material and energy balances.
- Build-on-demand property expressions used by downstream unit models.

The current implementation assumes a liquid-dominant approximation for derived
phase quantities.
"""

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
)
from idaes.core.util.initialization import fix_state_vars, revert_state_vars
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog


# PARAMETER BLOCK CLASS

@declare_process_block_class("ThermoParameterBlock")
class ThermoParameterData(PhysicalParameterBlock):
    """Thermophysical parameter block for the ASA component set.

    The block defines components, phases, reference conditions, and fixed global
    property constants used by associated state blocks.
    """

    def build(self):
        """Construct the parameter block and fixed global property variables.

        This method registers components and phases, defines reference pressure
        and temperature, and creates fixed per-component constants such as
        molecular weights, heat capacities, densities, and formation enthalpies.
        """
        super().build()
        self._state_block_class = ThermoStateBlock
        
        self.salicylic_acid = Component()
        self.acetic_anhydride = Component()
        self.sulfuric_acid = Component()
        self.aspirin = Component()
        self.acetic_acid = Component()
        self.water = Component()
        
        self.liquid = LiquidPhase()
        self.vapor = VaporPhase()
        self.solid = SolidPhase()
        
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
                "sulfuric_acid": 0.098079,
                "aspirin": 0.18016,
                "acetic_acid": 0.06005,
                "water": 0.01801528,
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
                "sulfuric_acid": 137.6,
                "aspirin": 270.0,
                "acetic_acid": 123.1,
                "water": 75.37,
            },
            units=pyunits.J / pyunits.mol / pyunits.K,
            doc="Constant molar heat capacity for liquid phase by component"
        )
        
        self.cp_mol_liq_comp.fix()
        
        self.cp_mol_vap_comp = Var(
            self.component_list,
            initialize={
                "salicylic_acid": 120.0,
                "acetic_anhydride": 99.50,
                "sulfuric_acid": 83.71,
                "aspirin": 149.0,
                "acetic_acid": 63.44,
                "water": 33.58,
            },
            units=pyunits.J/pyunits.mol/pyunits.K,
            doc="Constant molar heat capacity for vapor phase by component"
        )
        
        self.cp_mol_vap_comp.fix()
        
        self.cp_mol_sol_comp = Var(
            self.component_list,
            initialize={
                "salicylic_acid": 160.9,
                "acetic_anhydride": 168.2,
                "sulfuric_acid": 137.6,
                "aspirin": 217.8,
                "acetic_acid": 123.1,
                "water": 37.77,
            },
            units=pyunits.J/pyunits.mol/pyunits.K,
            doc="Constant molar heat capacity for solid phase by component"
        )
        
        self.cp_mol_sol_comp.fix()
        
        self.density_liq_comp = Var(
            self.component_list,
            initialize={
                "salicylic_acid": 1440.0,
                "acetic_anhydride": 1060.0,
                "sulfuric_acid": 1840.0,
                "aspirin": 1330.0,
                "acetic_acid": 1053.0,
                "water": 1000.0,
            },
            units=pyunits.kg / pyunits.m**3,
            doc="Density of pure liquid component at reference conditions (Some are estimated based on literature values at 25C, may need to be updated with more accurate values or temperature dependence)"
        )
        
        self.density_liq_comp.fix()
        
        self.density_sol_comp = Var(
            self.component_list,
            initialize={
                "salicylic_acid": 1443.0,
                "acetic_anhydride": 1346.0,
                "sulfuric_acid": 1940.0,
                "aspirin": 1335.0,
                "acetic_acid": 1266.0,
                "water": 916.7,
            },
            units=pyunits.kg / pyunits.m**3,
            doc="Density of pure solid component at reference conditions (Some are estimated based on literature values at 25C, may need to be updated with more accurate values or temperature dependence)"
        )
        
        self.density_sol_comp.fix()
        
        # STANDARD ENTHALPY OF FORMATION (J/mol)
        self.dh_form_liq_comp = Var(
            self.component_list,
            initialize={
                "salicylic_acid": -585260,
                "acetic_anhydride": -625000,
                "sulfuric_acid": -814000,
                "aspirin": -758200,
                "acetic_acid": -484000,
                "water": -286000,
            },
            units=pyunits.J / pyunits.mol,
            doc="Standard enthalpy of formation at reference conditions of 298.15K for liquid phase by component"
        )
        
        self.dh_form_liq_comp.fix()
        
        # NRTL Parameters:
        
        self.nrtl_pair_set = Set(
            initialize=[(i, j) for i in self.component_list for j in self.component_list],
            dimen=2,
        )
        
        self.tau_nrtl = Var(
            self.nrtl_pair_set,
            initialize=0.0,
            domain=Reals,
            units=pyunits.dimensionless
        )
        
        self.tau_nrtl.fix()
        
        self.alpha_nrtl = Var(
            self.nrtl_pair_set,
            initialize=0.3,
            domain=NonNegativeReals,
            units=pyunits.dimensionless
        )
        
        self.alpha_nrtl.fix()
        
        for c in self.component_list:
            self.tau_nrtl[c, c].fix(0.0)
            self.alpha_nrtl[c, c].fix(0.0)
        
        # -----------------------------
        # tau_ij seeds (dimensionless)
        # -----------------------------
        
        self.tau_nrtl["salicylic_acid", "acetic_anhydride"].fix(0.35)
        self.tau_nrtl["acetic_anhydride", "salicylic_acid"].fix(0.15)
        
        self.tau_nrtl["salicylic_acid", "sulfuric_acid"].fix(0.10)
        self.tau_nrtl["sulfuric_acid", "salicylic_acid"].fix(0.80)
        
        self.tau_nrtl["salicylic_acid", "aspirin"].fix(0.20)
        self.tau_nrtl["aspirin", "salicylic_acid"].fix(0.15)
        
        self.tau_nrtl["salicylic_acid", "acetic_acid"].fix(0.55)
        self.tau_nrtl["acetic_acid", "salicylic_acid"].fix(0.30)
        
        self.tau_nrtl["salicylic_acid", "water"].fix(1.80)
        self.tau_nrtl["water", "salicylic_acid"].fix(0.60)
        
        self.tau_nrtl["acetic_anhydride", "sulfuric_acid"].fix(0.20)
        self.tau_nrtl["sulfuric_acid", "acetic_anhydride"].fix(0.90)
        
        self.tau_nrtl["acetic_anhydride", "aspirin"].fix(0.45)
        self.tau_nrtl["aspirin", "acetic_anhydride"].fix(0.20)
        
        self.tau_nrtl["acetic_anhydride", "acetic_acid"].fix(0.25)
        self.tau_nrtl["acetic_acid", "acetic_anhydride"].fix(0.10)
        
        self.tau_nrtl["acetic_anhydride", "water"].fix(2.20)
        self.tau_nrtl["water", "acetic_anhydride"].fix(0.40)
        
        self.tau_nrtl["sulfuric_acid", "aspirin"].fix(0.20)
        self.tau_nrtl["aspirin", "sulfuric_acid"].fix(0.90)
        
        self.tau_nrtl["sulfuric_acid", "acetic_acid"].fix(-0.20)
        self.tau_nrtl["acetic_acid", "sulfuric_acid"].fix(0.80)
        
        self.tau_nrtl["sulfuric_acid", "water"].fix(-0.80)
        self.tau_nrtl["water", "sulfuric_acid"].fix(2.50)
        
        self.tau_nrtl["aspirin", "acetic_acid"].fix(0.40)
        self.tau_nrtl["acetic_acid", "aspirin"].fix(0.25)
        
        self.tau_nrtl["aspirin", "water"].fix(2.00)
        self.tau_nrtl["water", "aspirin"].fix(0.80)
        
        self.tau_nrtl["acetic_acid", "water"].fix(1.20)
        self.tau_nrtl["water", "acetic_acid"].fix(-0.35)
        
        # -----------------------------------
        # alpha_ij seeds (dimensionless)
        # -----------------------------------
        # Keep most pairs at 0.30, override selected pairs below
        
        self.alpha_nrtl["salicylic_acid", "water"].fix(0.35)
        self.alpha_nrtl["water", "salicylic_acid"].fix(0.35)
        
        self.alpha_nrtl["aspirin", "water"].fix(0.35)
        self.alpha_nrtl["water", "aspirin"].fix(0.35)
        
        self.alpha_nrtl["sulfuric_acid", "water"].fix(0.20)
        self.alpha_nrtl["water", "sulfuric_acid"].fix(0.20)
        
        self.alpha_nrtl["sulfuric_acid", "acetic_acid"].fix(0.20)
        self.alpha_nrtl["acetic_acid", "sulfuric_acid"].fix(0.20)
        
        self.alpha_nrtl["sulfuric_acid", "acetic_anhydride"].fix(0.20)
        self.alpha_nrtl["acetic_anhydride", "sulfuric_acid"].fix(0.20)
        
        self.alpha_nrtl["sulfuric_acid", "salicylic_acid"].fix(0.20)
        self.alpha_nrtl["salicylic_acid", "sulfuric_acid"].fix(0.20)
        
        self.alpha_nrtl["sulfuric_acid", "aspirin"].fix(0.20)
        self.alpha_nrtl["aspirin", "sulfuric_acid"].fix(0.20)
        
        

    @classmethod
    def define_metadata(cls, obj):
        """Declare package metadata for supported properties and default units.

        Args:
            cls: Parameter block class reference.
            obj: IDAES metadata object to populate.
        """
        obj.add_properties(
            {
                'flow_mol': {'method': None},
                'temperature': {'method': None},
                'pressure': {'method': None},
                'mole_frac_comp': {'method': None},
                'enth_mol': {'method': '_enth_mol'},
                'dens_mass': {'method': '_dens_mass'},
                'cp_mol': {'method': '_cp_mol'},
                'flow_mol_phase_comp': {'method': '_flow_mol_phase_comp'},
                'mole_frac_phase_comp': {'method': '_mole_frac_phase_comp'},
                'phase_frac': {'method': '_phase_frac'},
            }
        )
        
        obj.define_custom_properties(
            {
                'act_coeff_liq_comp': {'method': '_act_coeff_liq_comp'},
                'log_gamma_liq_comp': {'method': None},
                'G_nrtl': {'method': None},
                'S_nrtl': {'method': None},
                'N_nrtl': {'method': None},
                'Q_nrtl': {'method': None},
                'P_nrtl': {'method': None},
                'W_nrtl': {'method': None},
                'D_nrtl': {'method': None},
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


# STATE BLOCK CLASSES

class ThermoStateBlockMethods(StateBlock):
    """Mixin methods attached to the generated ThermoStateBlock class.

    These methods provide initialization and state-release utilities required by
    IDAES workflows.
    """

    def initialize(
        self,
        state_args=None,
        state_vars_fixed=False,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """Initialize state variables for one or more state block instances.

        Args:
            state_args: Optional dictionary of state variable initial guesses.
            state_vars_fixed: If True, assumes the caller already fixed state vars.
            hold_state: If True, keep variables fixed and return release flags.
            outlvl: IDAES logging level.
            solver: Unused placeholder for API compatibility.
            optarg: Unused placeholder for API compatibility.

        Returns:
            dict or None: A flags dictionary if hold_state is True, otherwise None.
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        
        # Fix state variables unless caller already did so
        flags = None
        if not state_vars_fixed:
            flags = fix_state_vars(self, state_args)
        else:
            flags = {}
        
        # No property solve needed yet: all derived properties are expressions
        init_log.info("Property initialization complete (no solve required).")
        
        if hold_state:
            return flags
        
        # Release state vars if we fixed them here
        if not state_vars_fixed:
            revert_state_vars(self, flags)
        
        _ = solver
        _ = optarg
        
        return None
    
    
    def fix_initialization_states(self):
        """Fix all state variables on all indexed state block members."""
        fix_state_vars(self)
    
    
    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        """Release state variables fixed during initialize with hold_state=True.

        Args:
            flags: Flags returned from initialize(..., hold_state=True).
            outlvl: IDAES logging level.

        Raises:
            ConfigurationError: If flags is None.
        """
        
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        
        if flags is None:
            raise ConfigurationError(
                "release_state called with flags=None. "
                "Pass flags returned by initialize(..., hold_state=True)."
            )
        
        revert_state_vars(self, flags)
        init_log.info("State variables released.")



@declare_process_block_class("ThermoStateBlock", block_class=ThermoStateBlockMethods)
class ThermoStateBlockData(StateBlockData):
    """State block data class for the custom ASA thermophysical package.

    This class defines state variables and build-on-demand property expressions
    used in material and energy balances under a liquid-focused approximation.
    """
    
    def build(self):
        """Construct state variables and optional mole-fraction closure equation.

        Creates flow, temperature, pressure, and composition state variables.
        When defined_state is False, adds a mole-fraction summation constraint.
        """
        super().build()
        
        # State variables
        self.flow_mol = Var(
            initialize=1.0,
            domain=NonNegativeReals,
            units=pyunits.mol / pyunits.s,
            doc="Total molar flowrate",
        )
        
        self.temperature = Var(
            initialize=298.15,
            domain=NonNegativeReals,
            units=pyunits.K,
            doc="State temperature",
        )
        
        self.pressure = Var(
            initialize=101325.0,
            domain=NonNegativeReals,
            units=pyunits.Pa,
            doc="State pressure",
        )
        
        self.mole_frac_comp = Var(
            self.component_list,
            initialize=1.0 / len(self.component_list),
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Overall mole fraction by component",
        )
        
        # Closure equation for when state is NOT fully specified (via state block configuration argument defined_state = False)
        if not self.config.defined_state:
            self.sum_mole_frac = Constraint(
                expr=sum(self.mole_frac_comp[component] for component in self.component_list) == 1
            )
    
    
    # "Required" methods for manual property packages
    
    def define_state_vars(self):
        """Return the dictionary of primary state variables.

        Returns:
            dict: Mapping of state variable names to Pyomo variable objects.
        """
        return {
            "flow_mol": self.flow_mol,
            "temperature": self.temperature,
            "pressure": self.pressure,
            "mole_frac_comp": self.mole_frac_comp,
            }
    
    
    def get_material_flow_basis(self):
        """Return the material flow basis used by this property package.

        Returns:
            MaterialFlowBasis: MaterialFlowBasis.molar.
        """
        return MaterialFlowBasis.molar
    
    
    def get_material_flow_terms(self, phase, component):
        """Return material flow term for a given phase-component pair.

        Args:
            phase: Phase identifier.
            component: Component identifier.

        Returns:
            Expression-like term: Liquid flow contribution; zero for other phases.
        """
        if phase == "liquid":
            return self.flow_mol * self.mole_frac_comp[component]
        else:
            return 0 * pyunits.mol / pyunits.s
    
    
    def get_material_density_terms(self, phase, component):
        """Return component material density term for balance equations.

        Args:
            phase: Phase identifier.
            component: Component identifier.

        Returns:
            Expression-like term: Component molar density in liquid phase and
            zero for non-liquid phases.

        LaTeX form (liquid phase):
            c_i = x_i \\rho_mix / \\overline{MW}
            \\overline{MW} = \\sum_k x_k MW_k
        """
        if phase != "liquid":
            return 0 * pyunits.mol / pyunits.m**3
        mw_mix = sum(
            self.mole_frac_comp[c] * self.params.mw_comp[c]
            for c in self.component_list
        )
        return self.mole_frac_comp[component] * self.dens_mass / mw_mix
    
    
    def get_enthalpy_flow_terms(self, phase):
        """Return phase enthalpy flow term used in energy balances.

        Args:
            phase: Phase identifier.

        Returns:
            Expression-like term: Liquid enthalpy flow and zero for others.
        """
        if phase == "liquid":
            return self.flow_mol * self.enth_mol
        else:
            return 0 * pyunits.J / pyunits.s
    
    
    def get_energy_density_terms(self, phase):
        """Return phase energy density term for control-volume formulations.

        Args:
            phase: Phase identifier.

        Returns:
            Expression-like term: Approximate liquid internal-energy density and
            zero for non-liquid phases.

        LaTeX form:
            h = u + pv
            u_density = c_mol h - p
        """
        
        if phase == "liquid":
            # Total liquid molar density (mol/m^3)
            c_mol_liq = sum(
                self.get_material_density_terms(phase, component)
                for component in self.component_list
            )
            
            # Internal-energy density from h = u + p*v  -> u_density = c*h - p
            return c_mol_liq * self.enth_mol - self.pressure
        
        else:
            return 0 * pyunits.J / pyunits.m**3
    
    
    def default_material_balance_type(self):
        """Return the default material balance type.

        Returns:
            MaterialBalanceType: componentTotal.
        """
        return MaterialBalanceType.componentTotal
    
    
    def default_energy_balance_type(self):
        """Return the default energy balance type.

        Returns:
            EnergyBalanceType: enthalpyTotal.
        """
        return EnergyBalanceType.enthalpyTotal
    
    
    def define_port_members(self):
        """Return variables exposed through Ports for unit model connectivity.

        Returns:
            dict: Port member mapping using state variables.
        """
        return self.define_state_vars()
    
    
    def define_display_vars(self):
        """Return variables used for model display and reporting.

        Returns:
            dict: Display variable mapping using state variables.
        """
        return self.define_state_vars()
    
    
    # Build-on-demand methods for properties mentioned in metadata
    
    def _enth_mol(self):
        """Build mixture molar enthalpy expression for the liquid approximation.

        The expression combines standard liquid formation enthalpy at reference
        temperature with constant liquid heat-capacity correction.

        LaTeX form:
            h = \\sum_i x_i \\left(\\Delta h_{f,i}^{ref} + C_{p,i}^{liq}(T - T_ref)\\right)
        """
        self.enth_mol = Expression(
            expr=sum(
                self.mole_frac_comp[component]
                * (
                    self.params.dh_form_liq_comp[component]
                    + self.params.cp_mol_liq_comp[component]
                    * (self.temperature - self.params.temperature_ref)
                    
                )
                for component in self.component_list
            ),
            doc="Mixture molar enthalpy (liquid-phase approximation)",
        )
    
    
    def _dens_mass(self):
        """Build mixture mass-density expression from idealized volume mixing.

        Uses component molecular weights and fixed liquid densities to estimate
        mixture density in the liquid phase.

        LaTeX form:
            \\overline{MW} = \\sum_i x_i MW_i
            \\overline{V}_m = \\sum_i x_i MW_i / \\rho_i
            \\rho_mix = \\overline{MW} / \\overline{V}_m
        """
        mixture_mw = sum(
            self.mole_frac_comp[component] * self.params.mw_comp[component]
            for component in self.component_list
        )
        mixture_molar_volume = sum(
            self.mole_frac_comp[component]
            * self.params.mw_comp[component]
            / self.params.density_liq_comp[component]
            for component in self.component_list
        )
        self.dens_mass = Expression(
            expr=mixture_mw / mixture_molar_volume,
            doc="Mixture mass density (liquid_phase approximation)",
        )
    
    
    def _cp_mol(self):
        """Build mixture molar heat-capacity expression for liquid phase.

        Uses mole-fraction weighted averaging of fixed component liquid heat
        capacities.

        LaTeX form:
            C_{p,mix}^{liq} = \\sum_i x_i C_{p,i}^{liq}
        """
        self.cp_mol = Expression(
            expr=sum(
                self.mole_frac_comp[component]
                * self.params.cp_mol_liq_comp[component]
                for component in self.component_list
            ),
            doc="Mixture molar heat capacity (liquid-phase approximation)",
        )
    
    
    def _flow_mol_phase_comp(self):
        """Build phase-component molar-flow expression.

        Assumes all material flow is in the liquid phase and sets non-liquid
        phase-component flows to zero.
        """
        def flow_mol_phase_comp_rule(b, phase, component):
            if phase == "liquid":
                return b.flow_mol * b.mole_frac_comp[component]
            else:
                return 0 * pyunits.mol / pyunits.s
        
        self.flow_mol_phase_comp = Expression(
            self.phase_list,
            self.component_list,
            rule=flow_mol_phase_comp_rule,
            doc="Phase-component molar flow (liquid-only approximation)",
        )
    
    
    def _mole_frac_phase_comp(self):
        """Build phase-component mole-fraction expression.

        Assumes liquid-phase composition equals overall composition and sets
        non-liquid phase fractions to zero.
        """
        def mole_frac_phase_comp_rule(b, phase, component):
            if phase == "liquid":
                return b.mole_frac_comp[component]
            else:
                return 0.0
        
        self.mole_frac_phase_comp = Expression(
            self.phase_list,
            self.component_list,
            rule=mole_frac_phase_comp_rule,
            doc="Phase-component mole fraction (liquid-only approximation)",
        )
    
    
    def _phase_frac(self):
        """Build phase-fraction expression for a liquid-only approximation.

        Sets liquid phase fraction to one and all non-liquid phase fractions
        to zero.
        """
        def phase_frac_rule(b, phase):
            if phase == "liquid":
                return 1.0
            else:
                return 0.0
        self.phase_frac = Expression(
            self.phase_list,
            rule=phase_frac_rule,
            doc="Phase fraction (liquid-only approximation)",
        )
    
    
    def _act_coeff_liq_comp(self):
        """Build NRTL liquid-phase activity coefficients.

        LaTeX form:
            G_{ij} = \\exp(-\\alpha_{ij}\\tau_{ij})
            S_i = \\sum_k x_k G_{ki}
            N_i = \\sum_j x_j \\tau_{ji} G_{ji}
            Q_j = \\sum_k x_k G_{kj}
            P_j = \\sum_m x_m \\tau_{mj} G_{mj}
            W_{ij} = x_j G_{ij} / Q_j
            D_{ij} = \\tau_{ij} - P_j / Q_j
            \\ln(\\gamma_i) = N_i / S_i + \\sum_j W_{ij} D_{ij}
            \\gamma_i = \\exp(\\ln(\\gamma_i))
        """
        eps = 1e-12
        
        self.G_nrtl = Expression(
            self.component_list,
            self.component_list,
            rule=lambda b, i, j: exp(
                -b.params.alpha_nrtl[i, j] * b.params.tau_nrtl[i, j]
            ),
            doc="NRTL G_ij = exp(-alpha_ij * tau_ij)",
        )
        
        self.S_nrtl = Expression(
            self.component_list,
            rule=lambda b, i: sum(
                b.mole_frac_comp[k] * b.G_nrtl[k, i]
                for k in b.component_list
            ) + eps,
            doc="NRTL S_i = sum_k x_k G_ki",
        )
        
        self.N_nrtl = Expression(
            self.component_list,
            rule=lambda b, i: sum(
                b.mole_frac_comp[j]
                * b.params.tau_nrtl[j, i]
                * b.G_nrtl[j, i]
                for j in b.component_list
            ),
            doc="NRTL N_i = sum_j x_j tau_ji G_ji",
        )
        
        self.Q_nrtl = Expression(
            self.component_list,
            rule=lambda b, j: sum(
                b.mole_frac_comp[k] * b.G_nrtl[k, j]
                for k in b.component_list
            ) + eps,
            doc="NRTL Q_j = sum_k x_k G_kj",
        )
        
        self.P_nrtl = Expression(
            self.component_list,
            rule=lambda b, j: sum(
                b.mole_frac_comp[m]
                * b.params.tau_nrtl[m, j]
                * b.G_nrtl[m, j]
                for m in b.component_list
            ),
            doc="NRTL P_j = sum_m x_m tau_mj G_mj",
        )
        
        self.W_nrtl = Expression(
            self.component_list,
            self.component_list,
            rule=lambda b, i, j: (
                b.mole_frac_comp[j] * b.G_nrtl[i, j] / b.Q_nrtl[j]
            ),
            doc="NRTL W_ij = x_j G_ij / Q_j",
        )
        
        self.D_nrtl = Expression(
            self.component_list,
            self.component_list,
            rule=lambda b, i, j: (
                b.params.tau_nrtl[i, j] - b.P_nrtl[j] / b.Q_nrtl[j]
            ),
            doc="NRTL D_ij = tau_ij - P_j/Q_j",
        )
        
        self.log_gamma_liq_comp = Expression(
            self.component_list,
            rule=lambda b, i: (
                b.N_nrtl[i] / b.S_nrtl[i]
                + sum(b.W_nrtl[i, j] * b.D_nrtl[i, j] for j in b.component_list)
            ),
            doc="NRTL ln(gamma_i)",
        )
        
        self.act_coeff_liq_comp = Expression(
            self.component_list,
            rule=lambda b, i: exp(b.log_gamma_liq_comp[i]),
            doc="Liquid activity coefficient gamma_i from NRTL",
        )
