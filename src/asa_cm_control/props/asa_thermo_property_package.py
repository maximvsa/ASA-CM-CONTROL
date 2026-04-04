# To-Do:
# - 

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

@declare_process_block_class("ASAThermoParameterBlock")
class ASAThermoParameterData(PhysicalParameterBlock):
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
        self._state_block_class = ASAThermoStateBlock  # pyright: ignore[reportUndefinedVariable]
        
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
                "aspirin": 0.18016,
                "acetic_acid": 0.06005,
                "water": 0.01801528,
                "H2SO4": 0.098079,
                "H_plus": 0.001008,
                "HSO4_minus": 0.097071,
                "SO4_2minus": 0.096063,
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
                "H2SO4": 137.6,
                "H_plus": 75.37,
                "HSO4_minus": 135.0,
                "SO4_2minus": 20.0,
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
                "aspirin": 149.0,
                "acetic_acid": 63.44,
                "water": 33.58,
                "H2SO4": 83.71,
                "H_plus": 0,
                "HSO4_minus": 0,
                "SO4_2minus": 0,
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
                "aspirin": 217.8,
                "acetic_acid": 123.1,
                "water": 37.77,
                "H2SO4": 137.6,
                "H_plus": 0,
                "HSO4_minus": 0,
                "SO4_2minus": 0,
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
                "aspirin": 1330.0,
                "acetic_acid": 1053.0,
                "water": 1000.0,
                "H2SO4": 1840.0,
                "H_plus": 0,
                "HSO4_minus": 0,
                "SO4_2minus": 0,
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
                "aspirin": 1335.0,
                "acetic_acid": 1266.0,
                "water": 916.7,
                "H2SO4": 1940.0,
                "H_plus": 0,
                "HSO4_minus": 0,
                "SO4_2minus": 0,
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
                "aspirin": -758200,
                "acetic_acid": -484000,
                "water": -286000,
                "H2SO4": -814000,
                "H_plus": 0,
                "HSO4_minus": -907500,
                "SO4_2minus": -882000,
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
            doc="",
        )
        
        self.charge_comp.fix()
        
        neutral_components = [
            "salicylic_acid",
            "acetic_anhydride",
            "aspirin",
            "acetic_acid",
            "water",
            "H2SO4",
        ]
        
        self.neutral_component_set = Set(initialize=neutral_components)
        
        self.molecular_component_set = Set(
            initialize=[
                "salicylic_acid",
                "acetic_anhydride",
                "aspirin",
                "acetic_acid",
                "water",
            ]
        )
        self.ionic_component_set = Set(
            initialize=[
                "H_plus",
                "HSO4_minus",
                "SO4_2minus",
            ]
        )
        self.neutral_pair_set = Set(
            dimen=2,
            initialize=[(i, j) for i in neutral_components for j in neutral_components],
        )
        
        self.tau_nrtl = Var(
            self.neutral_pair_set,
            initialize=0.0,
            domain=Reals,
            units=pyunits.dimensionless
        )
        
        self.tau_nrtl.fix()
        
        self.alpha_nrtl = Var(
            self.neutral_pair_set,
            initialize=0.3,
            domain=NonNegativeReals,
            units=pyunits.dimensionless
        )
        
        self.alpha_nrtl.fix()
        
        for c in self.neutral_component_set:
            self.tau_nrtl[c, c].fix(0.0)
            self.alpha_nrtl[c, c].fix(0.0)
        

        self.A_Davies = Var(
            initialize=0.509,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Davies correlation constant at ~25C",
        )

        self.A_Davies.fix()
        
        # -----------------------------
        # tau_ij seeds (dimensionless)
        # -----------------------------
        
        self.tau_nrtl["salicylic_acid", "acetic_anhydride"].fix(0.35)
        self.tau_nrtl["acetic_anhydride", "salicylic_acid"].fix(0.15)
        
        self.tau_nrtl["salicylic_acid", "H2SO4"].fix(0.10)
        self.tau_nrtl["H2SO4", "salicylic_acid"].fix(0.80)
        
        self.tau_nrtl["salicylic_acid", "aspirin"].fix(0.20)
        self.tau_nrtl["aspirin", "salicylic_acid"].fix(0.15)
        
        self.tau_nrtl["salicylic_acid", "acetic_acid"].fix(0.55)
        self.tau_nrtl["acetic_acid", "salicylic_acid"].fix(0.30)
        
        self.tau_nrtl["salicylic_acid", "water"].fix(1.80)
        self.tau_nrtl["water", "salicylic_acid"].fix(0.60)
        
        self.tau_nrtl["acetic_anhydride", "H2SO4"].fix(0.20)
        self.tau_nrtl["H2SO4", "acetic_anhydride"].fix(0.90)
        
        self.tau_nrtl["acetic_anhydride", "aspirin"].fix(0.45)
        self.tau_nrtl["aspirin", "acetic_anhydride"].fix(0.20)
        
        self.tau_nrtl["acetic_anhydride", "acetic_acid"].fix(0.25)
        self.tau_nrtl["acetic_acid", "acetic_anhydride"].fix(0.10)
        
        self.tau_nrtl["acetic_anhydride", "water"].fix(2.20)
        self.tau_nrtl["water", "acetic_anhydride"].fix(0.40)
        
        self.tau_nrtl["H2SO4", "aspirin"].fix(0.20)
        self.tau_nrtl["aspirin", "H2SO4"].fix(0.90)
        
        self.tau_nrtl["H2SO4", "acetic_acid"].fix(-0.20)
        self.tau_nrtl["acetic_acid", "H2SO4"].fix(0.80)
        
        self.tau_nrtl["H2SO4", "water"].fix(-0.80)
        self.tau_nrtl["water", "H2SO4"].fix(2.50)
        
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
        
        self.alpha_nrtl["H2SO4", "water"].fix(0.20)
        self.alpha_nrtl["water", "H2SO4"].fix(0.20)
        
        self.alpha_nrtl["H2SO4", "acetic_acid"].fix(0.20)
        self.alpha_nrtl["acetic_acid", "H2SO4"].fix(0.20)
        
        self.alpha_nrtl["H2SO4", "acetic_anhydride"].fix(0.20)
        self.alpha_nrtl["acetic_anhydride", "H2SO4"].fix(0.20)
        
        self.alpha_nrtl["H2SO4", "salicylic_acid"].fix(0.20)
        self.alpha_nrtl["salicylic_acid", "H2SO4"].fix(0.20)
        
        self.alpha_nrtl["H2SO4", "aspirin"].fix(0.20)
        self.alpha_nrtl["aspirin", "H2SO4"].fix(0.20)
        
        

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
                
                'act_coeff_true_comp': {'method': '_act_coeff_true_comp'},
                'activity_true_comp': {'method': '_activity_true_comp'},
                
                'ionic_strength': {'method': '_ionic_strength'},
                'log10_gamma_davies_ion_comp': {'method': '_log10_gamma_davies_ion_comp'},
                'gamma_davies_ion_comp': {'method': '_gamma_davies_ion_comp'},
                
                'a_H_plus': {'method': '_a_H_plus'},
                'x_H_plus': {'method': '_x_H_plus'},
                'gamma_H_plus': {'method': '_gamma_H_plus'},
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

class ASAThermoStateBlockMethods(StateBlock):
    
    def initialize(
        self,
        state_args=None,
        state_vars_fixed=False,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
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



@declare_process_block_class("ASAThermoStateBlock", block_class=ASAThermoStateBlockMethods)
class ASAThermoStateBlockData(StateBlockData):
    
    def build(self):
        super().build()
        
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
            initialize={
                c: (1.0 / len(self.params.neutral_component_set))
                if c in self.params.neutral_component_set
                else 0.0
                for c in self.component_list
            },
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="",
        )
        
        if not self.config.defined_state:
            self.sum_mole_frac = Constraint(
                expr=sum(self.mole_frac_comp[component] for component in self.component_list) == 1
            )
        
        self.mole_frac_comp_true = Expression(
            self.component_list,
            rule=lambda block, component: block.mole_frac_comp[component],
            doc="",
        )
    
    def define_state_vars(self):
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
            c_i = x_i \rho_{mix} / \overline{MW} \\\\
            \overline{MW} = \sum_k x_k MW_k
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
            h = u + pv \\\\
            u_{density} = c_{mol} h - p
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
            h = \sum_i x_i \left(\Delta h_{f,i}^{ref} + C_{p,i}^{liq}(T - T_{ref})\right)
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
            \overline{MW} = \sum_i x_i MW_i \\\\
            \overline{V}_m = \sum_i x_i MW_i / \rho_i \\\\
            \rho_{mix} = \overline{MW} / \overline{V}_m
        """
        mixture_mw = sum(
            self.mole_frac_comp[component] * self.params.mw_comp[component]
            for component in self.params.neutral_component_set
        )
        mixture_molar_volume = sum(
            self.mole_frac_comp[component]
            * self.params.mw_comp[component]
            / self.params.density_liq_comp[component]
            for component in self.params.neutral_component_set
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
            C_{p,mix}^{liq} = \sum_i x_i C_{p,i}^{liq}
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
            G_{ij} = \exp(-\alpha_{ij}\tau_{ij}) \\\\
            S_i = \sum_k x_k G_{ki} \\\\
            N_i = \sum_j x_j \tau_{ji} G_{ji} \\\\
            Q_j = \sum_k x_k G_{kj} \\\\
            P_j = \sum_m x_m \tau_{mj} G_{mj} \\\\
            W_{ij} = x_j G_{ij} / Q_j \\\\
            D_{ij} = \tau_{ij} - P_j / Q_j \\\\
            \ln(\gamma_i) = N_i / S_i + \sum_j W_{ij} D_{ij} \\\\
            \gamma_i = \exp(\ln(\gamma_i))
        """
        eps = 1e-12
        
        self.G_nrtl = Expression(
            self.params.neutral_pair_set,
            rule=lambda b, i, j: exp(
                -b.params.alpha_nrtl[i, j] * b.params.tau_nrtl[i, j]
            ),
            doc="NRTL G_ij = exp(-alpha_ij * tau_ij)",
        )
        
        self.S_nrtl = Expression(
            self.params.neutral_component_set,
            rule=lambda b, i: sum(
                b.mole_frac_comp_true[k] * b.G_nrtl[k, i]
                for k in b.params.neutral_component_set
            ) + eps,
            doc="NRTL S_i = sum_k x_k G_ki",
        )
        
        self.N_nrtl = Expression(
            self.params.neutral_component_set,
            rule=lambda b, i: sum(
                b.mole_frac_comp_true[j]
                * b.params.tau_nrtl[j, i]
                * b.G_nrtl[j, i]
                for j in b.params.neutral_component_set
            ),
            doc="NRTL N_i = sum_j x_j tau_ji G_ji",
        )
        
        self.Q_nrtl = Expression(
            self.params.neutral_component_set,
            rule=lambda b, j: sum(
                b.mole_frac_comp_true[k] * b.G_nrtl[k, j]
                for k in b.params.neutral_component_set
            ) + eps,
            doc="NRTL Q_j = sum_k x_k G_kj",
        )
        
        self.P_nrtl = Expression(
            self.params.neutral_component_set,
            rule=lambda b, j: sum(
                b.mole_frac_comp_true[m]
                * b.params.tau_nrtl[m, j]
                * b.G_nrtl[m, j]
                for m in b.params.neutral_component_set
            ),
            doc="NRTL P_j = sum_m x_m tau_mj G_mj",
        )
        
        self.W_nrtl = Expression(
            self.params.neutral_pair_set,
            rule=lambda b, i, j: (
                b.mole_frac_comp_true[j] * b.G_nrtl[i, j] / b.Q_nrtl[j]
            ),
            doc="NRTL W_ij = x_j G_ij / Q_j",
        )
        
        self.D_nrtl = Expression(
            self.params.neutral_pair_set,
            rule=lambda b, i, j: (
                b.params.tau_nrtl[i, j] - b.P_nrtl[j] / b.Q_nrtl[j]
            ),
            doc="NRTL D_ij = tau_ij - P_j/Q_j",
        )
        
        self.log_gamma_liq_comp = Expression(
            self.component_list,
            rule=lambda b, i: (
                b.N_nrtl[i] / b.S_nrtl[i]
                + sum(
                    b.W_nrtl[i, j] * b.D_nrtl[i, j]
                    for j in b.params.neutral_component_set
                )
            )
            if i in b.params.neutral_component_set
            else 0.0,
            doc="NRTL ln(gamma_i)",
        )
        
        self.act_coeff_liq_comp = Expression(
            self.component_list,
            rule=lambda b, i: exp(b.log_gamma_liq_comp[i])
            if i in b.params.neutral_component_set
            else 1.0,
            doc="Liquid activity coefficient gamma_i from NRTL",
        )
    
    def _ionic_strength(self):
        self.ionic_strength = Expression(
            expr=0.5 * sum(
                self.params.charge_comp[ion] ** 2 * self.mole_frac_comp_true[ion]
                for ion in self.params.ionic_component_set
            )
        )
    
    def _log10_gamma_davies_ion_comp(self):
        epsilon = 1e-8
        
        def rule(block, ion):
            sqrt_i = (block.ionic_strength + epsilon) ** 0.5
            davies_term = sqrt_i / (1 + sqrt_i) - 0.3 * block.ionic_strength
            return -block.params.A_Davies * (block.params.charge_comp[ion] ** 2) * davies_term
        
        self.log10_gamma_davies_ion_comp = Expression(
            self.params.ionic_component_set,
            rule=rule,
        )
    
    def _gamma_davies_ion_comp(self):
        self.gamma_davies_ion_comp = Expression(
            self.params.ionic_component_set,
            rule=lambda block, ion: 10 ** block.log10_gamma_davies_ion_comp[ion],
        )
    
    def _x_H_plus(self):
        self.x_H_plus = Expression(expr=self.mole_frac_comp_true["H_plus"])
    
    def _gamma_H_plus(self):
        self.gamma_H_plus = Expression(expr=self.gamma_davies_ion_comp["H_plus"])
    
    def _act_coeff_true_comp(self):
        self.act_coeff_true_comp = Expression(
            self.component_list,
            rule=lambda block, component: block.act_coeff_liq_comp[component]
            if component in block.params.neutral_component_set
            else (
                block.gamma_davies_ion_comp[component]
                if component in block.params.ionic_component_set
                else 1.0
            ),
        )
    
    def _activity_true_comp(self):
        self.activity_true_comp = Expression(
            self.component_list,
            rule=lambda block, component: block.act_coeff_true_comp[component] * block.mole_frac_comp_true[component],
        )
    
    def _a_H_plus(self):
        self.a_H_plus = Expression(expr=self.activity_true_comp["H_plus"])
