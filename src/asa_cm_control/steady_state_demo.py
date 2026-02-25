from pathlib import Path
import sys

SRC_DIR = Path(__file__).resolve().parents[1]
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pyomo.environ import (
    Constraint,
    ConcreteModel,
    Expression,
    Param,
    Var,
    check_optimal_termination,
    units as pyunits,
    value,
)
from pyomo.network import Arc, Port, SequentialDecomposition

from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)
from idaes.models.unit_models import (
    CSTR,
    Feed,
    Heater,
    Mixer,
    Product,
    Pump,
    Separator,
    StateJunction,
    StoichiometricReactor,
    Valve,
    ValveFunctionType,
)
from idaes.models.unit_models.separator import SplittingType
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from pyomo.environ import TransformationFactory

from asa_cm_control.props.reaction_config import configuration as reaction_config
from asa_cm_control.props.thermo_config import configuration as thermo_config


def _tank_unit(property_package, reaction_package, has_heat_transfer=False):
    return CSTR(
        property_package=property_package,
        reaction_package=reaction_package,
        has_equilibrium_reactions=False,
        has_phase_equilibrium=False,
        has_heat_of_reaction=False,
        has_heat_transfer=has_heat_transfer,
        has_pressure_change=False,
    )


def build_model() -> ConcreteModel:
    m = ConcreteModel(name="asa_cm_control_steady_state_demo")
    m.fs = FlowsheetBlock(dynamic=False)
    t0 = m.fs.time.first()

    m.fs.thermo_params = GenericParameterBlock(**thermo_config)
    m.fs.reaction_params = GenericReactionParameterBlock(
        property_package=m.fs.thermo_params,
        **reaction_config,
    )

    m.fs.sa_makeup_feed_100 = Feed(property_package=m.fs.thermo_params)
    m.fs.aa_makeup_feed_105 = Feed(property_package=m.fs.thermo_params)
    m.fs.h2so4_makeup_feed_145 = Feed(property_package=m.fs.thermo_params)
    m.fs.wash_makeup_feed_404 = Feed(property_package=m.fs.thermo_params)

    m.fs.sa_day_hopper_110 = _tank_unit(
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.reaction_params,
    )
    m.fs.sa_loss_in_weight_feeder_115 = StateJunction(property_package=m.fs.thermo_params)
    m.fs.sa_loss_in_weight_feeder_115.setpoint_flow_mol = Var(
        initialize=1.0,
        units=pyunits.mol / pyunits.s,
    )
    m.fs.sa_loss_in_weight_feeder_115.bias_flow_mol = Var(
        initialize=0.0,
        units=pyunits.mol / pyunits.s,
    )
    m.fs.sa_loss_in_weight_feeder_115.bias_flow_mol.fix(0.0)
    m.fs.sa_loss_in_weight_feeder_115.metering_eq = Constraint(
        expr=m.fs.sa_loss_in_weight_feeder_115.outlet.flow_mol[t0]
        == m.fs.sa_loss_in_weight_feeder_115.setpoint_flow_mol
        + m.fs.sa_loss_in_weight_feeder_115.bias_flow_mol
    )

    m.fs.aa_day_tank_feed_mixer_120 = Mixer(
        property_package=m.fs.thermo_params,
        inlet_list=["aa_makeup_in", "aa_recycle_in"],
    )
    m.fs.aa_day_tank_120 = _tank_unit(
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.reaction_params,
    )
    m.fs.aa_feed_pump_121 = Pump(
        property_package=m.fs.thermo_params,
        thermodynamic_assumption=ThermodynamicAssumption.isothermal,
    )
    m.fs.aa_inline_filter_122 = Heater(
        property_package=m.fs.thermo_params,
        has_pressure_change=False,
    )
    m.fs.aa_coriolis_ft_123 = Heater(
        property_package=m.fs.thermo_params,
        has_pressure_change=False,
    )
    m.fs.aa_coriolis_ft_123.measured_flow_mol = Var(
        initialize=1.0,
        units=pyunits.mol / pyunits.s,
    )
    m.fs.aa_coriolis_ft_123.bias_flow_mol = Var(
        initialize=0.0,
        units=pyunits.mol / pyunits.s,
    )
    m.fs.aa_coriolis_ft_123.measurement_eq = Constraint(
        expr=m.fs.aa_coriolis_ft_123.measured_flow_mol
        == m.fs.aa_coriolis_ft_123.outlet.flow_mol[t0]
        + m.fs.aa_coriolis_ft_123.bias_flow_mol
    )

    m.fs.aa_flow_control_valve_124 = Valve(
        property_package=m.fs.thermo_params,
        valve_function_callback=ValveFunctionType.linear,
    )
    m.fs.sa_aa_wetting_mixer_130 = Mixer(
        property_package=m.fs.thermo_params,
        inlet_list=["aa_in", "sa_in"],
    )
    m.fs.slurry_make_tank_140 = _tank_unit(
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.reaction_params,
    )
    m.fs.reactor_feed_pump_141 = Pump(
        property_package=m.fs.thermo_params,
        thermodynamic_assumption=ThermodynamicAssumption.isothermal,
    )
    m.fs.reactor_feed_ft_142 = Heater(
        property_package=m.fs.thermo_params,
        has_pressure_change=False,
    )
    m.fs.reactor_feed_ft_142.measured_flow_mol = Var(
        initialize=1.0,
        units=pyunits.mol / pyunits.s,
    )
    m.fs.reactor_feed_ft_142.bias_flow_mol = Var(
        initialize=0.0,
        units=pyunits.mol / pyunits.s,
    )
    m.fs.reactor_feed_ft_142.measurement_eq = Constraint(
        expr=m.fs.reactor_feed_ft_142.measured_flow_mol
        == m.fs.reactor_feed_ft_142.outlet.flow_mol[t0]
        + m.fs.reactor_feed_ft_142.bias_flow_mol
    )

    m.fs.h2so4_day_tank_150 = _tank_unit(
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.reaction_params,
    )
    m.fs.h2so4_metering_pump_151 = Pump(
        property_package=m.fs.thermo_params,
        thermodynamic_assumption=ThermodynamicAssumption.isothermal,
    )
    m.fs.h2so4_flow_ft_152 = Heater(
        property_package=m.fs.thermo_params,
        has_pressure_change=False,
    )
    m.fs.h2so4_flow_ft_152.measured_flow_mol = Var(
        initialize=1.0,
        units=pyunits.mol / pyunits.s,
    )
    m.fs.h2so4_flow_ft_152.bias_flow_mol = Var(
        initialize=0.0,
        units=pyunits.mol / pyunits.s,
    )
    m.fs.h2so4_flow_ft_152.measurement_eq = Constraint(
        expr=m.fs.h2so4_flow_ft_152.measured_flow_mol
        == m.fs.h2so4_flow_ft_152.outlet.flow_mol[t0]
        + m.fs.h2so4_flow_ft_152.bias_flow_mol
    )

    m.fs.reactor_inlet_static_mixer_210 = Mixer(
        property_package=m.fs.thermo_params,
        inlet_list=["slurry_in", "acid_in"],
    )
    m.fs.asa_acetylation_reactor_220 = StoichiometricReactor(
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.reaction_params,
        has_heat_of_reaction=True,
        has_heat_transfer=False,
        has_pressure_change=False,
    )

    m.fs.asa_seed_takeoff_splitter_305 = Separator(
        property_package=m.fs.thermo_params,
        outlet_list=["product_out", "seed_recycle_out"],
    )
    m.fs.asa_seed_slurry_skid_310 = Mixer(
        property_package=m.fs.thermo_params,
        inlet_list=["seed_makeup_in", "solvent_in"],
    )
    m.fs.seed_distribution_manifold_315 = Separator(
        property_package=m.fs.thermo_params,
        outlet_list=["seed1_out", "seed2_out"],
    )
    m.fs.cooled_crystallizer_1_320 = _tank_unit(
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.reaction_params,
        has_heat_transfer=True,
    )
    m.fs.cooled_crystallizer_2_330 = _tank_unit(
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.reaction_params,
        has_heat_transfer=True,
    )

    m.fs.cake_wash_solvent_feed_mixer_405 = Mixer(
        property_package=m.fs.thermo_params,
        inlet_list=["wash_makeup_in", "wash_recycle_in"],
    )
    m.fs.cake_wash_solvent_tank_405 = _tank_unit(
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.reaction_params,
    )

    m.fs.pusher_centrifuge_feed_mixer_410 = Mixer(
        property_package=m.fs.thermo_params,
        inlet_list=["slurry_in", "wash_in"],
    )
    m.fs.pusher_centrifuge_splitter_410 = Separator(
        property_package=m.fs.thermo_params,
        outlet_list=["wet_cake_out", "centrate_out"],
    )

    m.fs.mother_liquor_receiver_430 = _tank_unit(
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.reaction_params,
    )
    m.fs.solvent_recovery_column_feed_mixer_440 = Mixer(
        property_package=m.fs.thermo_params,
        inlet_list=["feed_in", "reflux_in"],
    )
    m.fs.solvent_recovery_column_splitter_440 = Separator(
        property_package=m.fs.thermo_params,
        outlet_list=["overhead_vapor_out", "bottoms_out"],
    )
    m.fs.overhead_condenser_441 = Heater(
        property_package=m.fs.thermo_params,
        has_pressure_change=False,
    )
    m.fs.reflux_drum_receiver_442 = Separator(
        property_package=m.fs.thermo_params,
        outlet_list=["reflux_out", "distillate_out"],
    )
    m.fs.recovered_volatiles_tank_450 = _tank_unit(
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.reaction_params,
    )
    m.fs.recycle_splitter_455 = Separator(
        property_package=m.fs.thermo_params,
        split_basis=SplittingType.componentFlow,
        outlet_list=["aa_recycle_out", "wash_recycle_out", "seed_solvent_out"],
    )
    m.fs.column_bottoms_purge_460 = Product(property_package=m.fs.thermo_params)

    m.fs.fluidized_bed_dryer_510 = Heater(
        property_package=m.fs.thermo_params,
        has_pressure_change=False,
    )
    m.fs.asa_product_bin_520 = Product(property_package=m.fs.thermo_params)

    m.fs.s100_110_sa_makeup = Arc(
        source=m.fs.sa_makeup_feed_100.outlet,
        destination=m.fs.sa_day_hopper_110.inlet,
    )
    m.fs.s110_115_sa = Arc(
        source=m.fs.sa_day_hopper_110.outlet,
        destination=m.fs.sa_loss_in_weight_feeder_115.inlet,
    )
    m.fs.s115_130_sa = Arc(
        source=m.fs.sa_loss_in_weight_feeder_115.outlet,
        destination=m.fs.sa_aa_wetting_mixer_130.sa_in,
    )

    m.fs.s105_120_aa_makeup = Arc(
        source=m.fs.aa_makeup_feed_105.outlet,
        destination=m.fs.aa_day_tank_feed_mixer_120.aa_makeup_in,
    )
    m.fs.s120mix_120tank_aa = Arc(
        source=m.fs.aa_day_tank_feed_mixer_120.outlet,
        destination=m.fs.aa_day_tank_120.inlet,
    )
    m.fs.s120_121_aa = Arc(
        source=m.fs.aa_day_tank_120.outlet,
        destination=m.fs.aa_feed_pump_121.inlet,
    )
    m.fs.s121_122_aa = Arc(
        source=m.fs.aa_feed_pump_121.outlet,
        destination=m.fs.aa_inline_filter_122.inlet,
    )
    m.fs.s122_123_aa = Arc(
        source=m.fs.aa_inline_filter_122.outlet,
        destination=m.fs.aa_coriolis_ft_123.inlet,
    )
    m.fs.s123_124_aa = Arc(
        source=m.fs.aa_coriolis_ft_123.outlet,
        destination=m.fs.aa_flow_control_valve_124.inlet,
    )
    m.fs.s124_130_aa = Arc(
        source=m.fs.aa_flow_control_valve_124.outlet,
        destination=m.fs.sa_aa_wetting_mixer_130.aa_in,
    )

    m.fs.s130_140_slurry = Arc(
        source=m.fs.sa_aa_wetting_mixer_130.outlet,
        destination=m.fs.slurry_make_tank_140.inlet,
    )
    m.fs.s140_141_slurry = Arc(
        source=m.fs.slurry_make_tank_140.outlet,
        destination=m.fs.reactor_feed_pump_141.inlet,
    )
    m.fs.s141_142_slurry = Arc(
        source=m.fs.reactor_feed_pump_141.outlet,
        destination=m.fs.reactor_feed_ft_142.inlet,
    )
    m.fs.s142_210_slurry = Arc(
        source=m.fs.reactor_feed_ft_142.outlet,
        destination=m.fs.reactor_inlet_static_mixer_210.slurry_in,
    )

    m.fs.s145_150_h2so4_makeup = Arc(
        source=m.fs.h2so4_makeup_feed_145.outlet,
        destination=m.fs.h2so4_day_tank_150.inlet,
    )
    m.fs.s150_151_h2so4 = Arc(
        source=m.fs.h2so4_day_tank_150.outlet,
        destination=m.fs.h2so4_metering_pump_151.inlet,
    )
    m.fs.s151_152_h2so4 = Arc(
        source=m.fs.h2so4_metering_pump_151.outlet,
        destination=m.fs.h2so4_flow_ft_152.inlet,
    )
    m.fs.s152_210_h2so4 = Arc(
        source=m.fs.h2so4_flow_ft_152.outlet,
        destination=m.fs.reactor_inlet_static_mixer_210.acid_in,
    )

    m.fs.s210_220_feed = Arc(
        source=m.fs.reactor_inlet_static_mixer_210.outlet,
        destination=m.fs.asa_acetylation_reactor_220.inlet,
    )
    m.fs.s220_320_rxn = Arc(
        source=m.fs.asa_acetylation_reactor_220.outlet,
        destination=m.fs.cooled_crystallizer_1_320.inlet,
    )

    m.fs.s320_330_slurry1 = Arc(
        source=m.fs.cooled_crystallizer_1_320.outlet,
        destination=m.fs.cooled_crystallizer_2_330.inlet,
    )
    m.fs.s330_410_slurry2 = Arc(
        source=m.fs.cooled_crystallizer_2_330.outlet,
        destination=m.fs.pusher_centrifuge_feed_mixer_410.slurry_in,
    )

    m.fs.s404_405_wash_makeup = Arc(
        source=m.fs.wash_makeup_feed_404.outlet,
        destination=m.fs.cake_wash_solvent_feed_mixer_405.wash_makeup_in,
    )
    m.fs.s405mix_405tank = Arc(
        source=m.fs.cake_wash_solvent_feed_mixer_405.outlet,
        destination=m.fs.cake_wash_solvent_tank_405.inlet,
    )
    m.fs.s405_410_wash = Arc(
        source=m.fs.cake_wash_solvent_tank_405.outlet,
        destination=m.fs.pusher_centrifuge_feed_mixer_410.wash_in,
    )
    m.fs.s410mix_410split = Arc(
        source=m.fs.pusher_centrifuge_feed_mixer_410.outlet,
        destination=m.fs.pusher_centrifuge_splitter_410.inlet,
    )
    m.fs.s410_510_wet_cake = Arc(
        source=m.fs.pusher_centrifuge_splitter_410.wet_cake_out,
        destination=m.fs.fluidized_bed_dryer_510.inlet,
    )
    m.fs.s410_430_centrate = Arc(
        source=m.fs.pusher_centrifuge_splitter_410.centrate_out,
        destination=m.fs.mother_liquor_receiver_430.inlet,
    )

    m.fs.s430_440_feed = Arc(
        source=m.fs.mother_liquor_receiver_430.outlet,
        destination=m.fs.solvent_recovery_column_feed_mixer_440.feed_in,
    )
    m.fs.s442_440_reflux = Arc(
        source=m.fs.reflux_drum_receiver_442.reflux_out,
        destination=m.fs.solvent_recovery_column_feed_mixer_440.reflux_in,
    )
    m.fs.s440mix_440split = Arc(
        source=m.fs.solvent_recovery_column_feed_mixer_440.outlet,
        destination=m.fs.solvent_recovery_column_splitter_440.inlet,
    )
    m.fs.s440_441_overhead = Arc(
        source=m.fs.solvent_recovery_column_splitter_440.overhead_vapor_out,
        destination=m.fs.overhead_condenser_441.inlet,
    )
    m.fs.s441_442_condensate = Arc(
        source=m.fs.overhead_condenser_441.outlet,
        destination=m.fs.reflux_drum_receiver_442.inlet,
    )
    m.fs.s442_450_distillate = Arc(
        source=m.fs.reflux_drum_receiver_442.distillate_out,
        destination=m.fs.recovered_volatiles_tank_450.inlet,
    )
    m.fs.s450_455_recycle = Arc(
        source=m.fs.recovered_volatiles_tank_450.outlet,
        destination=m.fs.recycle_splitter_455.inlet,
    )
    m.fs.s455_120_aa_recycle = Arc(
        source=m.fs.recycle_splitter_455.aa_recycle_out,
        destination=m.fs.aa_day_tank_feed_mixer_120.aa_recycle_in,
    )
    m.fs.s455_405_wash_recycle = Arc(
        source=m.fs.recycle_splitter_455.wash_recycle_out,
        destination=m.fs.cake_wash_solvent_feed_mixer_405.wash_recycle_in,
    )
    m.fs.s455_310_seed_solvent = Arc(
        source=m.fs.recycle_splitter_455.seed_solvent_out,
        destination=m.fs.asa_seed_slurry_skid_310.solvent_in,
    )
    m.fs.s440_460_bottoms = Arc(
        source=m.fs.solvent_recovery_column_splitter_440.bottoms_out,
        destination=m.fs.column_bottoms_purge_460.inlet,
    )

    m.fs.s510_305_dry_powder = Arc(
        source=m.fs.fluidized_bed_dryer_510.outlet,
        destination=m.fs.asa_seed_takeoff_splitter_305.inlet,
    )
    m.fs.s305_520_product = Arc(
        source=m.fs.asa_seed_takeoff_splitter_305.product_out,
        destination=m.fs.asa_product_bin_520.inlet,
    )
    m.fs.s305_310_seed_solids = Arc(
        source=m.fs.asa_seed_takeoff_splitter_305.seed_recycle_out,
        destination=m.fs.asa_seed_slurry_skid_310.seed_makeup_in,
    )
    m.fs.s310_315_seed_slurry = Arc(
        source=m.fs.asa_seed_slurry_skid_310.outlet,
        destination=m.fs.seed_distribution_manifold_315.inlet,
    )
    m.fs.s315_320_seed1 = Arc(
        source=m.fs.seed_distribution_manifold_315.seed1_out,
        destination=m.fs.cooled_crystallizer_1_320.inlet,
    )
    m.fs.s315_330_seed2 = Arc(
        source=m.fs.seed_distribution_manifold_315.seed2_out,
        destination=m.fs.cooled_crystallizer_2_330.inlet,
    )

    eps = 1e-8
    major = 1.0 - 4.0 * eps

    m.fs.sa_makeup_feed_100.outlet.temperature[t0].fix(298.15)
    m.fs.sa_makeup_feed_100.outlet.pressure[t0].fix(101325.0)
    m.fs.sa_makeup_feed_100.outlet.mole_frac_comp[t0, "AS"].fix(major)
    m.fs.sa_makeup_feed_100.outlet.mole_frac_comp[t0, "AA"].fix(eps)
    m.fs.sa_makeup_feed_100.outlet.mole_frac_comp[t0, "ASA"].fix(eps)
    m.fs.sa_makeup_feed_100.outlet.mole_frac_comp[t0, "AcOH"].fix(eps)
    m.fs.sa_makeup_feed_100.outlet.mole_frac_comp[t0, "H2SO4"].fix(eps)

    m.fs.sa_makeup_feed_100.outlet.flow_mol[t0].fix(1.0)

    m.fs.aa_makeup_feed_105.outlet.temperature[t0].fix(298.15)
    m.fs.aa_makeup_feed_105.outlet.pressure[t0].fix(101325.0)
    m.fs.aa_makeup_feed_105.outlet.mole_frac_comp[t0, "AS"].fix(eps)
    m.fs.aa_makeup_feed_105.outlet.mole_frac_comp[t0, "AA"].fix(major)
    m.fs.aa_makeup_feed_105.outlet.mole_frac_comp[t0, "ASA"].fix(eps)
    m.fs.aa_makeup_feed_105.outlet.mole_frac_comp[t0, "AcOH"].fix(eps)
    m.fs.aa_makeup_feed_105.outlet.mole_frac_comp[t0, "H2SO4"].fix(eps)

    m.fs.aa_makeup_feed_105.outlet.flow_mol[t0].fix(2.0)

    m.fs.h2so4_makeup_feed_145.outlet.flow_mol[t0].fix(0.01)
    m.fs.h2so4_makeup_feed_145.outlet.temperature[t0].fix(298.15)
    m.fs.h2so4_makeup_feed_145.outlet.pressure[t0].fix(101325.0)
    m.fs.h2so4_makeup_feed_145.outlet.mole_frac_comp[t0, "AS"].fix(eps)
    m.fs.h2so4_makeup_feed_145.outlet.mole_frac_comp[t0, "AA"].fix(eps)
    m.fs.h2so4_makeup_feed_145.outlet.mole_frac_comp[t0, "ASA"].fix(eps)
    m.fs.h2so4_makeup_feed_145.outlet.mole_frac_comp[t0, "AcOH"].fix(eps)
    m.fs.h2so4_makeup_feed_145.outlet.mole_frac_comp[t0, "H2SO4"].fix(major)

    m.fs.wash_makeup_feed_404.outlet.flow_mol[t0].fix(0.10)
    m.fs.wash_makeup_feed_404.outlet.temperature[t0].fix(298.15)
    m.fs.wash_makeup_feed_404.outlet.pressure[t0].fix(101325.0)
    m.fs.wash_makeup_feed_404.outlet.mole_frac_comp[t0, "AS"].fix(eps)
    m.fs.wash_makeup_feed_404.outlet.mole_frac_comp[t0, "AA"].fix(major)
    m.fs.wash_makeup_feed_404.outlet.mole_frac_comp[t0, "ASA"].fix(eps)
    m.fs.wash_makeup_feed_404.outlet.mole_frac_comp[t0, "AcOH"].fix(eps)
    m.fs.wash_makeup_feed_404.outlet.mole_frac_comp[t0, "H2SO4"].fix(eps)

    m.fs.recycle_splitter_455.split_fraction[t0, "aa_recycle_out", "AA"].fix(0.88)
    m.fs.recycle_splitter_455.split_fraction[t0, "wash_recycle_out", "AA"].fix(0.07)
    m.fs.recycle_splitter_455.split_fraction[t0, "seed_solvent_out", "AA"].fix(0.05)

    m.fs.recycle_splitter_455.split_fraction[t0, "aa_recycle_out", "AcOH"].fix(0.001)
    m.fs.recycle_splitter_455.split_fraction[t0, "wash_recycle_out", "AcOH"].fix(0.60)
    m.fs.recycle_splitter_455.split_fraction[t0, "seed_solvent_out", "AcOH"].fix(0.399)

    m.fs.recycle_splitter_455.split_fraction[t0, "aa_recycle_out", "AS"].fix(0.001)
    m.fs.recycle_splitter_455.split_fraction[t0, "wash_recycle_out", "AS"].fix(0.4995)
    m.fs.recycle_splitter_455.split_fraction[t0, "seed_solvent_out", "AS"].fix(0.4995)

    m.fs.recycle_splitter_455.split_fraction[t0, "aa_recycle_out", "ASA"].fix(0.001)
    m.fs.recycle_splitter_455.split_fraction[t0, "wash_recycle_out", "ASA"].fix(0.4995)
    m.fs.recycle_splitter_455.split_fraction[t0, "seed_solvent_out", "ASA"].fix(0.4995)

    m.fs.recycle_splitter_455.split_fraction[t0, "aa_recycle_out", "H2SO4"].fix(0.001)
    m.fs.recycle_splitter_455.split_fraction[t0, "wash_recycle_out", "H2SO4"].fix(0.4995)
    m.fs.recycle_splitter_455.split_fraction[t0, "seed_solvent_out", "H2SO4"].fix(0.4995)

    m.fs.aa_recycle_min_aa_frac = Constraint(
        expr=m.fs.recycle_splitter_455.aa_recycle_out_state[t0].mole_frac_comp["AA"] >= 0.99
    )
    m.fs.wash_recycle_min_aa_frac = Constraint(
        expr=m.fs.recycle_splitter_455.wash_recycle_out_state[t0].mole_frac_comp["AA"] >= 0.90
    )
    m.fs.seed_solvent_min_aa_frac = Constraint(
        expr=m.fs.recycle_splitter_455.seed_solvent_out_state[t0].mole_frac_comp["AA"] >= 0.90
    )

    m.fs.final_product_purity = Expression(
        expr=m.fs.asa_product_bin_520.inlet.mole_frac_comp[t0, "ASA"]
    )

    m.fs.mw_kg_per_mol = Param(
        m.fs.thermo_params.component_list,
        initialize={
            "AS": 0.13812,
            "AA": 0.10209,
            "ASA": 0.18016,
            "AcOH": 0.06005,
            "H2SO4": 0.09808,
        },
        mutable=True,
    )

    m.fs.raw_material_price_per_kg = Param(
        m.fs.thermo_params.component_list,
        initialize={
            "AS": 2.50,
            "AA": 1.60,
            "ASA": 0.00,
            "AcOH": 0.00,
            "H2SO4": 0.35,
        },
        mutable=True,
    )
    m.fs.asa_market_value_per_kg = Param(initialize=5.00, mutable=True)
    m.fs.electricity_price_per_kwh = Param(initialize=0.08, mutable=True)
    m.fs.thermal_utility_price_per_kwh = Param(initialize=0.03, mutable=True)

    m.fs.raw_material_cost_per_s = Expression(
        expr=sum(
            feed.outlet.flow_mol[t0]
            * feed.outlet.mole_frac_comp[t0, comp]
            * m.fs.mw_kg_per_mol[comp]
            * m.fs.raw_material_price_per_kg[comp]
            for feed in (
                m.fs.sa_makeup_feed_100,
                m.fs.aa_makeup_feed_105,
                m.fs.h2so4_makeup_feed_145,
                m.fs.wash_makeup_feed_404,
            )
            for comp in m.fs.thermo_params.component_list
        )
    )

    m.fs.electricity_cost_per_s = Expression(
        expr=(
            sum(
                pyunits.convert(pump.work_mechanical[t0], to_units=pyunits.W)
                for pump in (
                    m.fs.aa_feed_pump_121,
                    m.fs.reactor_feed_pump_141,
                    m.fs.h2so4_metering_pump_151,
                )
            )
            / 3.6e6
        )
        * m.fs.electricity_price_per_kwh
    )

    m.fs.thermal_utility_cost_per_s = Expression(
        expr=(
            sum(
                abs(pyunits.convert(unit.heat_duty[t0], to_units=pyunits.W))
                for unit in (
                    m.fs.aa_inline_filter_122,
                    m.fs.aa_coriolis_ft_123,
                    m.fs.reactor_feed_ft_142,
                    m.fs.h2so4_flow_ft_152,
                    m.fs.overhead_condenser_441,
                    m.fs.fluidized_bed_dryer_510,
                )
            )
            / 3.6e6
        )
        * m.fs.thermal_utility_price_per_kwh
    )

    m.fs.utility_cost_per_s = Expression(
        expr=m.fs.electricity_cost_per_s + m.fs.thermal_utility_cost_per_s
    )
    m.fs.total_operating_cost_per_s = Expression(
        expr=m.fs.raw_material_cost_per_s + m.fs.utility_cost_per_s
    )

    m.fs.asa_product_mass_flow_kg_per_s = Expression(
        expr=m.fs.asa_product_bin_520.inlet.flow_mol[t0]
        * m.fs.asa_product_bin_520.inlet.mole_frac_comp[t0, "ASA"]
        * m.fs.mw_kg_per_mol["ASA"]
    )
    m.fs.asa_product_value_per_s = Expression(
        expr=m.fs.asa_product_mass_flow_kg_per_s * m.fs.asa_market_value_per_kg
    )
    m.fs.net_operating_margin_per_s = Expression(
        expr=m.fs.asa_product_value_per_s - m.fs.total_operating_cost_per_s
    )

    m.fs.aa_feed_pump_121.deltaP[t0].fix(2.0e5)
    m.fs.reactor_feed_pump_141.deltaP[t0].fix(2.0e5)
    m.fs.h2so4_metering_pump_151.deltaP[t0].fix(2.0e5)

    m.fs.asa_seed_takeoff_splitter_305.split_fraction[t0, "product_out"].fix(0.95)
    m.fs.seed_distribution_manifold_315.split_fraction[t0, "seed1_out"].fix(0.50)
    m.fs.pusher_centrifuge_splitter_410.split_fraction[t0, "wet_cake_out"].fix(0.50)
    m.fs.solvent_recovery_column_splitter_440.split_fraction[t0, "overhead_vapor_out"].fix(0.50)
    m.fs.reflux_drum_receiver_442.split_fraction[t0, "reflux_out"].fix(0.50)

    TransformationFactory("network.expand_arcs").apply_to(m)

    for arc in m.component_data_objects(Arc, descend_into=True):
        propagate_state(arc)

    return m


def report_recycle_aa_fractions(m: ConcreteModel, time_point=None):
    t = m.fs.time.first() if time_point is None else time_point
    aa_recycle = value(m.fs.recycle_splitter_455.aa_recycle_out_state[t].mole_frac_comp["AA"])
    wash_recycle = value(m.fs.recycle_splitter_455.wash_recycle_out_state[t].mole_frac_comp["AA"])
    seed_solvent = value(m.fs.recycle_splitter_455.seed_solvent_out_state[t].mole_frac_comp["AA"])

    print(f"AA mole fraction - aa_recycle_out: {aa_recycle:.6f}")
    print(f"AA mole fraction - wash_recycle_out: {wash_recycle:.6f}")
    print(f"AA mole fraction - seed_solvent_out: {seed_solvent:.6f}")

    return {
        "aa_recycle_out": aa_recycle,
        "wash_recycle_out": wash_recycle,
        "seed_solvent_out": seed_solvent,
    }


def report_operating_costs(m: ConcreteModel):
    print(f"Raw material cost ($/s): {value(m.fs.raw_material_cost_per_s):.6f}")
    print(f"Electricity cost ($/s): {value(m.fs.electricity_cost_per_s):.6f}")
    print(f"Thermal utility cost ($/s): {value(m.fs.thermal_utility_cost_per_s):.6f}")
    print(f"Total utility cost ($/s): {value(m.fs.utility_cost_per_s):.6f}")
    print(f"Total operating cost ($/s): {value(m.fs.total_operating_cost_per_s):.6f}")
    print(f"ASA product value ($/s): {value(m.fs.asa_product_value_per_s):.6f}")
    print(f"Net operating margin ($/s): {value(m.fs.net_operating_margin_per_s):.6f}")

    return {
        "raw_material_cost_per_s": value(m.fs.raw_material_cost_per_s),
        "electricity_cost_per_s": value(m.fs.electricity_cost_per_s),
        "thermal_utility_cost_per_s": value(m.fs.thermal_utility_cost_per_s),
        "utility_cost_per_s": value(m.fs.utility_cost_per_s),
        "total_operating_cost_per_s": value(m.fs.total_operating_cost_per_s),
        "asa_product_value_per_s": value(m.fs.asa_product_value_per_s),
        "net_operating_margin_per_s": value(m.fs.net_operating_margin_per_s),
    }


def initialize_model_sequential(m: ConcreteModel, run_unit_initialization: bool = False):
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 5

    def _unfix_port_states(unit):
        for port in unit.component_data_objects(Port, descend_into=False):
            for member in port.vars.values():
                if hasattr(member, "is_indexed") and member.is_indexed():
                    variables = member.values()
                else:
                    variables = [member]
                for variable in variables:
                    if hasattr(variable, "fixed") and variable.fixed:
                        variable.unfix()

    def _initialize_unit(unit):
        if not run_unit_initialization:
            return
        if isinstance(unit, (Feed, Product, StateJunction)):
            return
        if unit is m.fs.recycle_splitter_455:
            return
        if hasattr(unit, "initialize"):
            try:
                unit.initialize(outlvl=idaeslog.WARNING, solver="ipopt", hold_state=False)
                _unfix_port_states(unit)
            except TypeError:
                try:
                    unit.initialize(outlvl=idaeslog.WARNING)
                    _unfix_port_states(unit)
                except Exception as err:
                    print(f"Warning: initialization skipped for {unit.name}: {err}")
            except Exception as err:
                print(f"Warning: initialization skipped for {unit.name}: {err}")

    graph = seq.create_graph(m)
    order = seq.calculation_order(graph)

    for _ in range(seq.options.iterLim):
        for blocks in order:
            for unit in blocks:
                _initialize_unit(unit)
                for _, _, edge_data in graph.out_edges(unit, data=True):
                    arc = edge_data.get("arc")
                    if arc is None:
                        continue
                    try:
                        propagate_state(arc, overwrite_fixed=True)
                    except TypeError:
                        propagate_state(arc)


def solve_model(m: ConcreteModel):
    solver = get_solver()
    return solver.solve(m, tee=False)


if __name__ == "__main__":
    model = build_model()
    print("Built steady-state unit model scaffold from YAML node set.")
    print(f"Degrees of freedom: {degrees_of_freedom(model)}")
    initialize_model_sequential(model)
    results = solve_model(model)
    print(f"Solver termination: {results.solver.termination_condition}")
    print(f"Solver optimal: {check_optimal_termination(results)}")
    print(f"Final product purity (ASA mole fraction): {value(model.fs.final_product_purity):.6f}")
    report_recycle_aa_fractions(model)
    report_operating_costs(model)