using PowerSystems
using PowerSimulations
using Dates
using Logging
using DataFrames
logger = configure_logging(console_level = Logging.Info)
const PSI = PowerSimulations
const PSY = PowerSystems
using DrWatson

using JuMP
using Xpress
solver = optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 1e-5, "OUTPUTLOG" => 1, "PRESOLVE" => 0,  "MAXTIME" => 300)
#using Cbc
#solver = optimizer_with_attributes(Cbc.Optimizer, "MIPRELSTOP" => 1e-5, "OUTPUTLOG" => 1, "PRESOLVE" => 0, "MAXTIME" => 300)

### Parsing Args
sys_name = ARGS[1]
interval = parse(Int, ARGS[2])
horizon = parse(Int, ARGS[3])
steps = parse(Int, ARGS[4])
battery = parse(Bool, ARGS[5])


params = (
    originalsim = replace(split(sys_name, "/")[end-1], "=" => "-"),
    interval = interval,
    horizon = horizon,
    steps = steps,
    battery = battery
)
println(params)
name = savename(params)
### Simulation Setup

template_uc = template_unit_commitment(network = NetworkModel(DCPPowerModel, duals = [NodalBalanceActiveConstraint], use_slacks = true))
set_device_model!(template_uc, ThermalStandard, ThermalStandardUnitCommitment)
if battery
    set_device_model!(template_uc, DeviceModel(GenericBattery, BookKeeping; duals=[EnergyBalanceConstraint]))
end
set_device_model!(template_uc, MonitoredLine, StaticBranchBounds)
set_device_model!(template_uc, Line, StaticBranch)
set_device_model!(template_uc, Transformer2W, StaticBranch)
set_device_model!(template_uc, TapTransformer, StaticBranch)
set_device_model!(template_uc, HVDCLine, HVDCLossless)

sys = System(datadir(sys_name); time_series_directory = "/tmp/scratch")
PSY.transform_single_time_series!(sys, horizon, Hour(interval))

models = SimulationModels(
    decision_models = [
        DecisionModel(
            template_uc,
            sys,
            name = "UC",
            optimizer = solver,
            optimizer_solve_log_print = true,
        ),
    ],
)

sequence = SimulationSequence(
    models = models,
    ini_cond_chronology = InterProblemChronology()
)

sim = Simulation(
    name = name,
    steps = steps,
    models = models,
    sequence = sequence,
    simulation_folder = datadir("powersim_results")
)
build!(sim)
execute!(sim)
