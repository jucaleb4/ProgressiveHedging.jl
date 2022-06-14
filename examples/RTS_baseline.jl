using Pkg
Pkg.activate("..")
empty!(DEPOT_PATH)
push!(DEPOT_PATH, "/lustre/eaglefs/projects/pvb/cju/.julia")
# using Debugger

# Distributed setup
using ProgressiveHedging
const PH = ProgressiveHedging
using Ipopt
using JuMP, Xpress, PowerSimulations, PowerSystems, PowerSystemCaseBuilder, TimeSeries
const PSI = PowerSimulations
const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder

# Local setup
using JuMP, Xpress, PowerSimulations, PowerSystems, PowerSystemCaseBuilder, InfrastructureSystems, Dates, TimeSeries
const PSI = PowerSimulations
const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder

solver = optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 1e-4, "OUTPUTLOG" => 1, "MAXTIME" => 500,)

sys_name = "/lustre/eaglefs/projects/pvb/cju/gen_systems/data/rts_with_battery_060822_sys.json"
system = PSY.System(sys_name, time_series_directory = "/tmp/scratch")
PSY.transform_single_time_series!(system, 672, Hour(168))
# PSY.transform_single_time_series!(system, #Horizon, #Interval)

template_ed = ProblemTemplate(PowerSimulations.DCPPowerModel)
# Device models
set_device_model!(template_ed, PowerSystemCaseBuilder.Line, StaticBranch)
set_device_model!(template_ed, PowerSystemCaseBuilder.Transformer2W, StaticBranch)
set_device_model!(template_ed, PowerSystemCaseBuilder.TapTransformer, StaticBranch)

# Other systems
set_device_model!(template_ed, PowerSystems.ThermalStandard, ThermalBasicUnitCommitment)
set_device_model!(template_ed, PowerSystems.RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_ed, PowerSystems.PowerLoad, StaticPowerLoad)
set_device_model!(template_ed, PowerSystems.HydroDispatch, FixedOutput)
battery_device_model = DeviceModel(
                           GenericBattery,
                           BookKeeping;
                           attributes=Dict{String, Any}("reservation" => false),
                           )
set_device_model!(template_ed, battery_device_model)
# set_device_model!(template_ed, PowerSystems.GenericBattery, BookKeeping)
set_device_model!(template_ed, PowerSystems.RenewableFix, FixedOutput)

# Ancillary services
# set_service_model!(template_ed, VariableReserve{ReserveUp}, RangeReserve)
# set_service_model!(template_ed, VariableReserve{ReserveDown}, RangeReserve)

problem = DecisionModel(
    template_ed, 
    system,
    horizon = 24*15,
    name = "SubProblem",
    optimizer = solver, 
    warm_start = false,
    export_pwl_vars = true,
    initialize_model = false,
    initial_time = DateTime("2020-01-01T00:00:00"),
)

directory ="./simulation_folder/"
build!(problem, output_dir = directory)
solve!(problem)

res = ProblemResults(problem)
opt_stats = get_optimizer_stats(res)
println("===============")
objval = get_objective_value(res)
println("Objective value: ", objval)
# read_variables(res)
