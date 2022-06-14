using Pkg
using Debugger
Pkg.activate("..")

using Distributed
const WORKERS = 1 # Change to > 1 to use parallel
if nworkers() < WORKERS
    diff = (nprocs() == nworkers() ? WORKERS : WORKERS - nworkers())
    println("Adding $diff worker processes.")
    Distributed.addprocs(diff)
    # Make sure these workers also have an environment with PH installed
    @everywhere using Pkg
    for w in workers()
        @spawnat(w, Pkg.activate(".."))
    end
end

#-------------------------
@everywhere using ProgressiveHedging
@everywhere const PH = ProgressiveHedging
@everywhere using Ipopt
@everywhere using JuMP, Xpress, PowerSimulations, PowerSystems, PowerSystemCaseBuilder, TimeSeries
@everywhere const PSI = PowerSimulations
@everywhere const PSY = PowerSystems
@everywhere const PSB = PowerSystemCaseBuilder

using JuMP, Xpress, PowerSimulations, PowerSystems, PowerSystemCaseBuilder, InfrastructureSystems, Dates, TimeSeries
const PSI = PowerSimulations
const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder
solver = optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 1e-4)

#-------------------------
system = PSB.build_system(SIIPExampleSystems, "5_bus_hydro_uc_sys"; force_build =true);
remove_time_series!(system, DeterministicSingleTimeSeries)
for load in get_components(PowerLoad, system)
    ts = get_time_series(SingleTimeSeries, load, "max_active_power");
    ts_val = get_time_series_values(SingleTimeSeries, load, "max_active_power")
    time_stamps = get_initial_timestamp(ts):Hour(1):get_initial_timestamp(ts)+Hour(335)
    remove_time_series!(system, SingleTimeSeries, load, "max_active_power")
    data = TimeArray(time_stamps, vcat(ts_val, ts_val + rand(168)*0.2))
    forecast = SingleTimeSeries("max_active_power", data;)
    add_time_series!(system, load, forecast)
end

for hyd in get_components(HydroDispatch, system), label in ["max_active_power"]
    ts = get_time_series(SingleTimeSeries, hyd, label);
    ts_val = get_time_series_values(SingleTimeSeries, hyd, label)
    time_stamps = get_initial_timestamp(ts):Hour(1):get_initial_timestamp(ts)+Hour(335)
    remove_time_series!(system, SingleTimeSeries, hyd, label)
    data = TimeArray(time_stamps, vcat(ts_val, ts_val + rand(168)*0.01))
    forecast = SingleTimeSeries(label, data;)
    add_time_series!(system, hyd, forecast)
end

labels = ["max_active_power", "inflow"]
for hyd in get_components(HydroEnergyReservoir, system), label in labels
    ts = get_time_series(SingleTimeSeries, hyd, label);
    ts_val = get_time_series_values(SingleTimeSeries, hyd, label)
    time_stamps = get_initial_timestamp(ts):Hour(1):get_initial_timestamp(ts)+Hour(335)
    remove_time_series!(system, SingleTimeSeries, hyd, label)
    data = TimeArray(time_stamps, vcat(ts_val, ts_val + rand(168)*0.01))
    forecast = SingleTimeSeries(label, data;)
    add_time_series!(system, hyd, forecast)
end
PSY.transform_single_time_series!(system, 168, Hour(168))

new_cost= Dict(
    "Alta" => ThreePartCost((0.0, 14.0), 0.0, 4.0, 2.0),
    "Park City" => ThreePartCost((0.0, 15.0), 0.0, 1.5, 0.75),
    "Solitude" => ThreePartCost((0.0, 30.0), 0.0, 3.0, 1.5),
    "Sundance" => ThreePartCost((0.0, 40.0), 0.0, 4.0, 2.0),
    "Brighton" => ThreePartCost((0.0, 10.0), 0.0, 1.5, 0.75),
    )
for th in get_components(ThermalStandard, system)
    set_operation_cost!(th, new_cost[th.name])
end

print(">>> Made it here...")
PSY.to_json(system, "5bus_hydro_2week.json"; force=true)
