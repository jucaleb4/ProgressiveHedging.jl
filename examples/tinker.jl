using JuMP, Xpress, PowerSimulations, PowerSystems, PowerSystemCaseBuilder, InfrastructureSystems, Dates, TimeSeries
const PSI = PowerSimulations
const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder
solver = optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 1e-4)
using Debugger

function build()
    system = PSB.build_system(SIIPExampleSystems, "5_bus_hydro_uc_sys_with_targets"; force_build =true, has_reserve=false);
    for serv in get_components(Service, system)
        remove_component!(system, serv)
    end
    remove_time_series!(system, DeterministicSingleTimeSeries)
    remove_time_series!(system, Deterministic)
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
        @show time_stamps
        remove_time_series!(system, SingleTimeSeries, hyd, label)
        data = TimeArray(time_stamps, vcat(ts_val, ts_val + rand(168)*0.01))
        forecast = SingleTimeSeries(label, data;)
        add_time_series!(system, hyd, forecast)
    end
    @enter PSY.transform_single_time_series!(system, 24, Hour(24))
end


build()
