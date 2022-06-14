using DrWatson

using PowerSystems
using Logging
using Unitful
const PSY = PowerSystems
logger = PSY.configure_logging(console_level = Logging.Info)

function copy_component(re::PSY.RenewableDispatch, sys, bus_name, name)
    bus = PSY.get_component(PSY.Bus, sys, bus_name)
    device = PSY.RenewableDispatch(;
        name = name,
        available=true,
        bus=bus,
        active_power=re.active_power,
        reactive_power=re.reactive_power,
        rating=re.rating,
        prime_mover=re.prime_mover,
        reactive_power_limits=re.reactive_power_limits,
        power_factor=re.power_factor,
        operation_cost=re.operation_cost,
        base_power=re.base_power,
    )
    return device
end

line_susceptance = Dict(
    "line_ab" => 0.00356,
    "line_ae" => 0.00329,
    "line_ad" => 0.01563,
    "line_bc" => 0.00926,
    "line_cd" => 0.00337,
    "line_de" => 0.00337,
)


function load_sys(fivebusdir)
    rawsys = PSY.PowerSystemTableData(
        fivebusdir,
        100.0,
        joinpath(fivebusdir, "user_descriptors.yaml"),
        generator_mapping_file = joinpath(fivebusdir, "generator_mapping.yaml"),
    )
    sys = PSY.System(rawsys; time_series_directory=datadir("timeseries"))

    for line in PSY.get_components(PSY.Line, sys)
        b = line_susceptance[line.name]
        PSY.set_b!(line, (from=b, to=b))
    end
    return sys
end

function go_30_to_45!(sys)
    re = PSY.get_component(PSY.RenewableGen, sys, "SolarPV2")
    re_3 = copy_component(re, sys, "node_b", "SolarPV3")
    PSY.add_component!(sys, re_3)
    PSY.copy_time_series!(re_3, re)

    wind = PSY.get_component(PSY.RenewableGen, sys, "Wind")
    wind_2 = copy_component(re, sys, "node_b", "Wind2")
    PSY.set_rating!(wind_2, wind.rating/2)
    PSY.add_component!(sys, wind_2)
    PSY.copy_time_series!(wind_2, wind)

    for name in ["Park_City", "Alta", "Sundance"]
        PSY.remove_component!(sys, PSY.get_component(PSY.ThermalGen, sys, name))
    end

    return sys
end

function go_45_to_60!(sys)
    re = PSY.get_component(PSY.RenewableGen, sys, "SolarPV1")
    re_4 = copy_component(re, sys, "node_e", "SolarPV4")
    PSY.add_component!(sys, re_4)
    PSY.copy_time_series!(re_4, re)

    wind = PSY.get_component(PSY.RenewableGen, sys, "Wind")
    wind_3 = copy_component(re, sys, "node_d", "Wind3")
    PSY.set_rating!(wind_3, wind.rating/2)
    PSY.add_component!(sys, wind_3)
    PSY.copy_time_series!(wind_3, wind)
    return sys
end

function make_battery(
    ; name,
    bus,
    hours,
    base_power=200u"MW",
    efficiency = 0.85
)
    # Everything is in units of base_power
    return PSY.GenericBattery(;
        name = name,
        available = true,
        bus = bus,
        prime_mover = PSY.PrimeMovers.BA,
        initial_energy = ustrip(u"hr", hours / 2),
        state_of_charge_limits = (min = 0.0, max = ustrip(u"hr", hours)),
        rating = 1.0,
        active_power = 0.0,
        input_active_power_limits = (min = 0.0, max = 1.0),
        output_active_power_limits = (min = 0.0, max = 1.0),
        efficiency = (in = efficiency, out = 1.0),
        reactive_power = 0.0,
        reactive_power_limits = nothing,
        base_power = ustrip(u"MW", base_power)
    )
end

function add_storage!(sys, storages)
    for hours in storages
        if hours == 24 || hours == 10
            node = "node_b"  # second largest demand
        else
            node = "node_a"  # no demand there
        end
        bat = make_battery(
            name="Storage_$hours", bus=PSY.get_component(PSY.Bus, sys, node),
            hours=hours * u"hr"
        )
        PSY.add_component!(sys, bat)
    end
    return sys
end

for storage_hours in [[2], [4], [10], [24], [4, 10], [2, 10], [4, 24], [2, 10]]
    @info "Loading rawsys from for no storage"
    sys = load_sys(datadir("initial_system", "inputs_no_bat"))
    add_storage!(sys, storage_hours)
    params = (pv=30, wind=30, storagehours=join(storage_hours, "+"))
    new_sys_path = datadir("new_systems", savename(params), "sys.json")
    @info "Saving $params system at $new_sys_path"
    PSY.to_json(sys, new_sys_path; force=true)

    params = (pv=45, wind=45, storagehours=join(storage_hours, "+"))
    @info "Creating $params system"
    go_30_to_45!(sys)
    new_sys_path = datadir("new_systems", savename(params), "sys.json")
    @info "Saving $new_sys_path"
    PSY.to_json(sys, new_sys_path; force = true)

    params = (pv=60, wind=60, storagehours=join(storage_hours, "+"))
    @info "Creating $params system"
    go_45_to_60!(sys)
    new_sys_path = datadir("new_systems", savename(params), "sys.json")
    @info "Saving $new_sys_path"
    PSY.to_json(sys, new_sys_path; force=true)
end
