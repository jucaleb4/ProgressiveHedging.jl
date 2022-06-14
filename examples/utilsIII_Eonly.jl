struct HydroDispatchPH <: PSI.AbstractHydroReservoirFormulation end
struct HydroDispatchReservoirPH <: PSI.AbstractHydroReservoirFormulation end

abstract type PHProblem <: PSI.PowerSimulationsOperationsProblem end
abstract type PHReservoirProblem <: PSI.PowerSimulationsOperationsProblem end
abstract type PHStochasticReservoirProblem <: PSI.PowerSimulationsOperationsProblem end

struct EnergyBudget <: PSI.VariableType end
PSI.make_variable_name(::Type{EnergyBudget}, ::Type{T}) where {T <: PSY.Device} = encode_symbol(T, "Pb")
const POWER_BUDGET = "Pb"

struct EnergyTarget <: PSI.VariableType end
PSI.make_variable_name(::Type{EnergyTarget}, ::Type{T}) where {T <: PSY.Device} = encode_symbol(T, "Eb")
const ENERGY_BUDGET = "Eb"

function PSI.construct_device!(
    optimization_container::PSI.OptimizationContainer,
    sys::PSY.System,
    model::PSI.DeviceModel{H, HydroDispatchPH},
    ::Type{S},
) where {
    H <: PSY.HydroGen,
    S <: PSI.PM.AbstractActivePowerModel,
}
    devices = PSI.get_available_components(H, sys)

    if !PSI.validate_available_devices(H, devices)
        return
    end

    # Variables
    PSI.add_variables!(optimization_container, PSI.ActivePowerVariable, devices, HydroDispatchPH())

    # Constraints
    PSI.add_constraints!(
        optimization_container,
        PSI.RangeConstraint,
        PSI.ActivePowerVariable,
        devices,
        model,
        S,
        PSI.get_feedforward(model),
    )
    PSI.feedforward!(optimization_container, devices, model, PSI.get_feedforward(model))

    # Cost Function
    PSI.cost_function!(optimization_container, devices, model, S, nothing)

    return
end

function PSI.construct_device!(
    optimization_container::PSI.OptimizationContainer,
    sys::PSY.System,
    model::PSI.DeviceModel{H, HydroDispatchReservoirPH},
    ::Type{S},
) where {H <: PSY.HydroEnergyReservoir, S <: PSI.PM.AbstractActivePowerModel}
    devices = PSI.get_available_components(H, sys)

    if !PSI.validate_available_devices(H, devices)
        return
    end

    # Variables

    PSI.add_variables!(
        optimization_container,
        PSI.ActivePowerVariable,
        devices,
        HydroDispatchReservoirPH(),
    )
    PSI.add_variables!(
        optimization_container,
        PSI.EnergyVariable,
        devices,
        HydroDispatchReservoirPH(),
    )
    PSI.add_variables!(
        optimization_container,
        PSI.SpillageVariable,
        devices,
        HydroDispatchReservoirPH(),
    )

    # Constraints
    PSI.add_constraints!(
        optimization_container,
        PSI.RangeConstraint,
        PSI.ActivePowerVariable,
        devices,
        model,
        S,
        PSI.get_feedforward(model),
    )

    # Initial Conditions
    PSI.storage_energy_initial_condition!(
        optimization_container,
        devices,
        HydroDispatchReservoirPH(),
    )
    # Energy Balance Constraint
    PSI.add_constraints!(
        optimization_container,
        PSI.EnergyBalanceConstraint,
        PSI.EnergyVariable,
        devices,
        model,
        S,
        PSI.get_feedforward(model),
    )
    PSI.feedforward!(optimization_container, devices, model, PSI.get_feedforward(model))

    # Cost Function
    PSI.cost_function!(optimization_container, devices, model, S, nothing)

    return
end


### Build PH OperationProblem

function build_operations_problem(
        system::PSY.System,
        initial_time,
        solver = JuMP.optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 1e-6),
        network = PSI.CopperPlatePowerModel,
        directory ="./simulation_folder/",
        template = nothing
    )

    if isnothing(template)
        template = PSI.OperationsProblemTemplate(network)
        PSI.set_device_model!(template, PSY.ThermalStandard, PSI.ThermalDispatchNoMin)
        PSI.set_device_model!(template, PSY.PowerLoad, PSI.StaticPowerLoad)
        PSI.set_device_model!(template, PSY.HydroEnergyReservoir, HydroDispatchReservoirPH)
        PSI.set_device_model!(template, PSY.RenewableDispatch, RenewableFullDispatch)
        PSI.set_device_model!(template, PSY.RenewableFix, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.HydroDispatch, PSI.FixedOutput)
    end
    
    problem = PSI.OperationsProblem(
        template, 
        system, 
        optimizer = solver,
        balance_slack_variables = true,
        initial_time = initial_time,
        PTDF = PSY.PTDF(system)
    )
    PSI.build!(problem, output_dir = directory)
    return problem
end
    


function populate_ext!(problem::PSI.OperationsProblem, scenario_no, inflow_data, no_node_2nd, no_node_3rd)
    problem.ext["scenario_no"] = scenario_no
    problem.ext["inflow_data"] = inflow_data
    problem.ext["no_node_2nd"] = no_node_2nd
    problem.ext["no_node_3rd"] = no_node_3rd

    return
end


function populate_PH_model!(problem)
    scenario_no = problem.ext["scenario_no"]
    inflow_data = problem.ext["inflow_data"]
    no_node_2nd = problem.ext["no_node_2nd"]
    no_node_3rd = problem.ext["no_node_3rd"]
    scen_no_2 =  Int(round((scenario_no+1)/no_node_3rd, RoundUp))
    if scenario_no <= (no_node_3rd-1)
        scen_no_3 = scenario_no +1
    else
        scen_no_3 =  Int((scenario_no+1) - ((scen_no_2-1)*no_node_3rd))
    end

    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    var_h = PSI.get_variable(optimization_container, :P__HydroEnergyReservoir)
    var_sp = PSI.get_variable(optimization_container, :Sp__HydroEnergyReservoir)
    var_e = PSI.get_variable(optimization_container, :E__HydroEnergyReservoir)
    energy_balance = PSI.get_constraint(optimization_container, :energy_capacity__HydroEnergyReservoir)
    energy_var = PSI.add_var_container!(
        optimization_container,
        :Eb_HydroEnergyReservoir,
        collect(keys(inflow_data)),
        1:no_node_2nd,
    )
    scen_energy_var = PSI.add_var_container!(
        optimization_container,
        :Ebs_HydroEnergyReservoir,
        collect(keys(inflow_data)),
        1:no_node_3rd,
    )
    
    jump_model = PSI.get_jump_model(optimization_container)
    for (name, inflow) in inflow_data
        for sec_1st in 1:no_node_2nd
            energy_var[name, sec_1st] = JuMP.@variable(
                jump_model,
                base_name = "Eb_HydroEnergyReservoir_{$(name)_{$(sec_1st)}}",
            )
        end
        for sec_2nd in 1:no_node_3rd
            scen_energy_var[name, sec_2nd] = JuMP.@variable(
                jump_model,
                base_name = "Ebs_HydroEnergyReservoir_{$(name)_{$(scen_no_2), $sec_2nd}}",
            )
        end
        if !(scenario_no <= no_node_3rd - 1)
            JuMP.delete(jump_model, energy_balance[name, 1])
            energy_balance[name, 1] = JuMP.@constraint(jump_model, 
                var_e[name, 1] == energy_var[name, scen_no_2-1] + inflow[(scen_no_2, scen_no_3)] -  var_sp[name, 1] - var_h[name, 1]
            )
        end

        JuMP.@constraint(jump_model, var_e[name, 24] == scen_energy_var[name, scen_no_3])
        JuMP.@constraint(jump_model, energy_var[name, scen_no_2] == sum(scen_energy_var[name, i]*(1/no_node_3rd) for i in 1:no_node_3rd))
    end
    return
end

function get_subproblem(problem)
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    
    first_stage_var =  JuMP.VariableRef[]
    secound_stage_var = JuMP.VariableRef[]
    third_stage_var = JuMP.VariableRef[]
    for (name, cont) in PSI.get_variables(problem) 
        if (name == :Eb_HydroEnergyReservoir)
            push!(first_stage_var, cont...)
        elseif (name == :Ebs_HydroEnergyReservoir)
            push!(secound_stage_var, cont...)
        else
            push!(third_stage_var, cont...)
        end
    end
    variable_map = Dict{PH.StageID, Vector{JuMP.VariableRef}}(
        [
            PH.stid(1) => first_stage_var,
            PH.stid(2) => secound_stage_var,
            PH.stid(3) => third_stage_var,
        ]
    )
    return jump_model, variable_map
end


function create_model_stochastic(
    scenario_id::PH.ScenarioID,
    problem;
    kwargs...
)
    inflow_data = problem.ext["inflow_data"]
    no_node_2nd = problem.ext["no_node_2nd"]
    no_node_3rd = problem.ext["no_node_3rd"]
    total_scenarios = no_node_2nd*no_node_3rd
    sim_initial_time = problem.internal.optimization_container.settings.initial_time.x
    initial_time = collect(sim_initial_time:Day(1):sim_initial_time+Day(total_scenarios))
    scenario_init_time = Dict(ix => init for (ix, init) in enumerate(initial_time));
    scen_no_2 =  Int(round((PH.value(scenario_id) + 1)/no_node_3rd, RoundUp))
    if PH.value(scenario_id) <= (no_node_3rd-1)
        scen_no_3 = PH.value(scenario_id) +1
    else
        scen_no_3 =  Int((PH.value(scenario_id) + 1) - (scen_no_2-1)*no_node_3rd)
    end
    systems = problem.ext["ph_systems"] 
    template = get(problem.ext, "ph_template", nothing)
    trans_template = isnothing(template) ? PSI.CopperPlatePowerModel : template.transmission
    sub_problem = build_operations_problem(
        PSY.System(joinpath("/home/sdalvi/work/hydro_simulations/data",systems[scen_no_3]), time_series_directory = "/tmp/scratch"),
        scenario_init_time[scen_no_2],
        JuMP.optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 1e-5, "MAXTIME" => 300, "MAXMEMORYSOFT" => 5000),
        trans_template,
        "./simulation_folder/",
        template,
    )
    populate_ext!(sub_problem, PH.value(scenario_id), inflow_data, no_node_2nd, no_node_3rd)
    populate_PH_model!(sub_problem)
    jump_model, variable_map = get_subproblem(sub_problem)
    return JuMPSubproblem(jump_model, scenario_id, variable_map)
end


### Custom Stage


add_to_ext!(p::PSI.OperationsProblem, key, data) = p.ext[key] = data

function generate_inflow_initial_value(problem)
    sim_initial_time = problem.internal.optimization_container.settings.initial_time.x
    sys = problem.sys
    systems = problem.ext["ph_systems"] 
    no_node_2nd = problem.ext["no_node_2nd"]
    no_node_3rd = problem.ext["no_node_3rd"]
    inflow_data = Dict()
    names = PSY.get_name.(PSY.get_components(HydroEnergyReservoir, sys))
    for h in  names
        inflow_data[h] = Dict()
        for sec in 1:no_node_2nd, sec2 in 1:no_node_3rd
            device = PSY.get_component(HydroEnergyReservoir, sys, h)
            ts = get_time_series_values(
                SingleTimeSeries, 
                device,
                "max_active_power", 
                start_time = sim_initial_time + Day(sec-1), 
                len = 24
            )
            mul = PSY.get_rating(device)
            inflow_data[h][(sec, sec2)] = ts[1] .* mul
        end
    end
  return inflow_data
end


function PSI.OperationsProblem(
    ::Type{PHStochasticReservoirProblem},
    template::OperationsProblemTemplate,
    sys::PSY.System,
    no_node_2nd,
    no_node_3rd,
    ph_systems,
    ph_template,
    jump_model::Union{Nothing, JuMP.AbstractModel} = nothing;
    kwargs...,
)
    problem = OperationsProblem{PHStochasticReservoirProblem}(template, sys, jump_model; kwargs...)

    add_to_ext!(problem, "no_node_2nd", no_node_2nd)
    add_to_ext!(problem, "no_node_3rd", no_node_3rd)
    add_to_ext!(problem, "ph_systems", ph_systems)
    add_to_ext!(problem, "ph_template", ph_template)

    return problem
end


function PSI._build!(problem::PSI.OperationsProblem{PHStochasticReservoirProblem}, serialize::Bool)
    ## TODO: Add Container for Storaging Budget
    PSI.build_pre_step!(problem)
    system = get_system(problem)
    inflow_data = generate_inflow_initial_value(problem)
    add_to_ext!(problem, "inflow_data", inflow_data)
    container = problem.internal.optimization_container
    devices = get_available_components(HydroEnergyReservoir, system)
    ## Need to fix ph_problem.internal.optimization_container.time_steps
    no_node_2nd = problem.ext["no_node_2nd"]
    no_node_3rd = problem.ext["no_node_3rd"]
    total_scenarios = no_node_2nd*no_node_3rd

    global_energy = PSI.container_spec(
        JuMP.VariableRef,
        [PSY.get_name(d) for d in devices],
        1:no_node_2nd,
    )

    PSI.assign_variable!(container, "Eb", PSY.HydroEnergyReservoir, global_energy)
    PSI.set_status!(problem, BuildStatus.BUILT)
    return BuildStatus.BUILT
end

function build_3_stage_tree(no_nodes_2nd, no_nodes_3rd)
    scen_tree = ScenarioTree()
    for s in 1:no_nodes_2nd
        node = add_node(scen_tree, root(scen_tree))
        for ss in 1:no_nodes_3rd
          add_leaf(scen_tree, node, 1/(no_nodes_2nd*no_nodes_3rd))  
        end
    end
    return scen_tree
end

function add_ph_results!(problem)
    solution = problem.ext["solution"]
    no_node_2nd = problem.ext["no_node_2nd"]
    inflow_data = problem.ext["inflow_data"]
    optimization_container = problem.internal.optimization_container
    jump_model = PSI.get_jump_model(optimization_container)
    g_energy = PSI.get_variable(optimization_container, Symbol("Eb__HydroEnergyReservoir"))
    names = collect(keys(inflow_data))

    Eb_mapping = Dict("Eb_HydroEnergyReservoir_{$(name)_{$(sec_1st)}}" => (name, sec_1st) for sec_1st in 1:no_node_2nd, name in names)
    Eb_variable_names = collect(keys(Eb_mapping))
    results = filter(row -> row.variable ∈ Eb_variable_names, solution)
    for row in eachrow(results)
        name, t = Eb_mapping[row.variable]
        g_energy[name, t] =  JuMP.@variable(
                jump_model,
                base_name = "Eb_HydroEnergyReservoir_{$(name)_{$(t)}",
            )
        JuMP.fix(g_energy[name, t], row.value)
    end
    JuMP.optimize!(jump_model)
end
  
  
function PSI.solve_impl(problem::PSI.OperationsProblem{PHStochasticReservoirProblem}; optimizer = nothing)

    no_node_2nd = problem.ext["no_node_2nd"]
    no_node_3rd = problem.ext["no_node_3rd"]
    st = build_3_stage_tree(no_node_2nd, no_node_3rd)
    (n, err, rerr, obj, soln, phd) = PH.solve(st,
                                              create_model_stochastic,
                                              PH.ScalarPenaltyParameter(45.0), problem,
                                              atol=1e-3, rtol=1e-4, max_iter=150,
                                              report=1, # print residual info every 10 iterations
                                              )
    println("iterations = ", n)
    println("error = ", err)
    println("relative error = ", rerr)
    println("Objective = ", obj)
    add_to_ext!(problem, "iterations", n)
    add_to_ext!(problem, "error", err)
    add_to_ext!(problem, "relative error ", rerr)
    add_to_ext!(problem, "Objective", obj)
    add_to_ext!(problem, "solution", soln)
    add_to_ext!(problem, "PHObject", phd)
    add_ph_results!(problem)
    PSI.set_run_status!(problem, RunStatus.SUCCESSFUL)
    return RunStatus.SUCCESSFUL
end


function PSI.write_problem_results!(
    step::Int,
    problem::PSI.OperationsProblem{T},
    start_time::Dates.DateTime,
    store::PSI.SimulationStore,
    exports,
) where {T<:Union{PHReservoirProblem, PHStochasticReservoirProblem, PHProblem}}

    PSI.write_model_results!(store, problem, start_time; exports = exports)
    return
end


struct BudgetLimitFF <: PSI.AbstractAffectFeedForward
    variable_source_problem::Symbol
    affected_variables::Vector{Symbol}
    cache::Union{Nothing, Type{<:PSI.AbstractCache}}
    function BudgetLimitFF(
        variable_source_problem::AbstractString,
        affected_variables::Vector{<:AbstractString},
        cache::Union{Nothing, Type{<:PSI.AbstractCache}},
    )
        new(Symbol(variable_source_problem), Symbol.(affected_variables), cache)
    end
end

function BudgetLimitFF(; variable_source_problem, affected_variables)
    return BudgetLimitFF(variable_source_problem, affected_variables, nothing)
end

PSI.get_variable_source_problem(p::BudgetLimitFF) = p.variable_source_problem

function PSI.feedforward!(
    optimization_container::PSI.OptimizationContainer,
    devices::PSI.IS.FlattenIteratorWrapper{T},
    ::PSI.DeviceModel{T, <:PSI.AbstractDeviceFormulation},
    ff_model::BudgetLimitFF,
) where {T <: PSY.StaticInjection}
    for prefix in PSI.get_affected_variables(ff_model)
        var_name = PSI.make_variable_name(prefix, T)
        source_var_name = PSI.make_variable_name(PSI.get_variable_source_problem(ff_model), T)
        parameter_ref = PSI.UpdateRef{JuMP.VariableRef}(source_var_name)
        budget_limit_ff(
            optimization_container,
            PSI.make_constraint_name(FEEDFORWARD_INTEGRAL_LIMIT, T),
            parameter_ref,
            var_name,
        )
    end
end

function budget_limit_ff(
    optimization_container::PSI.OptimizationContainer,
    cons_name::Symbol,
    param_reference::PSI.UpdateRef,
    var_name::Symbol,
)
    time_steps = PSI.model_time_steps(optimization_container)
    variable = PSI.get_variable(optimization_container, var_name)

    axes = JuMP.axes(variable)
    set_name = axes[1]

    @assert axes[2] == time_steps
    container_ub = PSI.add_param_container!(optimization_container, param_reference, set_name)
    param_ub = PSI.get_parameter_array(container_ub)
    multiplier_ub = PSI.get_multiplier_array(container_ub)
    con_ub = PSI.add_cons_container!(optimization_container, cons_name, set_name)

    for name in axes[1]
        value = JuMP.upper_bound(variable[name, 1])
        param_ub[name] = PSI.add_parameter(optimization_container.JuMPmodel, value)
        # default set to 1.0, as this implementation doesn't use multiplier
        multiplier_ub[name] = 1.0
        con_ub[name] = JuMP.@constraint(
            optimization_container.JuMPmodel,
            sum(variable[name, t] for t in time_steps)  <=
            param_ub[name] * multiplier_ub[name]
        )
    end
end

function get_hourly_budget_ts(hyd, utilization, length)
    ts_value = (get_active_power_limits(hyd).max * utilization) / get_storage_capacity(hyd)
    ts_vector = ones(length)*ts_value
    return ts_vector
end

