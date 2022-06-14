struct HydroDispatchPH <: PSI.AbstractHydroReservoirFormulation end
struct HydroDispatchReservoirPH <: PSI.AbstractHydroReservoirFormulation end

abstract type PHProblem <: PSI.PowerSimulationsOperationsProblem end
abstract type PHReservoirProblem <: PSI.PowerSimulationsOperationsProblem end

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
        PSI.set_device_model!(template, PSY.HydroEnergyReservoir, HydroDispatchPH)
        PSI.set_device_model!(template, PSY.RenewableDispatch, RenewableFullDispatch)
        PSI.set_device_model!(template, PSY.RenewableFix, PSI.FixedOutput)
        PSI.set_device_model!(template, PSY.HydroDispatch, PSI.FixedOutput)
    end
    
    problem = PSI.OperationsProblem(
        template, 
        system, 
        optimizer = solver,
        balance_slack_variables = true,
        warm_start = true,
        export_pwl_vars = true,
        initial_time = initial_time 
    )
    PSI.build!(problem, output_dir = directory)
    return problem
end
    
# function populate_ext!(problem::PSI.OperationsProblem, scenario_no, budget_data, total_scenarios)
function populate_ext!(problem::PSI.OperationsProblem, scenario_no, inflow_data, initial_condtions, total_scenarios)
    problem.ext["scenario_no"] = scenario_no
#     problem.ext["budget_data"] = budget_data
    problem.ext["total_scenarios"] = total_scenarios
    problem.ext["initial_condtions"] = initial_condtions
    problem.ext["inflow_data"] = inflow_data

    return
end

# function populate_PH_model!(problem::PSI.OperationsProblem)
#     no_scenarios = problem.ext["total_scenarios"]
#     scenario_no = problem.ext["scenario_no"]
#     budget_data = problem.ext["budget_data"]
    
#     optimization_container = PSI.get_optimization_container(problem)
#     time_steps = PSI.model_time_steps(optimization_container)
#     var_h = PSI.get_variable(optimization_container, :P__HydroEnergyReservoir)
    
#     constraint_budget = PSI.add_cons_container!(
#         optimization_container,
#         :Hydro_Budget_HydroEnergyReservoir,
#         collect(keys(budget_data)),
#     )
#     constraint_scenarios = PSI.add_cons_container!(
#         optimization_container,
#         :PH_scenario_HydroEnergyReservoir,
#         collect(keys(budget_data)),
#     )
#     budget_var = PSI.add_var_container!(
#         optimization_container,
#         :Pb_HydroEnergyReservoir,
#         collect(keys(budget_data)),
#         1:no_scenarios,
#     )
#     jump_model = PSI.get_jump_model(optimization_container)
#     for (name, budget) in budget_data
#         for sec in 1:no_scenarios
#             budget_var[name, sec] = JuMP.@variable(
#                 optimization_container.JuMPmodel,
#                 base_name = "Pb_HydroEnergyReservoir_{$(name)_{$(sec)}}",
#             )
#         end
#         constraint_budget[name] =
#             JuMP.@constraint(jump_model, sum(var_h[name, t] for t in time_steps) <= budget_var[name, scenario_no+1])
#         constraint_scenarios[name] = 
#             JuMP.@constraint(jump_model, sum(budget_var[name, sec] for sec in 1:no_scenarios ) == budget)
        
#     end
#     return
# end   

function populate_PH_model!(problem::PSI.OperationsProblem)
    no_scenarios = problem.ext["total_scenarios"]
    scenario_no = problem.ext["scenario_no"]
#     budget_data = problem.ext["budget_data"]
    inflow_data = problem.ext["inflow_data"]
    initial_condtions = problem.ext["initial_condtions"]
    
    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    var_h = PSI.get_variable(optimization_container, :P__HydroEnergyReservoir)
    var_sp = PSI.get_variable(optimization_container, :Sp__HydroEnergyReservoir)
    var_e = PSI.get_variable(optimization_container, :E__HydroEnergyReservoir)
    energy_balance = PSI.get_constraint(optimization_container, :energy_capacity__HydroEnergyReservoir)
    constraint_budget = PSI.add_cons_container!(
        optimization_container,
        :Hydro_Budget_HydroEnergyReservoir,
        collect(keys(inflow_data)),
    )
    constraint_balance = PSI.add_cons_container!(
        optimization_container,
        :PH_EnergyBalance_HydroEnergyReservoir,
        collect(keys(inflow_data)),
        1:no_scenarios,
    )
#     constraint_scenarios = PSI.add_cons_container!(
#         optimization_container,
#         :PH_scenario_HydroEnergyReservoir,
#         collect(keys(budget_data)),
#     )
    budget_var = PSI.add_var_container!(
        optimization_container,
        :Pb_HydroEnergyReservoir,
        collect(keys(inflow_data)),
        1:no_scenarios,
    )
    energy_var = PSI.add_var_container!(
        optimization_container,
        :Eb_HydroEnergyReservoir,
        collect(keys(inflow_data)),
        1:no_scenarios,
    )
    jump_model = PSI.get_jump_model(optimization_container)
    for (name, inflow) in inflow_data
        for sec in 1:no_scenarios
            budget_var[name, sec] = JuMP.@variable(
                jump_model,
                base_name = "Pb_HydroEnergyReservoir_{$(name)_{$(sec)}}",
            )
            energy_var[name, sec] = JuMP.@variable(
                jump_model,
                base_name = "Eb_HydroEnergyReservoir_{$(name)_{$(sec)}}",
            )
        end
        for sec in 1:no_scenarios
            if sec == 1
#                constraint_balance[name, sec] = JuMP.@constraint(
#                     jump_model, 
#                     energy_var[name, sec] == inflow_data[name][sec] + initial_condtions[name] - budget_var[name, sec] 
#                 )
#             else                
#                constraint_balance[name, sec] = JuMP.@constraint(
#                     jump_model, 
#                     energy_var[name, sec] == energy_var[name, sec-1] + inflow_data[name][sec] - budget_var[name, sec] 
#                 )
            end
        end
        if scenario_no != 0
            JuMP.delete(jump_model, energy_balance[name, 1])
            energy_balance[name, 1] = JuMP.@constraint(jump_model, 
                var_e[name, 1] == energy_var[name, scenario_no] + inflow[scenario_no+1] -  var_sp[name, 1] - var_h[name, 1]
            )
        end
        JuMP.@constraint(jump_model, var_e[name, 24] == energy_var[name, scenario_no+1])
        constraint_budget[name] =
            JuMP.@constraint(
            jump_model, sum(var_h[name, t] for t in 1:24) == budget_var[name, scenario_no+1]
        )
#         constraint_scenarios[name] = 
#             JuMP.@constraint(jump_model, sum(budget_var[name, sec] for sec in 1:no_scenarios ) == budget)
    end
    return
end   


function get_subproblem(problem::PSI.OperationsProblem)
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    
    first_stge_var =  JuMP.VariableRef[]
    secound_stge_var = JuMP.VariableRef[]
    for (name, cont) in PSI.get_variables(problem) 
        if (name == :Pb_HydroEnergyReservoir) || (name == :Eb_HydroEnergyReservoir)
            push!(first_stge_var, cont...)
        else
            push!(secound_stge_var, cont...)
        end
    end
    variable_map = Dict{PH.StageID, Vector{JuMP.VariableRef}}(
        [
            PH.stid(1) => first_stge_var,
            PH.stid(2) => secound_stge_var,
        ]
    )
    return jump_model, variable_map
end

function add_subproblem_variables!(problem, scenario_id, budget_vars)
    id = PH.value(scenario_id) + 1
    optimization_container = PSI.get_optimization_container(problem)
    container = PSI.get_variable(optimization_container, :Pb__HydroEnergyReservoir)
    for d in budget_vars.axes[1]
        container[d, id] = budget_vars[d, id]
    end
    return 
end

function add_subproblem_variables!(problem, scenario_id, budget_vars, energy_vars)
    id = PH.value(scenario_id) + 1
    optimization_container = PSI.get_optimization_container(problem)
    container = PSI.get_variable(optimization_container, :Pb__HydroEnergyReservoir)
    container_e = PSI.get_variable(optimization_container, :Eb__HydroEnergyReservoir)
    for d in budget_vars.axes[1]
        container[d, id] = budget_vars[d, id]
        container_e[d, id] = energy_vars[d, id]
    end
    return 
end

function create_model(
    scenario_id::PH.ScenarioID,
    problem;
    kwargs...
)
#     budget_data = problem.ext["hydro_budget"]
    inflow_data = problem.ext["inflow_data"]
    initial_condtions = problem.ext["initial_condtions"]
    total_scenarios = problem.ext["total_scenarios"]
    sim_initial_time = problem.internal.optimization_container.settings.initial_time.x
    initial_time = collect(sim_initial_time:Day(1):sim_initial_time+Day(total_scenarios))
    scenario_init_time = Dict(ix => init for (ix, init) in enumerate(initial_time));
    system = problem.ext["ph_system"] 
    template = get(problem.ext, "ph_template", nothing)
    sub_problem = build_operations_problem(
        system,
        scenario_init_time[PH.value(scenario_id) + 1],
        JuMP.optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 1e-5),
        PSI.CopperPlatePowerModel,
        "./simulation_folder/",
        template,
    )

#     populate_ext!(sub_problem, PH.value(scenario_id), budget_data, total_scenarios)
    populate_ext!(sub_problem, PH.value(scenario_id), inflow_data, initial_condtions, total_scenarios)
    populate_PH_model!(sub_problem)
    jump_model, variable_map = get_subproblem(sub_problem)
    budget_vars = get_variables(sub_problem)[:Pb_HydroEnergyReservoir]
    energy_vars = get_variables(sub_problem)[:Eb_HydroEnergyReservoir]
#     add_subproblem_variables!(problem, scenario_id, budget_vars)
    add_subproblem_variables!(problem, scenario_id, budget_vars, energy_vars)
    return JuMPSubproblem(jump_model, scenario_id, variable_map)
end


### Custom Stage


add_to_ext!(p::PSI.OperationsProblem, key, data) = p.ext[key] = data

function PSI.OperationsProblem(
    ::Type{PHProblem},
    template::OperationsProblemTemplate,
    sys::PSY.System,
    budget_data,
    total_scenarios,
    uc_sys, 
    jump_model::Union{Nothing, JuMP.AbstractModel} = nothing;
    kwargs...,
)
    problem = OperationsProblem{PHProblem}(template, sys, jump_model; kwargs...)
    add_to_ext!(problem, "hydro_budget", budget_data)
    add_to_ext!(problem, "total_scenarios", total_scenarios)
    add_to_ext!(problem, "ph_system", uc_sys)
    return problem
end

function PSI.OperationsProblem(
    ::Type{PHReservoirProblem},
    template::OperationsProblemTemplate,
    sys::PSY.System,
    inflow_data,
    initial_condtions,
    total_scenarios,
    ph_system,
    ph_template,
    jump_model::Union{Nothing, JuMP.AbstractModel} = nothing;
    kwargs...,
)
    problem = OperationsProblem{PHReservoirProblem}(template, sys, jump_model; kwargs...)
    add_to_ext!(problem, "inflow_data", inflow_data)
    add_to_ext!(problem, "total_scenarios", total_scenarios)
    add_to_ext!(problem, "ph_system", ph_system)
    add_to_ext!(problem, "ph_template", ph_template)
    add_to_ext!(problem, "initial_condtions", initial_condtions)
    return problem
end

function PSI._build!(problem::PSI.OperationsProblem{PHProblem}, serialize::Bool)
    ## TODO: Add Container for Storaging Budget
    PSI.build_pre_step!(problem)
    system = get_system(problem)
    devices = get_available_components(HydroEnergyReservoir, system)
    ## Need to fix ph_problem.internal.optimization_container.time_steps
    total_scenarios = problem.ext["total_scenarios"]
    container = PSI.container_spec(
        JuMP.VariableRef,
        [PSY.get_name(d) for d in devices],
        1:total_scenarios,
    )
    PSI.assign_variable!(problem.internal.optimization_container, "Pb", PSY.HydroEnergyReservoir, container)
    PSI.set_status!(problem, BuildStatus.BUILT)
    return BuildStatus.BUILT
end

function PSI._build!(problem::PSI.OperationsProblem{PHReservoirProblem}, serialize::Bool)
    ## TODO: Add Container for Storaging Budget
    PSI.build_pre_step!(problem)
    system = get_system(problem)
    devices = get_available_components(HydroEnergyReservoir, system)
    ## Need to fix ph_problem.internal.optimization_container.time_steps
    total_scenarios = problem.ext["total_scenarios"]
    container = PSI.container_spec(
        JuMP.VariableRef,
        [PSY.get_name(d) for d in devices],
        1:total_scenarios,
    )
    energy_container = PSI.container_spec(
        JuMP.VariableRef,
        [PSY.get_name(d) for d in devices],
        1:total_scenarios,
    )
    PSI.assign_variable!(problem.internal.optimization_container, "Pb", PSY.HydroEnergyReservoir, container)
    PSI.assign_variable!(problem.internal.optimization_container, "Eb", PSY.HydroEnergyReservoir, energy_container)

    PSI.set_status!(problem, BuildStatus.BUILT)
    return BuildStatus.BUILT
end

function PSI.solve_impl(problem::PSI.OperationsProblem{PHProblem}; optimizer = nothing)

    total_scenarios = problem.ext["total_scenarios"]
    st = PH.two_stage_tree(total_scenarios)
    (n, err, rerr, obj, soln, phd) = PH.solve(st,
                                              create_model,
                                              PH.ScalarPenaltyParameter(45.0), problem,
                                              atol=1e-8, rtol=1e-12, max_iter=500,
                                              report=100, # print residual info every 10 iterations
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
    ## 
    PSI.set_run_status!(problem, RunStatus.SUCCESSFUL)
    return RunStatus.SUCCESSFUL
end

function PSI.solve_impl(problem::PSI.OperationsProblem{PHReservoirProblem}; optimizer = nothing)

    total_scenarios = problem.ext["total_scenarios"]
    st = PH.two_stage_tree(total_scenarios)
    (n, err, rerr, obj, soln, phd) = PH.solve(st,
                                              create_model,
                                              PH.ScalarPenaltyParameter(45.0), problem,
                                              atol=1e-8, rtol=1e-12, max_iter=200,
                                              report=100, # print residual info every 10 iterations
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
    ## 
    PSI.set_run_status!(problem, RunStatus.SUCCESSFUL)
    return RunStatus.SUCCESSFUL
end

function PSI.write_problem_results!(
    step::Int,
    problem::PSI.OperationsProblem{PHProblem},
    start_time::Dates.DateTime,
    store::PSI.SimulationStore,
    exports,
)
#     stats = OptimizerStats(problem, step)
#     write_optimizer_stats!(store, get_name(problem), stats, start_time)
    PSI.write_model_results!(store, problem, start_time; exports = exports)
    return
end

function PSI.write_problem_results!(
    step::Int,
    problem::PSI.OperationsProblem{PHReservoirProblem},
    start_time::Dates.DateTime,
    store::PSI.SimulationStore,
    exports,
)
#     stats = OptimizerStats(problem, step)
#     write_optimizer_stats!(store, get_name(problem), stats, start_time)
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