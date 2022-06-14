#using Pkg
#Pkg.activate("..")

using Distributed
const WORKERS = 7 # Change to > 1 to use parallel
if nworkers() < WORKERS
    diff = (nprocs() == nworkers() ? WORKERS : WORKERS - nworkers())
    println("Adding $diff worker processes.")
    Distributed.addprocs(diff, exeflags="--project")
    # Make sure these workers also have an environment with PH installed
    #@everywhere using Pkg
    #for w in workers()
        #@spawnat(w, Pkg.activate(".."))
    #end
end

# Distributed setup
#@everywhere using Pkg
#@everywhere Pkg.activate("..")
# @everywhere ENV["XPRESSDIR"]="/nopt/nrel/apps/xpressmp/8.13.0"
@everywhere using Xpress
@everywhere using ProgressiveHedging
@everywhere const PH = ProgressiveHedging
@everywhere using Ipopt
@everywhere using JuMP, PowerSimulations, PowerSystems, PowerSystemCaseBuilder, TimeSeries
@everywhere const PSI = PowerSimulations
@everywhere const PSY = PowerSystems
@everywhere const PSB = PowerSystemCaseBuilder

# Local setup
using JuMP, Xpress, PowerSimulations, PowerSystems, PowerSystemCaseBuilder, InfrastructureSystems, Dates, TimeSeries
const PSI = PowerSimulations
const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder
solver = optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 1e-4)

# Setup the optimization problem
@everywhere function build_operations_problem(
        system::PSY.System,
        initial_time,
        solver,
        network = CopperPlatePowerModel,
        directory ="./simulation_folder/",
    )
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
    set_device_model!(template_ed, PowerSystems.GenericBattery, BookKeeping)
    set_device_model!(template_ed, PowerSystems.RenewableFix, FixedOutput)

    # Ancillary services
    # set_service_model!(template_ed, VariableReserve{ReserveUp}, RangeReserve)
    # set_service_model!(template_ed, VariableReserve{ReserveDown}, RangeReserve)

    problem = DecisionModel(
        template_ed, 
        system, 
        name = "SubProblem",
        optimizer = solver, 
        warm_start = false,
        export_pwl_vars = true,
        initialize_model = false,
        initial_time = initial_time,
    )
    build!(problem, output_dir = directory)
    return problem
end
   
@everywhere function populate_ext!(problem::DecisionModel{GenericOpProblem}, scenario_no, params)
    problem.ext["scenario_no"] = scenario_no
    problem.ext["efficiency_data"] = params["efficiency_data"]
    problem.ext["total_scenarios"] = params["total_scenarios"] 
    return
end

@everywhere struct EnergyPH <: PSI.VariableType end
@everywhere struct StorgeConsensusPH <: PSI.ConstraintType end

# function populate_PH_model!(problem::PSI.DecisionProblem)
@everywhere function populate_PH_model!(problem::DecisionModel{GenericOpProblem})
    no_scenarios = problem.ext["total_scenarios"]
    time_no = problem.ext["scenario_no"]
    time_no += 1
    efficiency_data = problem.ext["efficiency_data"]
    
    optimization_container = PSI.get_optimization_container(problem)
    # https://github.com/NREL-SIIP/PowerSimulations.jl/search?q=time_steps
    # (Depreciated?) time_steps = PSI.model_time_steps(optimization_container)
    time_steps = PSI.get_time_steps(optimization_container)
    var_pin = PSI.get_variable(optimization_container, ActivePowerInVariable(), GenericBattery)
    var_pout= PSI.get_variable(optimization_container, ActivePowerOutVariable(), GenericBattery)
    var_e   = PSI.get_variable(optimization_container, EnergyVariable(), GenericBattery)
    e_balance = PSI.get_constraint(optimization_container, EnergyBalanceConstraint(), GenericBattery)
    
    # the first element is the name of the battery, the second is the number of batteries in...
    device_names, timesteps = axes(var_e)
    num_hours = size(timesteps)[1] # this value is determined by @transform_single_time_series!
    
    storage_consensus = PSI.add_constraints_container!(
        optimization_container,
        StorgeConsensusPH(),
        GenericBattery,
        device_names,
    )

    var_ph = PSI.add_variable_container!(
        optimization_container,
        EnergyPH(), 
        GenericBattery,
        device_names,
        1:no_scenarios,
    )
    jump_model = PSI.get_jump_model(optimization_container)
    for name in device_names
        for sec in 1:no_scenarios
            var_ph[name, sec] = JuMP.@variable(
                optimization_container.JuMPmodel,
                base_name = "Energy_{$(name)_{$(sec)}}",
            )
        end
        
        if time_no > 1
            (eff_in, eff_out) = efficiency_data[name]
            JuMP.delete(jump_model, e_balance[name, 1])
            
            e_balance[name, 1] = JuMP.@constraint(jump_model, 
                var_e[name, 1] == var_ph[name, time_no-1] + eff_in * var_pin[name,1] 
                                  - var_pout[name,1]/eff_out
            )
        end
        
        storage_consensus[name] = 
            JuMP.@constraint(jump_model, var_e[name, num_hours] == var_ph[name, time_no])
        
    end
    return
end

# function get_subproblem(problem::PSI.DecisionProblem)
@everywhere function get_subproblem(problem::DecisionModel{GenericOpProblem})
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    
    first_stge_var =  JuMP.VariableRef[]
    secound_stge_var = JuMP.VariableRef[]
    for (var_key, cont) in PSI.get_variables(problem) 
        if var_key == PSI.VariableKey(EnergyPH, GenericBattery)
            push!(first_stge_var, cont...)
        else
            push!(secound_stge_var, cont...)
        end
    end
    variable_map = Dict{PH.StageID, Vector{JuMP.VariableRef}}(
        [
            PH.stid(1) => first_stge_var,
            PH.stid(2) => secound_stge_var,
        ])
    
    return jump_model, variable_map
end

@everywhere function create_model(
    scenario_id::PH.ScenarioID,
    params;
    kwargs...
)
    scen_len = params["scenario_len"] # think of as number of days in each sub-problem
    solver = params["solver"]
    sys_name = params["sys_name"]
    hour_len = 24 * scen_len
    initial_time = DateTime("2020-01-01T00:00:00") + Dates.Day(scen_len*PH.value(scenario_id))
    system = PSY.System(sys_name, time_series_directory = "/tmp/scratch")
    # system = build_system(PSITestSystems, "modified_RTS_GMLC_DA_sys", time_series_directory = "/tmp/scratch")
    # @enter PSY.transform_single_time_series!(system, 24, Hour(24))
    PSY.transform_single_time_series!(system, hour_len, Hour(hour_len))
    problem = build_operations_problem(
        system,
        initial_time,
        solver,
        CopperPlatePowerModel,
        "./simulation_folder/"
    )
    
    populate_ext!(problem, PH.value(scenario_id), params)
    populate_PH_model!(problem)
    jump_model, variable_map = get_subproblem(problem)
    return JuMPSubproblem(jump_model, scenario_id, variable_map)
end

function build_scen_tree(D)
    probs = [1.0/d for d=1:D]
    tree = PH.ScenarioTree()
    for l in 1:D
        PH.add_leaf(tree, tree.root, probs[l])
    end
    return tree
end

import Logging
Logging.disable_logging(Logging.Warn)

# Borrowed from previous notebook from Sourabh
sys_name = "/lustre/eaglefs/projects/pvb/cju/gen_systems/data/rts_with_battery_060822_sys.json"
system = PSY.System(sys_name, time_series_directory = "/tmp/scratch")
# system = build_system(PSITestSystems, "modified_RTS_GMLC_DA_sys", time_series_directory = "/tmp/scratch")
efficiency_data = Dict()
for h in  get_components(GenericBattery, system)
    name = PSY.get_name(h)
    efficiency_data[name] = PSY.get_efficiency(h)
end

params = Dict()
params["efficiency_data"] = efficiency_data
params["solver"] = solver
params["sys_name"] = sys_name

if(size(ARGS)[1] >= 2)
    num_scen = parse(Int64, ARGS[1])
    params["scenario_len"] = parse(Int64, ARGS[2])
else
    num_scen = 2
    params["scenario_len"] = 7
end
scen_len = params["scenario_len"]
params["total_scenarios"] = num_scen
println("$(num_scen*scen_len) broken into $(num_scen) groups of $(scen_len) days")

st = PH.two_stage_tree(num_scen)

#=
ef_model = PH.solve_extensive(st,
    create_model, 
    ()->Xpress.Optimizer(),
    params,
)
objval = JuMP.objective_value(ef_model)
println("[Ext] Partial obj val: $(objval) || full val: $(objval*num_scen)")

function lastline(file)
    last = ""
    open(file) do io
        while !eof(io)
            last = readline(io)
        end
    end
    return last
end
=#

max_iter = 1000
# max_iter = parse(Int64, lastline("holder.txt"))
println("Max iters: ", max_iter)

(n, err, rerr, obj, soln, phd) = PH.solve(st,
                                          create_model,
                                          PH.ScalarPenaltyParameter(25.0), 
                                          params,
                                          atol=1e-8, rtol=1e-12, max_iter=max_iter,
                                          report=1, # print residual info every 10 iterations
                                          # lower_bound=10, # lower bound from Gade
                                          )
println("[PH] Number of iterations: ", n)
println("[PH] L^2 error: ", err)
println("[PH] obj: ", obj)
