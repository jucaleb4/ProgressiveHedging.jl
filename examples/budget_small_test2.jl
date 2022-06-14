#=
NREL
@Author: Caleb Ju (cju@nrel.gov)

A preliminary implementation of using ProgressiveHedging (PH) to decompose a
multistage (or multi-day) problem into independent sub-problems.  We will
test on a fictionary economic dispatch problem where we have access to a
battery, and each day has different scenarios. We let T=4 be the number of
days we model.
=#

using Pkg
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

@everywhere using ProgressiveHedging
@everywhere const PH = ProgressiveHedging
@everywhere using Ipopt
@everywhere using JuMP
@everywhere using Xpress

using JuMP


#=
See Section 3.1 of the write-up. Here, we just do sample and time decomposition
all in one level.

TODO: Add time decomposition
=#
@everywhere function create_model(
        scenario_id::PH.ScenarioID,
        key_word_arg::String="Default key word argument", kwargs...
        )
    ed = Model(Xpress.Optimizer)
    T = 4 # Number of stages (we will use t to represent time)
    s = 4 # number of scenarios
    n = 2 # Number of generators

    #Define the economic dispatch (ED) model
    g_max = [1000, 1000];
    # Minimum power output of generators
    g_min = [0, 300];

    h_max = [100, 100];
    h_min = [0, 0];

    h_budget = 500 ;
    # Incremental cost of generators 
    c_g = [50, 100];
    # Fixed cost of generators
    c_g0 = [1000, 0]
    # Incremental cost of wind generators
    c_w = 50;
    # Total demand
    d = [1500 1600 1400 1500; 2200 2350 2100 2000; 1700 2100 2100 2150; 2500 1800 2000 2500];

    # Wind forecast
    w_f = [200.0, 150.0, 75.0, 150.0];

    # Define decision variables    
    stage1 = Vector{JuMP.VariableRef}()
    vref = @variable(ed, 0 <= b[t = 1:T+1] <= h_budget) # total budget remaining
    append!(stage1, vref)

    stage2 = Vector{JuMP.VariableRef}()
    # power generators (cost)
    vref_1 = @variable(ed, g[i = 1:n, t = 1:T]) 
    # power generators (free)
    vref_2 = @variable(ed, h[i = 1:n, t = 1:T]) 
    # wind power
    vref_3 = @variable(ed, w[i = 1:n, t = 1:T]) 
    # budget
    vref_4 = @variable(ed, h_b[i = 1:n, t = 1:T]) 
    append!(stage2, vref_1)
    append!(stage2, vref_2)
    append!(stage2, vref_3)
    append!(stage2, vref_4)


    # Define the constraint on the maximum and minimum power output of each generator
    @constraint(ed, [i = 1:n, t = 1:T], g_min[i] <= g[i, t] <= g_max[i]) 

    @constraint(ed, [i = 1:n, t = 1:T], h_min[i] <= h[i, t] <= h_max[i])

    # Define the constraint on the wind power injection
    @constraint(ed, [t = 1:n], w[t] <= w_f[t])

    # Budget constraint
    @constraint(ed, [t = 1:T], sum(h[i, t] for i in 1:n) <= h_b[t])

    # Budget constraint tying different time stages
    @constraint(ed, b[1] == h_budget)
    @constraint(ed, [t=1:T], b[t+1] - b[t] + h_b[t] == 0)
    
    # Define the constraint on the wind power injection and power balance 
    # if scenario_id == PH.scid(0)
    _id = PH.value(scenario_id)
    if 0 <= _id <= 15
        t = Int(floor(_id/T))+1
        _sc_id = (_id % T) + 1

        # Define the objective function
        @objective(ed, Min, sum(c_g[i]*g[i, t] for i in 1:n) + c_w *w[t])

        @constraint(ed, [t = 1:T], sum(g[:, t] + h[:, t]) + w[t] == d[t,_sc_id])
    else
        # TODO: print the wrong scenario ID
        println("Scenario id:", _id)
        @error "Incorrect scenario id"
    end

    vdict = Dict{PH.StageID, Vector{JuMP.VariableRef}}([PH.stid(1) => stage1,
                                                        PH.stid(2) => stage2,
                                                        ])

    return JuMPSubproblem(ed, scenario_id, vdict)
end

function build_scen_tree()
    s = 16
    probs = [1.0/s for i=1:s]
    tree = PH.ScenarioTree()
    
    for l in 1:s
        PH.add_leaf(tree, tree.root, probs[l])
    end
    return tree
end
;

ef_model = @time PH.solve_extensive(build_scen_tree(),
    create_model, 
    ()->Ipopt.Optimizer(),
    "Unused example string",
    opt_args=(print_level=0,)
)
println(ef_model)

println(JuMP.termination_status(ef_model))
println(JuMP.objective_value(ef_model))
for var in JuMP.all_variables(ef_model)
    println("$var = $(JuMP.value(var))")
end


st = build_scen_tree()
(n, err, rerr, obj, soln, phd) = PH.solve(st,
                                          create_model,
                                          PH.ScalarPenaltyParameter(25.0), "Passed to constructor.",
                                          atol=1e-8, rtol=1e-12, max_iter=1000,
                                          report=10, # print residual info every 10 iterations
                                          )
println("Number of iterations: ", n)
println("L^2 error: ", err)
println(obj)

aobj = PH.retrieve_aug_obj_value(phd)
println("Augmented Objective: ", aobj)
println("Difference: ", aobj - obj)
@show soln
