
#+ 

using Debugger
using Pkg
# Activate environment that has ProgressiveHedging installed
Pkg.activate("..")


#+ 

Pkg.instantiate()


#+ 

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


#+ 

using JuMP


#+ 

@everywhere function create_model(
        scenario_id::PH.ScenarioID,
        key_word_arg::String="Default key word argument", kwargs...
        )
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
    d = [1500, 1600];
    d1 = [2200, 2350];
    # Wind forecast
    w_f = [200.0, 150.0];

    scen1 = [1.0, 0.0]
    scen2 = [0.0, 1.0];
    
    ed = Model(Xpress.Optimizer)
    
    # Define decision variables    
    @variable(ed, 0 <= g[i = 1:2, t = 1:2] <= g_max[i]) # power output of generators
    @variable(ed, 0 <= h[i = 1:2, t = 1:2] <= h_max[i]) # power output of generators
    @variable(ed, 0 <= w[t = 1:2] <= w_f[t]) # wind power injection
    @variable(ed, 0 <= h_b[i = 1:2] <= h_budget) # power output of generators

    # Define the objective function
    @objective(ed, Min, sum(c_g[i]*(g[i, t] ) for i in 1:2, t in 1:2) + sum(c_w *w[t] for t in 1:2))

    # Define the constraint on the maximum and minimum power output of each generator
    @constraint(ed, [i = 1:2, t = 1:2], g[i, t] <= g_max[i]) #maximum
    @constraint(ed, [i = 1:2, t = 1:2], g[i, t] >= g_min[i]) #minimum

    @constraint(ed, [i = 1:2, t = 1:2], h[i, t] <= h_max[i]) #maximum
    @constraint(ed, [i = 1:2, t = 1:2], h[i, t] >= h_min[i]) #minimum
    
    # Define the constraint on the wind power injection
    if scenario_id == PH.scid(0)
        @constraint(ed, sum(h[i, t] for i in 1:2, t in 1:2) <= sum(h_b[i]*scen1[i] for i in 1:2))
    else
        @constraint(ed, sum(h[i, t] for i in 1:2, t in 1:2) <= sum(h_b[i]*scen2[i] for i in 1:2))
    end
    @constraint(ed, sum(h_b[i] for i in 1:2) == h_budget)
    
    # Define the constraint on the wind power injection
    @constraint(ed, [t = 1:2], w[t] <= w_f[t])

    # Define the power balance constraint
    if scenario_id == PH.scid(0)
        @constraint(ed, [t = 1:2], sum(g[:, t] + h[:, t]) + w[t] == d[t])
    else
        @constraint(ed, [t = 1:2], sum(g[:, t] + h[:, t]) + w[t] == d1[t])
    end
    
    vdict = Dict{PH.StageID, Vector{JuMP.VariableRef}}([PH.stid(1) => VariableRef[h_b...],
                                                        PH.stid(2) => vcat(VariableRef[g...], VariableRef[h...], VariableRef[w...]),
                                                        ])
    # vdict = Dict{PH.StageID, Vector{JuMP.VariableRef}}([PH.stid(1) => vcat(VariableRef[h_b...], VariableRef[g...]),
    #                                                     PH.stid(2) => vcat(VariableRef[h...], VariableRef[w...]),
    #                                                     ])
    return JuMPSubproblem(ed, scenario_id, vdict)
end

# Solve the economic dispatch problem


#+ 

function build_scen_tree()

    probs = [0.5, 0.5]
    
    tree = PH.ScenarioTree()
    
    for l in 1:2
        PH.add_leaf(tree, tree.root, probs[l])
    end
    return tree
end
;


#+ 

ef_model = @time PH.solve_extensive(build_scen_tree(),
    create_model, 
    ()->Ipopt.Optimizer(),
    "Unused example string",
    opt_args=(print_level=0,)
)
println(ef_model)


#+ 

println(JuMP.termination_status(ef_model))
println(JuMP.objective_value(ef_model))
for var in JuMP.all_variables(ef_model)
    println("$var = $(JuMP.value(var))")
end


#+ 

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


#+ 

aobj = PH.retrieve_aug_obj_value(phd)
println("Augmented Objective: ", aobj)
println("Difference: ", aobj - obj)


#+ 

soln


#+ 


