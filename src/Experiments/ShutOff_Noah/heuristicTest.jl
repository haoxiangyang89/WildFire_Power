## This test is design to compare the performance of solutions computed by Noah's paper and the algorithm

include("testData_RTS_New/gurobiTest.jl")



heuristicSolutionList = load("src/Experiments/ShutOff_Noah/heuristicSolutionList.jld2")["heuristicSolutionList"]
NoahSolutionList = load("src/Experiments/ShutOff_Noah/NoahSolutionList.jld2")["NoahSolutionList"]
riskValueBasedPlan = load("src/Experiments/ShutOff_Noah/riskValueBasedPlan.jld2")["riskValueBasedPlan"]
gurobiResultList(14, 0.1)   

## ------------------------------  Compute total costs by using different decisions  ----------------------------- ##
## scenario size = 200
totalCost = Dict{Float64, Float64}()
# costDistribution = Dict{Tuple{Float64, Int64}, NamedTuple{(:damageCost, :loadShedding), Tuple{Float64, Float64}}}()

indexSets = load("testData_RTS_New/indexSets.jld2")["indexSets"]
paramOPF = load("testData_RTS_New/paramOPF.jld2")["paramOPF"]
paramDemand = load("testData_RTS_New/paramDemand.jld2")["paramDemand"]
Ω_rv = load("testData_RTS_New/wholeSpace.jld2")["wholeSpace"]
indexSets.Ω = [1:101...]
prob = load("testData_RTS_New/prob.jld2")["prob"]

for alpha in keys(riskValueBasedPlan) 

    state_variable = riskValueBasedPlan[alpha].state_variable;
    state_variable = gurobiResultList[.1, 14].first_state_variable
    # state_variable = NoahSolutionList[alpha];

    ## -------------------------------- solve the first stage model -------------------------------- ##
    (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, indexSets.Ω);
    (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L);

    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                        "OutputFlag" => 1, 
                                        "Threads" =>0, 
                                        "MIPGap" => 5e-3, 
                                        "TimeLimit" => 600) 
                                        );
                                        
    @variable(Q, θ_angle[B, 1:T]);      ## phase angle of the bus i
    @variable(Q, P[L, 1:T]);       ## real power flow on line l; elements in L is Tuple (i, j)
    @variable(Q, s[G, 1:T]);       ## real power generation at generator g
    @variable(Q, 0 <= x[D, 1:T] <= 1);  ## load shedding
    # @variable(Q, slack_variable[B, 1:T])

    # constraint 3b 3c
    for l in L
        i = l[1];
        j = l[2];
        @constraint(Q, [t in 1:T], P[l, t] <= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmax * (1 - state_variable[:zl][l, t] ) ) );
        @constraint(Q, [t in 1:T], P[l, t] >= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmin * (1 - state_variable[:zl][l, t] ) ) );
    end

    ## constraint 1d
    @constraint(Q, [l in L, t in 1:T], P[l, t] >= - paramOPF.W[l] * state_variable[:zl][l, t] );
    @constraint(Q, [l in L, t in 1:T], P[l, t] <= paramOPF.W[l] * state_variable[:zl][l, t] );

    ## constraint 1e
    @constraint(Q, [i in B, t in 1:T], sum(s[g, t] for g in Gᵢ[i]) + 
                                        sum(P[(i, j), t] for j in out_L[i]) - 
                                        sum(P[(j, i), t] for j in in_L[i]) 
                                        .>=  sum(paramDemand.demand[t][d] * x[d, t] for d in Dᵢ[i]) );
    
    ## constraint 1f
    @constraint(Q, [g in G, t in 1:T], s[g, t] >= paramOPF.smin[g] * state_variable[:zg][g, t] );
    @constraint(Q, [g in G, t in 1:T], s[g, t] <= paramOPF.smax[g] * state_variable[:zg][g, t] );

    # constraint 1g h i j
    for i in B 
        @constraint(Q, [t in 1:T, d in Dᵢ[i]], state_variable[:zb][i, t] >= x[d, t] );
        # @constraint(Q, [t in 1:T, g in Gᵢ[i]], state_variable[:zb][i, t] >= state_variable[:zg][g, t]);
        # @constraint(Q, [t in 1:T, j in out_L[i]], state_variable[:zb][i, t] >= state_variable[:zl][(i, j), t] );
        # @constraint(Q, [t in 1:T, j in in_L[i]], state_variable[:zb][i, t] >= state_variable[:zl][(j, i), t] );
    end


    ## objective function
    @objective(Q, Min, sum( prob[ω] * 
                                    ( sum( sum(paramDemand.w[d] * (1 - x[d, t]) for d in D ) for t in 1:Ω_rv[ω].τ - 1 ) )
                                                                                                                                    for ω in Ω )  
            );
    optimize!(Q);
    state_value  = JuMP.objective_value(Q); 
                                                                                                                                    
    ## stage 2
    c = 0.0;
    for ω in indexSets.Ω
        @info "$alpha, $ω"
        randomVariables = Ω_rv[ω];
        forward2Info = forward_stage2_model!(indexSets, 
                                                paramDemand,
                                                paramOPF,
                                                randomVariables              
                                                );

        ẑ = Dict(   :zg => state_variable[:zg][:, randomVariables.τ - 1], 
                    :zb => state_variable[:zb][:, randomVariables.τ - 1], 
                    :zl => state_variable[:zl][:, randomVariables.τ - 1]
                    );
        ## modify the constraints according to the first stage state variables
        forward_stage2_modify_constraints!(indexSets, 
                                            forward2Info,
                                            1,
                                            ẑ,
                                            randomVariables
                                            );

        ####################################################### solve the model and display the result ###########################################################
        optimize!(forward2Info.model);
        # damageCost = sum(paramDemand.cb[i] * value(forward2Info.νb[i]) for i in indexSets.B) + 
        #                         sum(paramDemand.cg[g] * value(forward2Info.νg[g]) for g in indexSets.G) + 
        #                                         sum(paramDemand.cl[l] * value(forward2Info.νl[l]) for l in indexSets.L);
        # loadShedding = sum( sum(paramDemand.w[d] * (1 - value(x[d, t])) for d in D ) for t in 1:Ω_rv[ω].τ - 1 ) + (JuMP.objective_value(forward2Info.model) - damageCost);
        c = c + prob[ω] * JuMP.objective_value(forward2Info.model);
        # costDistribution[alpha, ω] = (damageCost = damageCost, loadShedding = loadShedding)
    end

    # sum(costDistribution[alpha, ω].damageCost for ω in 1:100)
    # sum(costDistribution[alpha, ω].loadShedding for ω in 1:100)


    u = state_value + c;
    totalCost[alpha] = u;
    save("src/Experiments/ShutOff_Noah/totalCost.jld2", "totalCost", totalCost)
end

save("src/Experiments/ShutOff_Noah/totalCost.jld2", "totalCost", totalCost)

totalCost = load("src/Experiments/ShutOff_Noah/totalCost.jld2")["totalCost"]

save("src/Experiments/ShutOff_Noah/costDistribution.jld2", "costDistribution", costDistribution)

