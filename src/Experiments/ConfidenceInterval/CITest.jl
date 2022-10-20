using CSV, DataFrames, Printf, Gurobi, JuMP
using JLD2, FileIO

# totalCost = load("src/Experiments/ScenariosSizeTest/totalCost.jld2")["totalCost"]  # find the solution with smallest cost 
gurobiResultList = load("src/Experiments/ScenariosSizeTest/gurobiResultList.jld2")["gurobiResultList"]
# (13, 200)

zstar = gurobiResultList[200,13]


## ------------------------------  Compute total costs by using different decisions  ----------------------------- ##
## scenario size = 200
totalCost = Dict{Tuple{Int64, Int64}, Float64}()

indexSets = load("testData_RTS/indexSets.jld2")["indexSets"]
paramOPF = load("testData_RTS/paramOPF.jld2")["paramOPF"]
paramDemand = load("testData_RTS/paramDemand.jld2")["paramDemand"]
Ω_rvList = load("testData_RTS/Ω_rvList.jld2")["Ω_rvList"]

for k in 1:20 
  for scenarioSize in [20, 50, 100, 200, 500]
      Ω_rv = Ω_rvList[scenarioSize, k]
      indexSets.Ω = [1:scenarioSize...]
      prob = Dict{Int64, Float64}();
      for ω in 1:scenarioSize 
        prob[ω] = 1/scenarioSize;
      end

      state_variable = zstar.first_state_variable;
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
                                            .>= sum(paramDemand.demand[t][d] * x[d, t] for d in Dᵢ[i]));
      
      ## constraint 1f
      @constraint(Q, [g in G, t in 1:T], s[g, t] >= paramOPF.smin[g] * state_variable[:zg][g, t] );
      @constraint(Q, [g in G, t in 1:T], s[g, t] <= paramOPF.smax[g] * state_variable[:zg][g, t] );

      ## constraint 1g h i j
      for i in B 
        @constraint(Q, [t in 1:T, d in Dᵢ[i]], state_variable[:zb][i, t] >= x[d, t] );
        @constraint(Q, [t in 1:T, g in Gᵢ[i]], state_variable[:zb][i, t] >= state_variable[:zg][g, t]);
        @constraint(Q, [t in 1:T, j in out_L[i]], state_variable[:zb][i, t] >= state_variable[:zl][(i, j), t] );
        @constraint(Q, [t in 1:T, j in in_L[i]], state_variable[:zb][i, t] >= state_variable[:zl][(j, i), t] );
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
        @info "$k, $ω, $scenarioSize"
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
        c = c + prob[ω] * JuMP.objective_value(forward2Info.model);
      end
      u = state_value + c;
      totalCost[(k, scenarioSize)] = u;
  end
end

## cost of 5000-dataset when using the best first-stage solution generated in ScenariosSizeTest
save("src/Experiments/ConfidenceInterval/totalCost.jld2", "totalCost", totalCost)
totalCost = load("src/Experiments/ConfidenceInterval/totalCost.jld2")["totalCost"]




## ------------------------- 10,000 sample size ----------------------------- ##
## scenario size = 10,000
totalCost10000 = Dict{Tuple{Int64, Int64}, Float64}()

indexSets = load("testData_RTS/indexSets.jld2")["indexSets"]
paramOPF = load("testData_RTS/paramOPF.jld2")["paramOPF"]
paramDemand = load("testData_RTS/paramDemand.jld2")["paramDemand"]
indexSets.Ω = [1:10000...]
Ω_rv = load("testData_RTS/Ω_rv10000.jld2")["Ω_rv"] 

for k in 1:1
  for scenarioSize in [10000]
      indexSets.Ω = [1:scenarioSize...]
      prob = Dict{Int64, Float64}();
      for ω in 1:scenarioSize 
        prob[ω] = 1/scenarioSize;
      end

      state_variable = zstar.first_state_variable;
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
                                            .>= sum(paramDemand.demand[t][d] * x[d, t] for d in Dᵢ[i]));
      
      ## constraint 1f
      @constraint(Q, [g in G, t in 1:T], s[g, t] >= paramOPF.smin[g] * state_variable[:zg][g, t] );
      @constraint(Q, [g in G, t in 1:T], s[g, t] <= paramOPF.smax[g] * state_variable[:zg][g, t] );

      ## constraint 1g h i j
      for i in B 
        @constraint(Q, [t in 1:T, d in Dᵢ[i]], state_variable[:zb][i, t] >= x[d, t] );
        @constraint(Q, [t in 1:T, g in Gᵢ[i]], state_variable[:zb][i, t] >= state_variable[:zg][g, t]);
        @constraint(Q, [t in 1:T, j in out_L[i]], state_variable[:zb][i, t] >= state_variable[:zl][(i, j), t] );
        @constraint(Q, [t in 1:T, j in in_L[i]], state_variable[:zb][i, t] >= state_variable[:zl][(j, i), t] );
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
        @info "$ω"
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
        c = c + prob[ω] * JuMP.objective_value(forward2Info.model);
      end
      u = state_value + c;
      totalCost10000[(k, scenarioSize)] = u;
  end
end
save("src/Experiments/ConfidenceInterval/totalCost10000.jld2", "totalCost10000", totalCost10000)
totalCost10000 = load("src/Experiments/ConfidenceInterval/totalCost10000.jld2")["totalCost10000"]

