## This test is design to figure out the cost when put in the solutions generated by different scenario sizes.
#  That is, scenario size = [10, 20, 50, 80, 100, 120], we obtain the first-stage decision for each case
#  And test this decision to a new data set scenario size = 100

## -----------------------------------  Generate the first-stage decision  ---------------------------------------- ##
gurobiResultList = Dict{Any, Any}()

indexSets = load("testData_RTS/indexSets.jld2")["indexSets"]
paramOPF = load("testData_RTS/paramOPF.jld2")["paramOPF"]
paramDemand = load("testData_RTS/paramDemand.jld2")["paramDemand"]
Ω_rvList = load("testData_RTS/Ω_rvList.jld2")["Ω_rvList"]

Ω = 50; 
pList = [.5, .9, 0.95, 0.99];
probList = Dict{Float64, Any}()
for p in pList
    prob = Dict{Int64, Float64}();
    for ω in 1:Ω
        prob[ω] = (1 - p)/Ω;
    end
    probList[p] = prob
end

for p in pList
  for i in 1:20 
    @info "$p, $i"
    Ω_rv = Ω_rvList[(Ω, i)]; indexSets.Ω = [1:Ω...];

    # timelimit = Ω >= 100 ? 14400 : 6000
    gurobiResultList[(p, i)] = gurobiOptimize!(indexSets, 
                                        paramDemand, 
                                        paramOPF, 
                                        Ω_rv,
                                        probList[p]; 
                                        mipGap = 2e-2, timelimit = 1800); 
  end
  save("src/Experiments/WeightTest/gurobiResultList.jld2", "gurobiResultList", gurobiResultList)
end

save("src/Experiments/WeightTest/gurobiResultList.jld2", "gurobiResultList", gurobiResultList)

# gurobiResultList = load("src/Experiments/WeightTest/gurobiResultList.jld2")["gurobiResultList"]



## ------------------------------  Compute total costs by using different decisions  ----------------------------- ##
## scenario size = 200
totalCost = Dict{Tuple{Int64, Float64}, Float64}()

indexSets = load("testData_RTS/indexSets.jld2")["indexSets"]
paramOPF = load("testData_RTS/paramOPF.jld2")["paramOPF"]
paramDemand = load("testData_RTS/paramDemand.jld2")["paramDemand"]
Ω_rv = load("testData_RTS/Ω_rv5000.jld2")["Ω_rv"]
indexSets.Ω = [1:5000...]
prob = Dict{Int64, Float64}();
for ω in 1:5000 
    prob[ω] = 1/5000;
end

gurobiResultList = load("src/Experiments/WeightTest/gurobiResultList.jld2")["gurobiResultList"]

for k in 1:20 
  for probEmpty in [.5, .9, 0.95, 0.99]

      state_variable = gurobiResultList[(probEmpty, k)].first_state_variable;
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
        @info "$probEmpty, $k, $ω"
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
      totalCost[(k, probEmpty)] = u;
  end
  save("src/Experiments/WeightTest/totalCost.jld2", "totalCost", totalCost)
end
save("src/Experiments/WeightTest/totalCost.jld2", "totalCost", totalCost)
totalCost = load("src/Experiments/WeightTest/totalCost.jld2")["totalCost"]

