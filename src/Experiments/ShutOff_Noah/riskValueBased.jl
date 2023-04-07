indexSets = load("testData_RTS_New/indexSets.jld2")["indexSets"]
paramOPF = load("testData_RTS_New/paramOPF.jld2")["paramOPF"]
paramDemand = load("testData_RTS_New/paramDemand.jld2")["paramDemand"]
riskValueDict = load("src/Experiments/ShutOff_Noah/riskValueDict.jld2")["riskValueDict"];



## ------------------------------ Multi-period Risk Value-based Model ------------------------------------------ ## 
function Risk_Value_Based_model!( ; indexSets::IndexSets = indexSets, 
                                    paramDemand::ParamDemand = paramDemand, 
                                    paramOPF::ParamOPF = paramOPF, 
                                    riskValueDict::Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 2}} = riskValueDict,
                                    outputFlag::Int64 = 0, timelimit::Real = 120, α::Float64 = .1)

    (D, G, L, B, T) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T) 
    (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L) 

    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                          "OutputFlag" => outputFlag, 
                                          "Threads" => 0, 
                                          "MIPGap" => 1e-3, 
                                          "TimeLimit" => timelimit)
                                          ) 
                                          
    @variable(Q, θ_angle[B, 1:T])      ## phase angle of the bus i
    @variable(Q, P[L, 1:T])            ## real power flow on line l; elements in L is Tuple (i, j)
    @variable(Q, s[G, 1:T] ≥ 0)            ## real power generation at generator g
    @variable(Q, 0 ≤ x[D, 1:T] ≤ 1)    ## load shedding


    @variable(Q, zg[G, 1:T], Bin)      ## binary status indicator of generator g
    @variable(Q, zb[B, 1:T], Bin)      ## binary status indicator of bus i
    @variable(Q, zl[L, 1:T], Bin)      ## binary status indicator of line l

    # constraint 3b 3c
    for l in L
      i = l[1]
      j = l[2]
      @constraint(Q, [t in 1:T], P[l, t] ≤ - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmax * (1 - zl[l, t] ) ) )
      @constraint(Q, [t in 1:T], P[l, t] ≥ - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmin * (1 - zl[l, t] ) ) )
    end

    ## constraint 1d
    @constraint(Q, [l in L, t in 1:T], P[l, t] ≥ - paramOPF.W[l] * zl[l, t] )
    @constraint(Q, [l in L, t in 1:T], P[l, t] ≤ paramOPF.W[l] * zl[l, t] )

    ## constraint 1e
    @constraint(Q, [i in B, t in 1:T], sum(s[g, t] for g in Gᵢ[i]) -
                                          sum(P[(i, j), t] for j in out_L[i]) + 
                                            sum(P[(j, i), t] for j in in_L[i]) 
                                              .== sum(paramDemand.demand[t][d] * x[d, t] for d in Dᵢ[i]) )
    
    ## constraint 1f
    @constraint(Q, [g in G, t in 1:T], s[g, t] ≥ paramOPF.smin[g] * zg[g, t] )
    @constraint(Q, [g in G, t in 1:T], s[g, t] ≤ paramOPF.smax[g] * zg[g, t] )

    ## constraint 1g h i j
    for i in B 
      @constraint(Q, [t in 1:T, d in Dᵢ[i]], zb[i, t] ≥ x[d, t] )
      @constraint(Q, [t in 1:T, g in Gᵢ[i]], zb[i, t] ≥ zg[g, t])
      @constraint(Q, [t in 1:T, j in out_L[i]], zb[i, t] ≥ zl[(i, j), t] )
      @constraint(Q, [t in 1:T, j in in_L[i]], zb[i, t] ≥ zl[(j, i), t] )
    end


    ## constraint 1k
    @constraint(Q, [i in B, t in 1:T-1], zb[i, t] ≥ zb[i, t+1] )
    @constraint(Q, [g in G, t in 1:T-1], zg[g, t] ≥ zg[g, t+1] )
    @constraint(Q, [l in L, t in 1:T-1], zl[l, t] ≥ zl[l, t+1] )


    ## Risk Value 
    @variable(Q, Risk[1:T] ≥ 0) 
    # for each period
    @constraint(Q, [t in 1:T], Risk[t] ≥ sum(zl[l, t] * riskValueDict[:zl][l, t] / 10 for l in L) 
                                                    )

    
    ## objective function
    @objective(Q, Min, sum( α * Risk[t] - (1 - α) * sum(paramDemand.w[d] * x[d, t] for d in D )
                                                                                                        for t in 1:T )  
              )
                                                                    
    ###################################################### solve the model and display the result ###########################################################  
    optimize!(Q)

    state_variable = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 2}}(:zg => round.(JuMP.value.(zg)), 
                                                                                :zb => round.(JuMP.value.(zb)), 
                                                                                :zl => round.(JuMP.value.(zl)))

    return (state_variable = state_variable, 
              obj_value = JuMP.objective_value(Q))  ## returen [state_variable, first_stage value, objective_value(Q)]

end

riskValueBasedPlan = Dict()
for α in [0.0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5]
    riskValueBasedPlan[α] = Risk_Value_Based_model!(;α = α, outputFlag = 1, timelimit = 500)
end

save("src/Experiments/ShutOff_Noah/riskValueBasedPlan.jld2", "riskValueBasedPlan", riskValueBasedPlan)
