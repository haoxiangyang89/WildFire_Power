## just the second stage model, have no any modification
function banchmark2_oracle!(randomVariables::RandomVariables;                          ## realization of the random time
                                    indexSets::IndexSets = indexSets, 
                                    paramDemand::ParamDemand = paramDemand, 
                                    paramOPF::ParamOPF = paramOPF, 
                                    )

    (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, indexSets.Ω);
    (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L);
    τ = randomVariables.τ;

    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                "OutputFlag" => 0, 
                "Threads" => 0) 
                )

    @variable(Q, θ_angle[B, 1:T])      ## phase angle of the bus i
    @variable(Q, P[L, 1:T] >= 0)       ## real power flow on line l; elements in L is Tuple (i, j)
    @variable(Q, s[G, 1:T] >= 0)       ## real power generation at generator g
    @variable(Q, 0 <= x[D, 1:T] <= 1)  ## load shedding

    @variable(Q, zb[B, 1:T], Bin)
    @variable(Q, zg[G, 1:T], Bin)
    @variable(Q, zl[L, 1:T], Bin)

    @variable(Q, νb[B], Bin)
    @variable(Q, νg[G], Bin)
    @variable(Q, νl[L], Bin)

    @variable(Q, slack_variable_b >= 0)
    @variable(Q, slack_variable_c <= 0)    

    ## constraint 5k l m 
    @constraint(Q, [t in 1:T-1, i in B], zb[i, t] >= zb[i, t+1] )
    @constraint(Q, [t in 1:T-1, g in G], zg[g, t] >= zg[g, t+1] )
    @constraint(Q, [t in 1:T-1, l in L], zl[l, t] >= zl[l, t+1] )

    @constraint(Q, [i in B], zb[i, τ] <= 1- νb[i] )
    @constraint(Q, [g in G], zg[g, τ] <= 1- νg[g] )
    @constraint(Q, [l in L], zl[l, τ] <= 1- νl[l] )

    @constraint(Q, [i in B], νb[i] >= randomVariables.vb[i] )
    @constraint(Q, [g in G], νg[g] >= randomVariables.vg[g] )
    @constraint(Q, [l in L], νl[l] >= randomVariables.vl[l] )

   
    for i in B 
      ## constraint 5e
      @constraint(Q, [t in 1:T], sum(s[g, t] for g in Gᵢ[i]) + 
                                                        sum(P[(i, j), t] for j in out_L[i]) - 
                                                            sum(P[(j, i), t] for j in in_L[i]) 
                                                              .== sum(paramDemand.demand[t][d] * x[d, t] for d in Dᵢ[i]) )

      ## constraint 5g h i j
      @constraint(Q, [t in 1:T, d in Dᵢ[i]], zb[i, t] >= x[d, t] )
      @constraint(Q, [t in 1:T, g in Gᵢ[i]], zb[i, t] >= zg[g, t])
      @constraint(Q, [t in 1:T, j in out_L[i]], zb[i, t] >= zl[(i, j), t] )
      @constraint(Q, [t in 1:T, j in in_L[i]], zb[i, t] >= zl[(j, i), t] )

      ## constraint n
      @constraint(Q, [j in unique(randomVariables.Ibb[i])], νb[j] >= randomVariables.ub[i] * zb[i, τ] )
      @constraint(Q, [j in unique(randomVariables.Ibg[i])], νg[j] >= randomVariables.ub[i] * zb[i, τ] )
      @constraint(Q, [j in unique(randomVariables.Ibl[i])], νl[j] >= randomVariables.ub[i] * zb[i, τ] )
    end

    for g in G 
      ## constraint n
      @constraint(Q, [j in unique(randomVariables.Igb[g])], νb[j] >= randomVariables.ug[g] * zg[g, τ] )
      @constraint(Q, [j in unique(randomVariables.Igg[g])], νg[j] >= randomVariables.ug[g] * zg[g, τ] )
      @constraint(Q, [j in unique(randomVariables.Igl[g])], νl[j] >= randomVariables.ug[g] * zg[g, τ] )

       ## constraint 5f
      @constraint(Q, [t in 1:T], s[g, t] >= paramOPF.smin[g] * zg[g, t])
      @constraint(Q, [t in 1:T], s[g, t] <= paramOPF.smax[g] * zg[g, t])
    end

    for l in L 
      ## constraint n
      @constraint(Q, [j in unique(randomVariables.Ilb[l])], νb[j] >= randomVariables.ul[l] * zl[l, τ] )
      @constraint(Q, [j in unique(randomVariables.Ilg[l])], νg[j] >= randomVariables.ul[l] * zl[l, τ] )
      @constraint(Q, [j in unique(randomVariables.Ill[l])], νl[j] >= randomVariables.ul[l] * zl[l, τ] )


      ## constraint 5b 5c
      i = l[1]
      j = l[2]
      @constraint(Q, [t in 1:T], P[l, t] <= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmax * (1 - zl[l, t] ) ) + slack_variable_b )
      @constraint(Q, [t in 1:T], P[l, t] >= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmin * (1 - zl[l, t] ) ) + slack_variable_c )

      ## constraint 5d
      @constraint(Q, [t in 1:T], P[l, t] >= - paramOPF.W[l] * zl[l, t] )
      @constraint(Q, [t in 1:T], P[l, t] <= paramOPF.W[l] * zl[l, t] )
    end 

    ## objective function
    @objective(Q, Min,  
            sum( sum(paramDemand.w[d] * (1 - x[d, t]) for d in D ) for t in 1:T) +
            sum(paramDemand.cb[i] * νb[i] for i in B) + 
            sum(paramDemand.cg[g] * νg[g] for g in G) + 
            sum(paramDemand.cl[l] * νl[l] for l in L) + paramDemand.penalty * slack_variable_b - paramDemand.penalty * slack_variable_c
            )
    ####################################################### solve the model and display the result ###########################################################
    optimize!(Q)
    # JuMP.value.(slack_variable_b)
    # JuMP.value.(slack_variable_c)
    # st = termination_status(Q)
    state_variable = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 2}}(:zg => round.(JuMP.value.(zg)), :zb => round.(JuMP.value.(zb)), :zl => round.(JuMP.value.(zl)))
    state_value    = JuMP.objective_value(Q)

    return (state_variable = state_variable, 
                        state_value = state_value) 
end


function benchmarkOracle!(; Ω_rv::Dict{Int64, RandomVariables} = Ω_rv, 
                        indexSets::IndexSets = indexSets, 
                        paramDemand::ParamDemand = paramDemand, 
                        paramOPF::ParamOPF = paramOPF)
        
    costOracle = Dict{Int64, Float64}()
    for ω in indexSets.Ω
        randomVariables = Ω_rv[ω];

        (state_variable, 
                            state_value) = banchmark2_oracle!(randomVariables);

        costOracle[ω] = state_value
    end



    return costOracle

end

