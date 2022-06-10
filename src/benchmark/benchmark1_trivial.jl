function benchmark_stage1_model!(; indexSets::IndexSets = indexSets, 
                                    paramDemand::ParamDemand = paramDemand, 
                                    paramOPF::ParamOPF = paramOPF)


    (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, indexSets.Ω) 
    (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L) 

    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                          "OutputFlag" => 0, 
                                          "Threads" => 0) 
                                          )
                                          
    @variable(Q, θ_angle[B, 1:T])      ## phase angle of the bus i
    @variable(Q, P[L, 1:T] >= 0)       ## real power flow on line l; elements in L is Tuple (i, j)
    @variable(Q, s[G, 1:T] >= 0)       ## real power generation at generator g
    @variable(Q, 0 <= x[D, 1:T] <= 1)  ## load shedding


    @variable(Q, zg[G, 1:T], Bin)      ## binary status indicator of generator g
    @variable(Q, zb[B, 1:T], Bin)      ## binary status indicator of bus i
    @variable(Q, zl[L, 1:T], Bin)      ## binary status indicator of line l

    # constraint 1b 1c
    for l in L
      i = l[1]
      j = l[2]
      @constraint(Q, [t in 1:T], P[l, t] <= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmax * (1 - zl[l, t] ) ) )
      @constraint(Q, [t in 1:T], P[l, t] >= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmin * (1 - zl[l, t] ) ) )
    end

    ## constraint 1d
    @constraint(Q, [l in L, t in 1:T], P[l, t] >= - paramOPF.W[l] * zl[l, t] )
    @constraint(Q, [l in L, t in 1:T], P[l, t] <= paramOPF.W[l] * zl[l, t] )

    ## constraint 1e
    @constraint(Q, [i in B, t in 1:T], sum(s[g, t] for g in Gᵢ[i]) + 
                                          sum(P[(i, j), t] for j in out_L[i]) - 
                                              sum(P[(j, i), t] for j in in_L[i]) 
                                                .== sum(paramDemand.demand[t][d] * x[d, t] for d in Dᵢ[i]) )
    
    ## constraint 1f
    @constraint(Q, [g in G, t in 1:T], s[g, t] >= paramOPF.smin[g] * zg[g, t] )
    @constraint(Q, [g in G, t in 1:T], s[g, t] <= paramOPF.smax[g] * zg[g, t] )

    ## constraint 1g h i j
    for i in B 
      @constraint(Q, [t in 1:T, d in Dᵢ[i]], zb[i, t] >= x[d, t] )
      @constraint(Q, [t in 1:T, g in Gᵢ[i]], zb[i, t] >= zg[g, t])
      @constraint(Q, [t in 1:T, j in out_L[i]], zb[i, t] >= zl[(i, j), t] )
      @constraint(Q, [t in 1:T, j in in_L[i]], zb[i, t] >= zl[(j, i), t] )
    end


    ## constraint 1k
    @constraint(Q, [i in B, t in 1:T-1], zb[i, t] >= zb[i, t+1] )
    @constraint(Q, [g in G, t in 1:T-1], zg[g, t] >= zg[g, t+1] )
    @constraint(Q, [l in L, t in 1:T-1], zl[l, t] >= zl[l, t+1] )
    
    ## objective function
    @objective(Q, Min, sum( 
                            sum( paramDemand.w[d] * (1 - x[d, t]) for d in D )
                                                                                                 
                        for t in 1:T)
                                                                                    
                )                                                                    
    ###################################################### solve the model and display the result ###########################################################  
    optimize!(Q)

    state_variable = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 2}}(
                                                                                :zg => round.(JuMP.value.(zg)), 
                                                                                :zb => round.(JuMP.value.(zb)), 
                                                                                :zl => round.(JuMP.value.(zl))
                            )
    local_variable = JuMP.value.(x)

    return (state_variable = state_variable, 
                    local_variable = local_variable)  ## returen [state_variable, first_stage value, objective_value(Q)]
end



## just the second stage model, have no any modification
function benchmark_stage2_model!(ẑ::Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}},                                     
                                    randomVariables::RandomVariables ;                         ## realization of the random time
                                    indexSets::IndexSets = indexSets, 
                                    paramDemand::ParamDemand = paramDemand, 
                                    paramOPF::ParamOPF = paramOPF, 
                                    )

    (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, indexSets.Ω)
    (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L) 


    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                "OutputFlag" => 0, 
                "Threads" => 0) 
                )

    @variable(Q, θ_angle[B, 1:T])      ## phase angle of the bus i
    @variable(Q, P[L, 1:T] >= 0)       ## real power flow on line l; elements in L is Tuple (i, j)
    @variable(Q, s[G, 1:T] >= 0)       ## real power generation at generator g
    @variable(Q, 0 <= x[D, 1:T] <= 1)  ## load shedding

    @variable(Q, yb[B], Bin)
    @variable(Q, yg[G], Bin)
    @variable(Q, yl[L], Bin)

    @variable(Q, νb[B], Bin)
    @variable(Q, νg[G], Bin)
    @variable(Q, νl[L], Bin)

    @variable(Q, slack_variable_b >= 0)
    @variable(Q, slack_variable_c <= 0)    

    ## constraint k l m 
    @constraint(Q, [i in B], yb[i] <= ẑ[:zb][i] )
    @constraint(Q, [g in G], yg[g] <= ẑ[:zg][g] )
    @constraint(Q, [l in L], yl[l] <= ẑ[:zl][l] )

    @constraint(Q, [i in B], yb[i] <= 1- νb[i] )
    @constraint(Q, [g in G], yg[g] <= 1- νg[g] )
    @constraint(Q, [l in L], yl[l] <= 1- νl[l] )

    @constraint(Q, [i in B], νb[i] >= randomVariables.vb[i] )
    @constraint(Q, [g in G], νg[g] >= randomVariables.vg[g] )
    @constraint(Q, [l in L], νl[l] >= randomVariables.vl[l] )

   
    for i in B 
      ## constraint 3e
      @constraint(Q, [t in randomVariables.τ:T], sum(s[g, t] for g in Gᵢ[i]) + 
                                                        sum(P[(i, j), t] for j in out_L[i]) - 
                                                            sum(P[(j, i), t] for j in in_L[i]) 
                                                                .== sum(paramDemand.demand[t][d] * x[d, t] for d in Dᵢ[i]) )

      ## constraint g h i j
      @constraint(Q, [t in randomVariables.τ:T, d in Dᵢ[i]], yb[i] >= x[d, t] )
      @constraint(Q, [g in Gᵢ[i]], yb[i] >= yg[g])
      @constraint(Q, [j in out_L[i]], yb[i] >= yl[(i, j)] )
      @constraint(Q, [j in in_L[i]], yb[i] >= yl[(j, i)] )

      ## constraint n
      @constraint(Q, [j in unique(randomVariables.Ibb[i])], νb[j] >= randomVariables.ub[i] * ẑ[:zb][i] )
      @constraint(Q, [j in unique(randomVariables.Ibg[i])], νg[j] >= randomVariables.ub[i] * ẑ[:zb][i] )
      @constraint(Q, [j in unique(randomVariables.Ibl[i])], νl[j] >= randomVariables.ub[i] * ẑ[:zb][i] )
    end

    for g in G 
      ## constraint n
      @constraint(Q, [j in unique(randomVariables.Igb[g])], νb[j] >= randomVariables.ug[g] * ẑ[:zg][g] )
      @constraint(Q, [j in unique(randomVariables.Igg[g])], νg[j] >= randomVariables.ug[g] * ẑ[:zg][g] )
      @constraint(Q, [j in unique(randomVariables.Igl[g])], νl[j] >= randomVariables.ug[g] * ẑ[:zg][g] )

       ## constraint 3f
      @constraint(Q, [t in randomVariables.τ:T], s[g, t] >= paramOPF.smin[g] * yg[g])
      @constraint(Q, [t in randomVariables.τ:T], s[g, t] <= paramOPF.smax[g] * yg[g])
    end

    for l in L 
      ## constraint n
      @constraint(Q, [j in unique(randomVariables.Ilb[l])], νb[j] >= randomVariables.ul[l] * ẑ[:zl][l] )
      @constraint(Q, [j in unique(randomVariables.Ilg[l])], νg[j] >= randomVariables.ul[l] * ẑ[:zl][l] )
      @constraint(Q, [j in unique(randomVariables.Ill[l])], νl[j] >= randomVariables.ul[l] * ẑ[:zl][l] )


      ## constraint 3b 3c
      i = l[1]
      j = l[2]
      @constraint(Q, [t in randomVariables.τ:T], P[l, t] <= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmax * (1 - yl[l] ) ) + slack_variable_b )
      @constraint(Q, [t in randomVariables.τ:T], P[l, t] >= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmin * (1 - yl[l] ) ) + slack_variable_c )

      ## constraint 3d
      @constraint(Q, [t in randomVariables.τ:T], P[l, t] >= - paramOPF.W[l] * yl[l] )
      @constraint(Q, [t in randomVariables.τ:T], P[l, t] <= paramOPF.W[l] * yl[l] )
    end 

    ## objective function
    @objective(Q, Min,  
            sum( sum(paramDemand.w[d] * (1 - x[d, t]) for d in D ) for t in randomVariables.τ:T) +
            sum(paramDemand.cb[i] * νb[i] for i in B) + 
            sum(paramDemand.cg[g] * νg[g] for g in G) + 
            sum(paramDemand.cl[l] * νl[l] for l in L) + paramDemand.penalty * slack_variable_b - paramDemand.penalty * slack_variable_c
            )
    ####################################################### solve the model and display the result ###########################################################
    optimize!(Q)
    state_variable = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}(:zg => round.(JuMP.value.(yg)), :zb => round.(JuMP.value.(yb)), :zl => round.(JuMP.value.(yl)))
    state_value    = JuMP.objective_value(Q)

    return (state_variable = state_variable, 
                        state_value = state_value)
end




function benchmarkTrivial!(; Ω_rv::Dict{Int64, RandomVariables} = Ω_rv, 
                        indexSets::IndexSets = indexSets, 
                        paramDemand::ParamDemand = paramDemand, 
                        paramOPF::ParamOPF = paramOPF)
    (state1_variable, 
            local1_variable) = benchmark_stage1_model!();
    
    costTrivial = Dict{Int64, Float64}()
    for ω in indexSets.Ω
        randomVariables = Ω_rv[ω];

        ẑ = Dict(   :zg => state1_variable[:zg][:, randomVariables.τ - 1], 
                    :zb => state1_variable[:zb][:, randomVariables.τ - 1], 
                    :zl => state1_variable[:zl][:, randomVariables.τ - 1]
                    );
        ## modify the constraints according to the first stage state variables
        (state2_variable, 
                        state_value) = benchmark_stage2_model!(ẑ, randomVariables);
        costTrivial[ω] = sum( sum( paramDemand.w[d] * paramDemand.demand[t][d] * (1 - local1_variable[d, t]) for d in indexSets.D ) for t in 1:randomVariables.τ - 1) + 
                                                                                state_value
    end
    return costTrivial
end


