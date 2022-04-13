#############################################################################################
###################################  function: forward pass #################################
#############################################################################################
#     S = [(1, 1, 1), (N, N, N)]
# @variable(model, x1[i=1:N, j=1:N, k=1:N; (i, j, k) in S])

"""     forward_stage1_optimize!
1. Solve the forward problem Stage 1
2. return the decision value -- z and its state value

"""


function forward_stage1_optimize!(indexSets::IndexSets, 
                                    paramDemand::ParamDemand, 
                                    paramOPF::ParamOPF, 
                                    Ω_rv::Dict{Int64, RandomVariables},
                                    prob::Dict{Int64, Float64},
                                    cut_collection::Dict{Int64, CutCoefficient};  ## the index is ω
                                    θ_bound::Real = 0.0)


    (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, IndexSets.Ω) 
    (_D, _G, in_L, out_L) = (indexSets._D, indexSets._G, indexSets.in_L, indexSets.out_L) 

    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                          "OutputFlag" => 0, 
                                          "Threads" => 1) 
                                          )

    @variable(Q, θ_angle[B, 1:T]) 
    @variable(Q, P[L, 1:T] >= 0) ## elements in L is Tuple (i, j)
    @variable(Q, s[G, 1:T] >= 0)
    @variable(Q, 0 <= x[1:T, D] <= 1)

    @variable(Q, zg[G, 1:T], Bin)
    @variable(Q, zb[B, 1:T], Bin)
    @variable(Q, zl[L, 1:T], Bin)

    @variable(Q, θ[Ω] >= θ_bound)

    ## constraint 1b 1c
    for l in L
      i = l[1]
      j = l[2]
      @constraint(Q, [t in 1:T], P[l, t] <= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + ParamOPF.θmax * (1 - zl[l, t] ) ) )
      @constraint(Q, [t in 1:T], P[l, t] >= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + ParamOPF.θmin * (1 - zl[l, t] ) ) )
    end

    ## constraint 1d
    @constraint(Q, [l in L, t in 1:T], - ParamOPF.W[l] * zl[l, t] <= P[l, t] <= ParamOPF.W[l] * zl[l, t] )

    ## constraint 1e
    @constraint(Q, [i in B, t in 1:T], sum(s[g, t] for g in _G[i]) + sum(P[(i, j), t] for j in out_L[i] ) .== sum(paramDemand.demand[t][d] * x[t, d] for d in _D[i]) )
    
    ## constraint 1f
    @constraint(Q, [g in G, t in 1:T], ParamOPF.smin * zg[g, t] <= s[g, t] <= ParamOPF.smax * zg[g, t] )

    ## constraint g h i j
    for i in B
      if _D[i][1] != 0
        @constraint(Q, [t in 1:T, d in _D[i]], zb[i, t] >= x[t, d] )
      end
      if _G[i][1] != 0
        @constraint(Q, [t in 1:T, d in _G[i]], zb[i, t] >= zg[g, t])
      end
      
      if out_L[i][1] != 0
        @constraint(Q, [t in 1:T, j in out_L[i]], zb[i, t] >= zl[(i, j), t] )
      end

      if in_L[i][1] != 0
        @constraint(Q, [t in 1:T, j in in_L[i]], zb[i, t] >= zl[(j, i), t] )
      end
    end

    ## constraint k l m
    @constraint(Q, [i in B, t in 1:T-1], zb[i, t] >= zb[i, t+1] )
    @constraint(Q, [g in G, t in 1:T-1], zg[g, t] >= zg[g, t+1] )
    @constraint(Q, [l in L, t in 1:T-1], zl[l, t] >= zl[l, t+1] )


    ## add cut for value functions
    for ω in Ω
      cut_coefficient = cut_collection[ω]
      iter = length(keys(cut_coefficient.v))  ## iter num
      k = length(keys(cut_coefficient.v[1]))  ## scenario num

      @constraint(Q, cut[i in 1:iter-1, m in 1:k in Ω], θ[ω] >= cut_coefficient.v[i][m] + 
                                                cut_coefficient.πb[i][m]' * zb +
                                                cut_coefficient.πg[i][m]' * zg +
                                                cut_coefficient.πl[i][m]' * zl
                                                )

    end



    
    ## objective function
    @objecive(Q, Min, sum( prob[ω] * ( sum( sum(paramDemand.w[d] * paramDemand.demand[t, d] * (1 - x[t, d]) for d in D ) for t in 1:Ω_rv[ω].τ - 1 ) + θ[ω] ) for ω in Ω)  )

    optimize!(Q)

    state_variable = Dict{:Symbol, Vector{Int64}}(:zg => round.(JuMP.value.(zg)), :zb => round.(JuMP.value.(zb)), :zl => round.(JuMP.value.(zl)))
    state_value    = JuMP.objective_value(Q) - sum(prob[ω] * JuMP.value(θ[ω]) for ω in Ω)       ## 1a first term

    return [state_variable, state_value, JuMP.objective_value(Q)]  ## returen [Lt, y, θ, f]
end




"""     forward_stage2_optimize!
1. Solve the forward problem Stage 2
2. return the decision value -- y and its state value

"""
function forward_stage2_optimize!(indexSets::IndexSets, 
                                    paramDemand::ParamDemand, 
                                    paramOPF::ParamOPF, 
                                    ẑ::Dict{Symbol, Vector{Int64}},
                                    randomVariables::RandomVariables                          ## realization of the random time
                                    )

    (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, IndexSets.Ω)
    (_D, _G, in_L, out_L) = (indexSets._D, indexSets._G, indexSets.in_L, indexSets.out_L) 


    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                "OutputFlag" => 0, 
                "Threads" => 1) 
                )

    @variable(Q, θ_angle[B, 1:T]) 
    @variable(Q, P[L, 1:T] >= 0) ## elements in L is Tuple (i, j)
    @variable(Q, s[G, 1:T] >= 0)
    @variable(Q, 0 <= x[1:T, D] <= 1)

    @variable(Q, yb[B], Bin)
    @variable(Q, yg[G], Bin)
    @variable(Q, yl[L], Bin)

    @variable(Q, νb[B], Bin)
    @variable(Q, νg[G], Bin)
    @variable(Q, νl[L], Bin)

    ## constraint 3b 3c
    for l in L
      i = l[1]
      j = l[2]
      @constraint(Q, [t in randomVariables.τ:T], P[l, t] <= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + ParamOPF.θmax * (1 - yl[l] ) ) )
      @constraint(Q, [t in randomVariables.τ:T], P[l, t] >= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + ParamOPF.θmin * (1 - yl[l] ) ) )
    end

    ## constraint 3d
    @constraint(Q, [l in L, t in randomVariables.τ:T], - ParamOPF.W[l] * yl[l] <= P[l, t] <= ParamOPF.W[l] * yl[l] )

    ## constraint 3e
    @constraint(Q, [i in B, t in randomVariables.τ:T], sum(s[g, t] for g in _G[i]) + sum(P[(i, j), t] for j in out_L[i] ) .== sum(paramDemand.demand[t][d] * x[t, d] for d in _D[i]) )

    ## constraint 3f
    @constraint(Q, [g in G, t in randomVariables.τ:T], ParamOPF.smin * yg[g] <= s[g, t] <= ParamOPF.smax * yg[g])

    ## constraint g h i j
    for i in B
      if _D[i][1] != 0
        @constraint(Q, [t in randomVariables.τ:T, d in _D[i]], yb[i] >= x[t, d])
      end
      if _G[i][1] != 0
        @constraint(Q, [g in _G[i]], yb[i] >= yg[g])
      end
      if out_L[i][1] != 0
        @constraint(Q, [j in out_L[i]], yb[i] >= yl[(i, j)])
      end
      if in_L[i][1] != 0
        @constraint(Q, [j in in_L[i]], yb[i] >= yl[(j, i)])
      end
    end
    ## constraint k l m 
    @constraint(Q, [i in B], yb[i] <= ẑ[:zb][i] )
    @constraint(Q, [g in G], yg[g] <= ẑ[:zg][g] )
    @constraint(Q, [l in L], yl[l] <= ẑ[:zl][l] )

    @constraint(Q, [i in B], yb[i] <= 1- νb[b] )
    @constraint(Q, [g in G], yg[g] <= 1- νg[g] )
    @constraint(Q, [l in L], yl[l] <= 1- νl[l] )

    @constraint(Q, [i in B], νb[i] >= randomVariables.vb[i] )
    @constraint(Q, [g in G], νg[g] >= randomVariables.vg[g] )
    @constraint(Q, [l in L], νl[l] >= randomVariables.vl[l] )

    ## constraint n
    @constraint(Q, [i in B, j in randomVariables.Ibb[i]], νb[j] >= randomVariables.ub[i] * ẑ[:zb][i] )
    @constraint(Q, [i in B, j in randomVariables.Ibg[i]], νg[j] >= randomVariables.ub[i] * ẑ[:zb][i] )
    @constraint(Q, [i in B, j in randomVariables.Ibl[i]], νl[j] >= randomVariables.ub[i] * ẑ[:zb][i] )

    @constraint(Q, [i in G, j in randomVariables.Igb[i]], νb[j] >= randomVariables.ug[i] * ẑ[:zg][i] )
    @constraint(Q, [i in G, j in randomVariables.Igg[i]], νg[j] >= randomVariables.ug[i] * ẑ[:zg][i] )
    @constraint(Q, [i in G, j in randomVariables.Igl[i]], νl[j] >= randomVariables.ug[i] * ẑ[:zg][i] )

    @constraint(Q, [i in L, j in randomVariables.Ilb[i]], νb[j] >= randomVariables.ul[i] * ẑ[:zl][i] )
    @constraint(Q, [i in L, j in randomVariables.Ilg[i]], νg[j] >= randomVariables.ul[i] * ẑ[:zl][i] )
    @constraint(Q, [i in L, j in randomVariables.Ill[i]], νl[j] >= randomVariables.ul[i] * ẑ[:zl][i] )


    ## objective function
    @objecive(Q, Min,  
            sum( sum(paramDemand.w[d] * paramDemand.demand[t, d] * (1 - x[t, d]) for d in D ) for t in randomVariables.τ:T) +
            sum(ParamDemand.cb[i] * νb[i] for i in B) + 
            sum(ParamDemand.cg[g] * νg[g] for g in G) + 
            sum(ParamDemand.cl[l] * νl[l] for l in L)
            )

    optimize!(Q)
    state_variable = Dict{:Symbol, Vector{Int64}}(:zg => round.(JuMP.value.(zg)), :zb => round.(JuMP.value.(zb)), :zl => round.(JuMP.value.(zl)))
    state_value    = JuMP.objective_value(Q)

    return [state_variable, state_value]  ## returen [Lt, y, θ, f]
end