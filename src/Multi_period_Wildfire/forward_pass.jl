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


    (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, indexSets.Ω) 
    (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L) 

    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                          "OutputFlag" => 0, 
                                          "Threads" => 1) 
                                          )
                                          
    @variable(Q, θ_angle[B, 1:T])      ## phase angle of the bus i
    @variable(Q, P[L, 1:T] >= 0)       ## real power flow on line l; elements in L is Tuple (i, j)
    @variable(Q, s[G, 1:T] >= 0)       ## real power generation at generator g
    @variable(Q, 0 <= x[D, 1:T] <= 1)  ## load shedding


    @variable(Q, zg[G, 1:T], Bin)      ## binary status indicator of generator g
    @variable(Q, zb[B, 1:T], Bin)      ## binary status indicator of bus i
    @variable(Q, zl[L, 1:T], Bin)      ## binary status indicator of line l

    @variable(Q, θ[Ω] >= θ_bound)      ## sddp, estimated value function

    # constraint 3b 3c
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
    @constraint(Q, [i in B, t in 1:T], sum(s[g, t] for g in Gᵢ[i]) + sum(P[(i, j), t] for j in out_L[i] ) .== sum(paramDemand.demand[t][d] * x[d, t] for d in Dᵢ[i]) )
    
    ## constraint 1f
    @constraint(Q, [g in G, t in 1:T], s[g, t] >= paramOPF.smin[g] * zg[g, t] )
    @constraint(Q, [g in G, t in 1:T], s[g, t] <= paramOPF.smax[g] * zg[g, t] )

    ## constraint g h i j
    @constraint(Q, [i in B, t in 1:T, d in Dᵢ[i]], zb[i, t] >= x[d, t] )
    @constraint(Q, [i in B, t in 1:T, g in Gᵢ[i]], zb[i, t] >= zg[g, t])
    @constraint(Q, [i in B, t in 1:T, j in out_L[i]], zb[i, t] >= zl[(i, j), t] )
    @constraint(Q, [i in B, t in 1:T, j in in_L[i]], zb[i, t] >= zl[(j, i), t] )


    ## constraint k l m
    @constraint(Q, [i in B, t in 1:T-1], zb[i, t] >= zb[i, t+1] )
    @constraint(Q, [g in G, t in 1:T-1], zg[g, t] >= zg[g, t+1] )
    @constraint(Q, [l in L, t in 1:T-1], zl[l, t] >= zl[l, t+1] )


    # add cut for value functions
    for ω in Ω
      cut_coefficient = cut_collection[ω]
      iter = length(keys(cut_coefficient.v))  ## iter num
      k = length(keys(cut_coefficient.v[1]))  ## scenario num

      @constraint(Q, [i in 1:iter-1, m in 1:k, ω in Ω], θ[ω] .>= cut_coefficient.v[i][m] + 
                                                sum(cut_coefficient.πb[i][m]' * zb[:, Ω_rv[ω].τ] +
                                                cut_coefficient.πg[i][m]' * zg[:, Ω_rv[ω].τ] +
                                                cut_coefficient.πl[i][m]' * zl[:, Ω_rv[ω].τ] for t in 1:T)
                                                )

    end



    
    ## objective function
    @objective(Q, Min, sum( prob[ω] * ( sum( sum(paramDemand.w[d] * paramDemand.demand[t][d] * (1 - x[d, t]) for d in D )
                                                                                  for t in 1:Ω_rv[ω].τ - 1 ) + θ[ω] )
                                                                                    for ω in Ω)  
                                                                                )

    optimize!(Q)

    state_variable = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 2}}(:zg => round.(JuMP.value.(zg)), :zb => round.(JuMP.value.(zb)), :zl => round.(JuMP.value.(zl)))
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
                                    ẑ::Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}},
                                    randomVariables::RandomVariables                          ## realization of the random time
                                    )

    (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, indexSets.Ω)
    (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L) 


    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                "OutputFlag" => 0, 
                "Threads" => 1) 
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

    ## constraint 1b 1c
    for l in L
      i = l[1]
      j = l[2]
      @constraint(Q, [t in randomVariables.τ:T], P[l, t] <= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmax * (1 - yl[l] ) ) + slack_variable_b )
      @constraint(Q, [t in randomVariables.τ:T], P[l, t] >= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmin * (1 - yl[l] ) ) + slack_variable_c )
    end

    ## constraint 3d
    @constraint(Q, [l in L, t in randomVariables.τ:T], P[l, t] >= - paramOPF.W[l] * yl[l] )
    @constraint(Q, [l in L, t in randomVariables.τ:T], P[l, t] <= paramOPF.W[l] * yl[l] )

    ## constraint 3e
    @constraint(Q, [i in B, t in randomVariables.τ:T], sum(s[g, t] for g in Gᵢ[i]) + sum(P[(i, j), t] for j in out_L[i] ) .== sum(paramDemand.demand[t][d] * x[d, t] for d in Dᵢ[i]) )

    ## constraint 3f
    @constraint(Q, [g in G, t in randomVariables.τ:T], s[g, t] >= paramOPF.smin[g] * yg[g])
    @constraint(Q, [g in G, t in randomVariables.τ:T], s[g, t] <= paramOPF.smax[g] * yg[g])

    ## constraint g h i j
    @constraint(Q, [i in B, t in randomVariables.τ:T, d in Dᵢ[i]], yb[i] >= x[d, t] )
    @constraint(Q, [i in B, g in Gᵢ[i]], yb[i] >= yg[g])
    @constraint(Q, [i in B, j in out_L[i]], yb[i] >= yl[(i, j)] )
    @constraint(Q, [i in B, j in in_L[i]], yb[i] >= yl[(j, i)] )

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

    ## constraint n
    @constraint(Q, [i in B, j in unique(randomVariables.Ibb[i])], νb[j] >= randomVariables.ub[i] * ẑ[:zb][i] )
    @constraint(Q, [i in B, j in unique(randomVariables.Ibg[i])], νg[j] >= randomVariables.ub[i] * ẑ[:zb][i] )
    @constraint(Q, [i in B, j in unique(randomVariables.Ibl[i])], νl[j] >= randomVariables.ub[i] * ẑ[:zb][i] )

    @constraint(Q, [i in G, j in unique(randomVariables.Igb[i])], νb[j] >= randomVariables.ug[i] * ẑ[:zg][i] )
    @constraint(Q, [i in G, j in unique(randomVariables.Igg[i])], νg[j] >= randomVariables.ug[i] * ẑ[:zg][i] )
    @constraint(Q, [i in G, j in unique(randomVariables.Igl[i])], νl[j] >= randomVariables.ug[i] * ẑ[:zg][i] )

    @constraint(Q, [i in L, j in unique(randomVariables.Ilb[i])], νb[j] >= randomVariables.ul[i] * ẑ[:zl][i] )
    @constraint(Q, [i in L, j in unique(randomVariables.Ilg[i])], νg[j] >= randomVariables.ul[i] * ẑ[:zl][i] )
    @constraint(Q, [i in L, j in unique(randomVariables.Ill[i])], νl[j] >= randomVariables.ul[i] * ẑ[:zl][i] )


    ## objective function
    @objective(Q, Min,  
            sum( sum(paramDemand.w[d] * paramDemand.demand[t][d] * (1 - x[d, t]) for d in D ) for t in randomVariables.τ:T) +
            sum(paramDemand.cb[i] * νb[i] for i in B) + 
            sum(paramDemand.cg[g] * νg[g] for g in G) + 
            sum(paramDemand.cl[l] * νl[l] for l in L) + paramDemand.penalty * slack_variable_b - paramDemand.penalty * slack_variable_c
            )

    optimize!(Q)
    state_variable = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}(:zg => round.(JuMP.value.(yg)), :zb => round.(JuMP.value.(yb)), :zl => round.(JuMP.value.(yl)))
    state_value    = JuMP.objective_value(Q)

    return [state_variable, state_value]  ## returen [Lt, y, θ, f]
end

