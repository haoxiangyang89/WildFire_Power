#############################################################################################
###################################  function: forward pass #################################
#############################################################################################
#     S = [(1, 1, 1), (N, N, N)]
# @variable(model, x1[i=1:N, j=1:N, k=1:N; (i, j, k) in S])

"""     forward_stage1_optimize!
1. Solve the forward problem Stage 1
2. return the decision value -- z and its state value

"""


function forward_stage1_model!(indexSets::IndexSets, 
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
                                          "Threads" => 0) 
                                          )
                                          
    @variable(Q, θ_angle[B, 1:T])      ## phase angle of the bus i
    @variable(Q, P[L, 1:T] >= 0)       ## real power flow on line l; elements in L is Tuple (i, j)
    @variable(Q, s[G, 1:T] >= 0)       ## real power generation at generator g
    @variable(Q, 0 <= x[D, 1:T] <= 1)  ## load shedding


    @variable(Q, zg[G, 1:T], Bin)      ## binary status indicator of generator g
    @variable(Q, zb[B, 1:T], Bin)      ## binary status indicator of bus i
    @variable(Q, zl[L, 1:T], Bin)      ## binary status indicator of line l

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

    ## constraint 1g h i j
    for i in B 
      @constraint(Q, [t in 1:T, d in Dᵢ[i]], zb[i, t] >= x[d, t] )
      @constraint(Q, [t in 1:T, g in Gᵢ[i]], zb[i, t] >= zg[g, t])
      @constraint(Q, [t in 1:T, j in out_L[i]], zb[i, t] >= zl[(i, j), t] )
      @constraint(Q, [t in 1:T, j in in_L[i]], zb[i, t] >= zl[(j, i), t] )
    end


    ## constraint 1k l m
    @constraint(Q, [i in B, t in 1:T-1], zb[i, t] >= zb[i, t+1] )
    @constraint(Q, [g in G, t in 1:T-1], zg[g, t] >= zg[g, t+1] )
    @constraint(Q, [l in L, t in 1:T-1], zl[l, t] >= zl[l, t+1] )


    # add cut for value functions & objective function
    @variable(Q, θ[Ω] >= θ_bound)      ## sddp, estimated value function

    # # add cut for value functions
    # for ω in Ω
    #   cut_coefficient = cut_collection[ω]
    #   iter = length(keys(cut_coefficient.v))  ## iter num
    #   k = length(keys(cut_coefficient.v[1]))  ## scenario num
    #   τ = Ω_rv[ω].τ
    #   @constraint(Q, [i in 1:iter-1, m in 1:k], θ[ω] .>= cut_coefficient.v[i][m] + 
    #                                             cut_coefficient.πb[i][m]' * zb[:, τ-1] +
    #                                             cut_coefficient.πg[i][m]' * zg[:, τ-1] +
    #                                             cut_coefficient.πl[i][m]' * zl[:, τ-1]
    #                                             )
    # end
    
    ## objective function
    @objective(Q, Min, sum( prob[ω] * ( sum( sum(paramDemand.w[d] * paramDemand.demand[t][d] * (1 - x[d, t]) for d in D )
                                                                                  for t in 1:Ω_rv[ω].τ - 1 ) + θ[ω] )
                                                                                    for ω in Ω)  
                                                                                )

    forwardInfo = ForwardInfo(Q, θ, zg, zb, zl)                                                                       
    ####################################################### solve the model and display the result ###########################################################  
    # optimize!(Q)

    # state_variable = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 2}}(:zg => round.(JuMP.value.(zg)), :zb => round.(JuMP.value.(zb)), :zl => round.(JuMP.value.(zl)))
    # state_value    = JuMP.objective_value(Q) - sum(prob[ω] * JuMP.value(θ[ω]) for ω in Ω)       ## 1a first term

    # return (state_variable = state_variable, 
    #           state_value = state_value, 
    #           obj_value = JuMP.objective_value(Q))  ## returen [state_variable, first_stage value, objective_value(Q)]

    return forwardInfo
end




"""     forward_stage2_optimize!
1. Solve the forward problem Stage 2
2. return the decision value -- y and its state value

"""

function forward_stage2_model!(indexSets::IndexSets, 
                                    paramDemand::ParamDemand, 
                                    paramOPF::ParamOPF, 
                                    randomVariables::RandomVariables                          ## realization of the random time
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
    # @constraint(Q, [i in B], yb[i] <= ẑ[:zb][i] )
    # @constraint(Q, [g in G], yg[g] <= ẑ[:zg][g] )
    # @constraint(Q, [l in L], yl[l] <= ẑ[:zl][l] )

    @constraint(Q, [i in B], yb[i] <= 1- νb[i] )
    @constraint(Q, [g in G], yg[g] <= 1- νg[g] )
    @constraint(Q, [l in L], yl[l] <= 1- νl[l] )

    @constraint(Q, [i in B], νb[i] >= randomVariables.vb[i] )
    @constraint(Q, [g in G], νg[g] >= randomVariables.vg[g] )
    @constraint(Q, [l in L], νl[l] >= randomVariables.vl[l] )

   
    for i in B 
      ## constraint 3e
      @constraint(Q, [t in randomVariables.τ:T], sum(s[g, t] for g in Gᵢ[i]) + sum(P[(i, j), t] for j in out_L[i] ) .== sum(paramDemand.demand[t][d] * x[d, t] for d in Dᵢ[i]) )

      ## constraint g h i j
      @constraint(Q, [t in randomVariables.τ:T, d in Dᵢ[i]], yb[i] >= x[d, t] )
      @constraint(Q, [g in Gᵢ[i]], yb[i] >= yg[g])
      @constraint(Q, [j in out_L[i]], yb[i] >= yl[(i, j)] )
      @constraint(Q, [j in in_L[i]], yb[i] >= yl[(j, i)] )

      ## constraint n
      # @constraint(Q, [j in unique(randomVariables.Ibb[i])], νb[j] >= randomVariables.ub[i] * ẑ[:zb][i] )
      # @constraint(Q, [j in unique(randomVariables.Ibg[i])], νg[j] >= randomVariables.ub[i] * ẑ[:zb][i] )
      # @constraint(Q, [j in unique(randomVariables.Ibl[i])], νl[j] >= randomVariables.ub[i] * ẑ[:zb][i] )
    end

    for g in G 
      ## constraint n
      # @constraint(Q, [j in unique(randomVariables.Igb[g])], νb[j] >= randomVariables.ug[g] * ẑ[:zg][g] )
      # @constraint(Q, [j in unique(randomVariables.Igg[g])], νg[j] >= randomVariables.ug[g] * ẑ[:zg][g] )
      # @constraint(Q, [j in unique(randomVariables.Igl[g])], νl[j] >= randomVariables.ug[g] * ẑ[:zg][g] )

       ## constraint 3f
      @constraint(Q, [t in randomVariables.τ:T], s[g, t] >= paramOPF.smin[g] * yg[g])
      @constraint(Q, [t in randomVariables.τ:T], s[g, t] <= paramOPF.smax[g] * yg[g])
    end

    for l in L 
      ## constraint n
      # @constraint(Q, [j in unique(randomVariables.Ilb[l])], νb[j] >= randomVariables.ul[l] * ẑ[:zl][l] )
      # @constraint(Q, [j in unique(randomVariables.Ilg[l])], νg[j] >= randomVariables.ul[l] * ẑ[:zl][l] )
      # @constraint(Q, [j in unique(randomVariables.Ill[l])], νl[j] >= randomVariables.ul[l] * ẑ[:zl][l] )


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
            sum( sum(paramDemand.w[d] * paramDemand.demand[t][d] * (1 - x[d, t]) for d in D ) for t in randomVariables.τ:T) +
            sum(paramDemand.cb[i] * νb[i] for i in B) + 
            sum(paramDemand.cg[g] * νg[g] for g in G) + 
            sum(paramDemand.cl[l] * νl[l] for l in L) + paramDemand.penalty * slack_variable_b - paramDemand.penalty * slack_variable_c
            )
    ####################################################### solve the model and display the result ###########################################################
    # optimize!(Q)
    # state_variable = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}(:zg => round.(JuMP.value.(yg)), :zb => round.(JuMP.value.(yb)), :zl => round.(JuMP.value.(yl)))
    # state_value    = JuMP.objective_value(Q)

    # return [state_variable, state_value]  ## returen [Lt, y, θ, f]
    forward2Info = Forward2Info(Q, yb, yg, yl, νb, νg, νl)

    return forward2Info

end





function forward_stage2_modify_constraints!(indexSets::IndexSets, 
                                    forward2Info::Forward2Info,
                                    iter::Int64,
                                    ẑ::Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}},
                                    randomVariables::RandomVariables                          ## realization of the random time
                                    )

    ## modify the constraints according to the first stage state variables
    if iter > 1
        for i in indexSets.B
            delete(forward2Info.model, forward2Info.model[:k3][i])
            for j in unique(randomVariables.Ibb[i])
                delete(forward2Info.model, forward2Info.model[:n3B1][i, j])
            end
            for j in unique(randomVariables.Ibg[i])
                delete(forward2Info.model, forward2Info.model[:n3B2][i, j])
            end
            for j in unique(randomVariables.Ibl[i])
                delete(forward2Info.model, forward2Info.model[:n3B3][i, j])
            end
        end
        unregister(forward2Info.model, :k3)
        unregister(forward2Info.model, :n3B1)
        unregister(forward2Info.model, :n3B2)
        unregister(forward2Info.model, :n3B3)


        for g in indexSets.G
            delete(forward2Info.model, forward2Info.model[:l3][g])
            for j in unique(randomVariables.Igb[g])
                delete(forward2Info.model, forward2Info.model[:n3G1][g, j])
            end
            for j in unique(randomVariables.Igg[g])
                delete(forward2Info.model, forward2Info.model[:n3G2][g, j])
            end
            for j in unique(randomVariables.Igl[g])
                delete(forward2Info.model, forward2Info.model[:n3G3][g, j])
            end
        end
        unregister(forward2Info.model, :l3)
        unregister(forward2Info.model, :n3G1)
        unregister(forward2Info.model, :n3G2)
        unregister(forward2Info.model, :n3G3)

        for l in indexSets.L
            delete(forward2Info.model, forward2Info.model[:m3][l])
            for j in unique(randomVariables.Ilb[l])
              delete(forward2Info.model, forward2Info.model[:n3L1][l, j])
            end
            for j in unique(randomVariables.Ilg[l])
              delete(forward2Info.model, forward2Info.model[:n3L2][l, j])
            end
            for j in unique(randomVariables.Ill[l])
              delete(forward2Info.model, forward2Info.model[:n3L3][l, j])
            end
        end
        unregister(forward2Info.model, :m3)
        unregister(forward2Info.model, :n3L1)
        unregister(forward2Info.model, :n3L2)
        unregister(forward2Info.model, :n3L3)
    end
    # constraint k l m 
    @constraint(forward2Info.model, k3[i in indexSets.B], forward2Info.yb[i] <= ẑ[:zb][i] )
    @constraint(forward2Info.model, l3[g in indexSets.G], forward2Info.yg[g] <= ẑ[:zg][g] )
    @constraint(forward2Info.model, m3[l in indexSets.L], forward2Info.yl[l] <= ẑ[:zl][l] )
    # constraint n
    @constraint(forward2Info.model, n3L1[l in indexSets.L, j in unique(randomVariables.Ilb[l])], forward2Info.νb[j] >= randomVariables.ul[l] * ẑ[:zl][l] )
    @constraint(forward2Info.model, n3L2[l in indexSets.L, j in unique(randomVariables.Ilg[l])], forward2Info.νg[j] >= randomVariables.ul[l] * ẑ[:zl][l] )
    @constraint(forward2Info.model, n3L3[l in indexSets.L, j in unique(randomVariables.Ill[l])], forward2Info.νl[j] >= randomVariables.ul[l] * ẑ[:zl][l] )

    @constraint(forward2Info.model, n3G1[g in indexSets.G, j in unique(randomVariables.Igb[g])], forward2Info.νb[j] >= randomVariables.ug[g] * ẑ[:zg][g] )
    @constraint(forward2Info.model, n3G2[g in indexSets.G, j in unique(randomVariables.Igg[g])], forward2Info.νg[j] >= randomVariables.ug[g] * ẑ[:zg][g] )
    @constraint(forward2Info.model, n3G3[g in indexSets.G, j in unique(randomVariables.Igl[g])], forward2Info.νl[j] >= randomVariables.ug[g] * ẑ[:zg][g] )

    @constraint(forward2Info.model, n3B1[i in indexSets.B, j in unique(randomVariables.Ibb[i])], forward2Info.νb[j] >= randomVariables.ub[i] * ẑ[:zb][i] )
    @constraint(forward2Info.model, n3B2[i in indexSets.B, j in unique(randomVariables.Ibg[i])], forward2Info.νg[j] >= randomVariables.ub[i] * ẑ[:zb][i] )
    @constraint(forward2Info.model, n3B3[i in indexSets.B, j in unique(randomVariables.Ibl[i])], forward2Info.νl[j] >= randomVariables.ub[i] * ẑ[:zb][i] )


end

