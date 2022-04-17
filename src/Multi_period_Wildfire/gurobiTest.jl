using JuMP, Gurobi, Random

const GRB_ENV = Gurobi.Env()


include("data_struct.jl")

## input data

# include("generationTest.jl")
# include("runtests_small2.jl")
# include("runtests_small3.jl")


########################################################################################################################################################
############################################  auxiliary function: nonanticipativity for multistage problem #############################################
########################################################################################################################################################
function recursion_scenario_constraint(pathList::Vector{Int64}, P::Float64, scenario_sequence::Dict{Int64, Dict{Int64, Any}}, t::Int64;   
                    Ω::Dict{Int64,Dict{Int64,RandomVariables}} = Ω, prob::Dict{Int64,Vector{Float64}} = prob, T::Int64 = 2, gurobiModelInfo::GurobiModelInfo = gurobiModelInfo)

    if t <= T
        for ω_key in keys(Ω[t])

            pathList_copy = copy(pathList)
            P_copy = copy(P)

            push!(pathList_copy, ω_key)
            P_copy = P_copy * prob[t][ω_key]

            ## nonanticipativity for multi-stage problem
            if t < T
                if haskey(scenario_sequence, 1)
                    first =  maximum(keys(scenario_sequence)) + 1
                    last  =  maximum(keys(scenario_sequence)) + gurobiModelInfo.num_Ω^(T-t)
                    @constraint(gurobiModelInfo.model, [i = first, j in (first + 1): last], gurobiModelInfo.x[:, t, i] .== gurobiModelInfo.x[:, t, j]) 
                    @constraint(gurobiModelInfo.model, [i = first, j in (first + 1): last], gurobiModelInfo.y[:, t, i] .== gurobiModelInfo.y[:, t, j]) 
                    @constraint(gurobiModelInfo.model, [i = first, j in (first + 1): last], gurobiModelInfo.slack[t, i] .== gurobiModelInfo.slack[t, j]) 

                else
                    @constraint(gurobiModelInfo.model, [i = 1, j in 2:gurobiModelInfo.num_Ω^(T-t)], gurobiModelInfo.x[:, t, i] .== gurobiModelInfo.x[:, t, j]) 
                    @constraint(gurobiModelInfo.model, [i = 1, j in 2:gurobiModelInfo.num_Ω^(T-t)], gurobiModelInfo.y[:, t, i] .== gurobiModelInfo.y[:, t, j]) 
                    @constraint(gurobiModelInfo.model, [i = 1, j in 2:gurobiModelInfo.num_Ω^(T-t)], gurobiModelInfo.slack[t, i] .== gurobiModelInfo.slack[t, j]) 
                end
            end
            recursion_scenario_constraint(pathList_copy, P_copy, scenario_sequence, t+1, Ω = Ω, prob = prob, T = T, gurobiModelInfo = gurobiModelInfo)
        end
    else
        if haskey(scenario_sequence, 1)
            scenario_sequence[maximum(keys(scenario_sequence))+1] = Dict(1 => pathList, 2 => P)
        else
            scenario_sequence[1] = Dict(1 => pathList, 2 => P)
        end
        return scenario_sequence
    end

end







################################################################################################################################################
############################################################     Gurobi function   #############################################################
################################################################################################################################################

function gurobiOptimize!(Ω::Dict{Int64,Dict{Int64,RandomVariables}}, prob::Dict{Int64,Vector{Float64}}, StageCoefficient::Dict{Int64, StageData}; 
                        binaryInfo::BinaryInfo = binaryInfo)

    (A, n, d) = (binaryInfo.A, binaryInfo.n, binaryInfo.d)
    T = length(Ω); num_Ω = length(Ω[1]);
    W = num_Ω^(T-1) # number of scenarios


    model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                                "OutputFlag" => 1, 
                                                "Threads" => 1) 
                                                )
    @variable(model, x[i = 1:d, t = 1:T, ω in 1:W] >= 0, Int)   ## for current state, x is the number of generators will be built in this stage
    @variable(model, y[i = 1:d, t = 1:T, ω in 1:W] >= 0)        ## amount of electricity
    @variable(model, slack[t = 1:T, ω in 1:W] >=0 )
    @variable(model, S[i = 1:d, t = 1:T, ω in 1:W] >=0 )

    @constraint(model, [t in 1:T, ω in 1:W], S[:, t, ω] .== sum(x[:, j, ω] for j in 1:t ) )
    @constraint(model, [t in 1:T], S[:, t, :] .<= StageCoefficient[t].ū)
    @constraint(model, [t in 1:T, ω in 1:W], y[:,t, ω] .<= StageCoefficient[t].h * StageCoefficient[t].N * (S[:, t, ω] + StageCoefficient[t].s₀))

    ## nonanticipativity for multistage problem 
    @constraint(model, [ω in 2:W], x[:, 1, 1] .== x[:, 1, ω])  ## nonanticipativity for 2-stage problem
    gurobiModelInfo = GurobiModelInfo(model, x, y, slack, num_Ω)
    scenario_sequence = Dict{Int64, Dict{Int64, Any}}()  ## the first index is for scenario index, the second one is for stage
    pathList = Vector{Int64}()
    push!(pathList, 1)

    recursion_scenario_constraint(pathList, 1.0, scenario_sequence, 2, Ω = Ω, prob = prob, T = T, gurobiModelInfo = gurobiModelInfo)
    
    scenario_tree = scenario_sequence
    #############################################################################################################################
    @constraint(model, [t in 1:T, ω in 1:W], sum(y[i, t, ω] for i in 1:d) + slack[t, ω] >= Ω[t][scenario_tree[ω][1][t]].d[1] )

    @objective(model, Min, sum( sum( scenario_tree[ω][2] * (StageCoefficient[t].c1' * x[:, t, ω] + StageCoefficient[t].c2' * y[:, t, ω] + StageCoefficient[t].penalty * slack[t, ω]) for t in 1:T ) for ω in 1:W) )
    optimize!(model)
    ####################################################### solve the model and display the result ###########################################################

    # JuMP.objective_value(model)
    # JuMP.value.(x[:, 1, :])
    # JuMP.value.(x[:, 2, :])
    # JuMP.value.(x[:, 3, :])

    # JuMP.value.(y[:, 2, :])

    # JuMP.value.(slack)
    gurobiResult = Dict(1=> JuMP.objective_value(model), 2 =>  JuMP.value.(x[:, 1, 1]), 3=>  JuMP.value.(y[:, 1, 1]))

    return gurobiResult
end



# gurobiOptimize!(Ω, prob, StageCoefficient,
#                         binaryInfo = binaryInfo)








function gurobiOptimize!(indexSets::IndexSets, 
                            paramDemand::ParamDemand, 
                            paramOPF::ParamOPF, 
                            Ω_rv::Dict{Int64, RandomVariables},
                            prob::Dict{Int64, Float64})  


    (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, indexSets.Ω) 
    (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L) 

    model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                          "OutputFlag" => 0, 
                                          "Threads" => 1) 
                                          )
    
    ## the first stage variables
    @variable(model, θ_angle[B, 1:T])      ## phase angle of the bus i
    @variable(model, P[L, 1:T] >= 0)       ## real power flow on line l; elements in L is Tuple (i, j)
    @variable(model, s[G, 1:T] >= 0)       ## real power generation at generator g
    @variable(model, 0 <= x[D, 1:T] <= 1)  ## load shedding


    @variable(model, zg[G, 1:T], Bin)      ## binary status indicator of generator g
    @variable(model, zb[B, 1:T], Bin)      ## binary status indicator of bus i
    @variable(model, zl[L, 1:T], Bin)      ## binary status indicator of line l

    ## the second stage variables
    @variable(model, θω[B, 1:T, Ω])             ## phase angle of the bus i
    @variable(model, Pω[L, 1:T, Ω] >= 0)        ## real power flow on line l; elements in L is Tuple (i, j)  
    @variable(model, sω[G, 1:T, Ω] >= 0)        ## real power generation at generator g
    @variable(model, 0 <= xω[D, 1:T, Ω] <= 1)   ## load shedding

    @variable(model, yb[B, Ω], Bin)
    @variable(model, yg[G, Ω], Bin)
    @variable(model, yl[L, Ω], Bin)

    @variable(model, νb[B, Ω], Bin)
    @variable(model, νg[G, Ω], Bin)
    @variable(model, νl[L, Ω], Bin)

    @variable(model, slack_variable_b[Ω] >= 0)
    @variable(model, slack_variable_c[Ω] <= 0)

    # constraint 1b 1c
    for l in L
      i = l[1]
      j = l[2]
      @constraint(model, [t in 1:T], P[l, t] <= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmax * (1 - zl[l, t] ) ) )
      @constraint(model, [t in 1:T], P[l, t] >= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmin * (1 - zl[l, t] ) ) )

      ## constraints 3b 3c
      @constraint(model, [ω in Ω, t in Ω_rv[ω].τ:T], Pω[l, t, ω] <= - paramOPF.b[l] * (θω[i, t, ω] - θω[j, t, ω] + paramOPF.θmax * (1 - yl[l, ω] ) ) + slack_variable_b[ω] )
      @constraint(model, [ω in Ω, t in Ω_rv[ω].τ:T], Pω[l, t, ω] >= - paramOPF.b[l] * (θω[i, t, ω] - θω[j, t, ω] + paramOPF.θmin * (1 - yl[l, ω] ) ) + slack_variable_c[ω] )
    end

    ## constraint 1d
    @constraint(model, [l in L, t in 1:T], P[l, t] >= - paramOPF.W[l] * zl[l, t] )
    @constraint(model, [l in L, t in 1:T], P[l, t] <= paramOPF.W[l] * zl[l, t] )

    ## constraint 1e
    @constraint(model, [i in B, t in 1:T], sum(s[g, t] for g in Gᵢ[i]) + sum(P[(i, j), t] for j in out_L[i] ) .== sum(paramDemand.demand[t][d] * x[d, t] for d in Dᵢ[i]) )
    
    ## constraint 1f
    @constraint(model, [g in G, t in 1:T], s[g, t] >= paramOPF.smin[g] * zg[g, t] )
    @constraint(model, [g in G, t in 1:T], s[g, t] <= paramOPF.smax[g] * zg[g, t] )

    ## constraint 1g h i j
    @constraint(model, [i in B, t in 1:T, d in Dᵢ[i]], zb[i, t] >= x[d, t] )
    @constraint(model, [i in B, t in 1:T, g in Gᵢ[i]], zb[i, t] >= zg[g, t])
    @constraint(model, [i in B, t in 1:T, j in out_L[i]], zb[i, t] >= zl[(i, j), t] )
    @constraint(model, [i in B, t in 1:T, j in in_L[i]], zb[i, t] >= zl[(j, i), t] )


    ## constraint 1k l m
    @constraint(model, [i in B, t in 1:T-1], zb[i, t] >= zb[i, t+1] )
    @constraint(model, [g in G, t in 1:T-1], zg[g, t] >= zg[g, t+1] )
    @constraint(model, [l in L, t in 1:T-1], zl[l, t] >= zl[l, t+1] )



    ## constraint 3d
    @constraint(model, [ω in Ω, l in L, t in Ω_rv[ω].τ:T], Pω[l, t, ω] >= - paramOPF.W[l] * yl[l, ω] )
    @constraint(model, [ω in Ω, l in L, t in Ω_rv[ω].τ:T], Pω[l, t, ω] <= paramOPF.W[l] * yl[l, ω] )

    ## constraint 3e
    @constraint(model, [ω in Ω, i in B, t in Ω_rv[ω].τ:T], sum(sω[g, t, ω] for g in Gᵢ[i]) + sum(Pω[(i, j), t, ω] for j in out_L[i] ) .== sum(paramDemand.demand[t][d] * xω[d, t, ω] for d in Dᵢ[i]) )

    ## constraint 3f
    @constraint(model, [ω in Ω, g in G, t in Ω_rv[ω].τ:T], sω[g, t, ω] >= paramOPF.smin[g] * yg[g, ω])
    @constraint(model, [ω in Ω, g in G, t in Ω_rv[ω].τ:T], sω[g, t, ω] <= paramOPF.smax[g] * yg[g, ω])

    ## constraint 3g h i j
    @constraint(model, [ω in Ω, i in B, t in Ω_rv[ω].τ:T, d in Dᵢ[i]], yb[i, ω] >= xω[d, t, ω] )
    @constraint(model, [ω in Ω, i in B, g in Gᵢ[i]], yb[i, ω] >= yg[g, ω])
    @constraint(model, [ω in Ω, i in B, j in out_L[i]], yb[i, ω] >= yl[(i, j), ω] )
    @constraint(model, [ω in Ω, i in B, j in in_L[i]], yb[i, ω] >= yl[(j, i), ω] )

    ## constraint 3k l m 
    @constraint(model, [ω in Ω, i in B], yb[i, ω] <= zb[i, Ω_rv[ω].τ-1] ) 
    @constraint(model, [ω in Ω, g in G], yg[g, ω] <= zg[g, Ω_rv[ω].τ-1] ) 
    @constraint(model, [ω in Ω, l in L], yl[l, ω] <= zl[l, Ω_rv[ω].τ-1] )

    @constraint(model, [ω in Ω, i in B], yb[i, ω] <= 1- νb[i, ω] )
    @constraint(model, [ω in Ω, g in G], yg[g, ω] <= 1- νg[g, ω] )
    @constraint(model, [ω in Ω, l in L], yl[l, ω] <= 1- νl[l, ω] )

    @constraint(model, [ω in Ω, i in B], νb[i, ω] >= Ω_rv[ω].vb[i] )
    @constraint(model, [ω in Ω, g in G], νg[g, ω] >= Ω_rv[ω].vg[g] )
    @constraint(model, [ω in Ω, l in L], νl[l, ω] >= Ω_rv[ω].vl[l] )

    ## constraint 3n
    @constraint(model, [ω in Ω, i in B, j in unique(Ω_rv[ω].Ibb[i])], νb[j, ω] >= Ω_rv[ω].ub[i] * zb[i, Ω_rv[ω].τ-1] )
    @constraint(model, [ω in Ω, i in B, j in unique(Ω_rv[ω].Ibg[i])], νg[j, ω] >= Ω_rv[ω].ub[i] * zb[i, Ω_rv[ω].τ-1] )
    @constraint(model, [ω in Ω, i in B, j in unique(Ω_rv[ω].Ibl[i])], νl[j, ω] >= Ω_rv[ω].ub[i] * zb[i, Ω_rv[ω].τ-1] )

    @constraint(model, [ω in Ω, i in G, j in unique(Ω_rv[ω].Igb[i])], νb[j, ω] >= Ω_rv[ω].ug[i] * zg[i, Ω_rv[ω].τ-1] )
    @constraint(model, [ω in Ω, i in G, j in unique(Ω_rv[ω].Igg[i])], νg[j, ω] >= Ω_rv[ω].ug[i] * zg[i, Ω_rv[ω].τ-1] )
    @constraint(model, [ω in Ω, i in G, j in unique(Ω_rv[ω].Igl[i])], νl[j, ω] >= Ω_rv[ω].ug[i] * zg[i, Ω_rv[ω].τ-1] )

    @constraint(model, [ω in Ω, i in L, j in unique(Ω_rv[ω].Ilb[i])], νb[j, ω] >= Ω_rv[ω].ul[i] * zl[i, Ω_rv[ω].τ-1] )
    @constraint(model, [ω in Ω, i in L, j in unique(Ω_rv[ω].Ilg[i])], νg[j, ω] >= Ω_rv[ω].ul[i] * zl[i, Ω_rv[ω].τ-1] )
    @constraint(model, [ω in Ω, i in L, j in unique(Ω_rv[ω].Ill[i])], νl[j, ω] >= Ω_rv[ω].ul[i] * zl[i, Ω_rv[ω].τ-1] )

    
    ## objective function 1a & 3a
    @objective(model, Min, sum( prob[ω] * ( sum( sum(paramDemand.w[d] * paramDemand.demand[t][d] * (1 - x[d, t]) for d in D)
                for t in 1:Ω_rv[ω].τ - 1 ) + ## 3a
                sum( sum(paramDemand.w[d] * paramDemand.demand[t][d] * (1 - xω[d, t, ω]) for d in D) for t in Ω_rv[ω].τ:T) + 
                sum(paramDemand.cb[i] * νb[i, ω] for i in B) + 
                sum(paramDemand.cg[g] * νg[g, ω] for g in G) + 
                sum(paramDemand.cl[l] * νl[l, ω] for l in L) + 
                + paramDemand.penalty * slack_variable_b[ω] - paramDemand.penalty * slack_variable_c[ω]) 
                for ω in Ω)  
                    )

    @objective(Q, Min,  sum( sum(paramDemand.w[d] * paramDemand.demand[t][d] * (1 - x[d, t]) for d in D) for t in randomVariables.τ:T) +
                        sum(paramDemand.cb[i] * νb[i] for i in B) + 
                        sum(paramDemand.cg[g] * νg[g] for g in G) + 
                        sum(paramDemand.cl[l] * νl[l] for l in L) + paramDemand.penalty * slack_variable_b - paramDemand.penalty * slack_variable_c
                        )
                                                                                     
    ####################################################### solve the model and display the result ###########################################################
    optimize!(model) 
    first_state_variable = Dict(:zg => round.(JuMP.value.(zg)), 
                                    :zb => round.(JuMP.value.(zb)), 
                                    :zl => round.(JuMP.value.(zl))
                                    )
    second_state_variable = Dict(:zg => round.(JuMP.value.(yg)), 
                                    :zb => round.(JuMP.value.(yb)), 
                                    :zl => round.(JuMP.value.(yl))
                                    )

    gurobiResult = Dict("OPT"=> JuMP.objective_value(model), "1st"=>  first_state_variable, "2nd"=>  second_state_variable )

    return gurobiResult
end














