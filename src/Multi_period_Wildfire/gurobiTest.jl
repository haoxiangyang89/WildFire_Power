using JuMP, Gurobi, Random

const GRB_ENV = Gurobi.Env()


include("data_struct.jl")

## input data

# include("generationTest.jl")
# include("runtests_small2.jl")
# include("runtests_small3.jl")

################################################################################################################################################
############################################################     Gurobi function   #############################################################
################################################################################################################################################

function gurobiOptimize!(indexSets::IndexSets, 
                            paramDemand::ParamDemand, 
                            paramOPF::ParamOPF, 
                            Ω_rv::Dict{Int64, RandomVariables},
                            prob::Dict{Int64, Float64})  


    (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, indexSets.Ω) 
    (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L) 

    model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                          "OutputFlag" => 0, 
                                          "Threads" =>0) 
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

    ## second stage constraints
    for ω in Ω 
        ## constraint 3d
        @constraint(model, [l in L, t in Ω_rv[ω].τ:T], Pω[l, t, ω] >= - paramOPF.W[l] * yl[l, ω] )
        @constraint(model, [l in L, t in Ω_rv[ω].τ:T], Pω[l, t, ω] <= paramOPF.W[l] * yl[l, ω] )

        ## constraint 3e
        @constraint(model, [i in B, t in Ω_rv[ω].τ:T], sum(sω[g, t, ω] for g in Gᵢ[i]) + sum(Pω[(i, j), t, ω] for j in out_L[i] ) .== sum(paramDemand.demand[t][d] * xω[d, t, ω] for d in Dᵢ[i]) )

        ## constraint 3f
        @constraint(model, [g in G, t in Ω_rv[ω].τ:T], sω[g, t, ω] >= paramOPF.smin[g] * yg[g, ω])
        @constraint(model, [g in G, t in Ω_rv[ω].τ:T], sω[g, t, ω] <= paramOPF.smax[g] * yg[g, ω])

        ## constraint 3g h i j
        @constraint(model, [i in B, t in Ω_rv[ω].τ:T, d in Dᵢ[i]], yb[i, ω] >= xω[d, t, ω] )
        @constraint(model, [i in B, g in Gᵢ[i]], yb[i, ω] >= yg[g, ω])
        @constraint(model, [i in B, j in out_L[i]], yb[i, ω] >= yl[(i, j), ω] )
        @constraint(model, [i in B, j in in_L[i]], yb[i, ω] >= yl[(j, i), ω] )

        ## constraint 3k l m 
        @constraint(model, [i in B], yb[i, ω] <= zb[i, Ω_rv[ω].τ-1] ) 
        @constraint(model, [g in G], yg[g, ω] <= zg[g, Ω_rv[ω].τ-1] ) 
        @constraint(model, [l in L], yl[l, ω] <= zl[l, Ω_rv[ω].τ-1] )

        @constraint(model, [i in B], yb[i, ω] <= 1- νb[i, ω] )
        @constraint(model, [g in G], yg[g, ω] <= 1- νg[g, ω] )
        @constraint(model, [l in L], yl[l, ω] <= 1- νl[l, ω] )

        @constraint(model, [i in B], νb[i, ω] >= Ω_rv[ω].vb[i] )
        @constraint(model, [g in G], νg[g, ω] >= Ω_rv[ω].vg[g] )
        @constraint(model, [l in L], νl[l, ω] >= Ω_rv[ω].vl[l] )

        ## constraint 3n
        @constraint(model, [i in B, j in unique(Ω_rv[ω].Ibb[i])], νb[j, ω] >= Ω_rv[ω].ub[i] * zb[i, Ω_rv[ω].τ-1] )
        @constraint(model, [i in B, j in unique(Ω_rv[ω].Ibg[i])], νg[j, ω] >= Ω_rv[ω].ub[i] * zb[i, Ω_rv[ω].τ-1] )
        @constraint(model, [i in B, j in unique(Ω_rv[ω].Ibl[i])], νl[j, ω] >= Ω_rv[ω].ub[i] * zb[i, Ω_rv[ω].τ-1] )

        @constraint(model, [i in G, j in unique(Ω_rv[ω].Igb[i])], νb[j, ω] >= Ω_rv[ω].ug[i] * zg[i, Ω_rv[ω].τ-1] )
        @constraint(model, [i in G, j in unique(Ω_rv[ω].Igg[i])], νg[j, ω] >= Ω_rv[ω].ug[i] * zg[i, Ω_rv[ω].τ-1] )
        @constraint(model, [i in G, j in unique(Ω_rv[ω].Igl[i])], νl[j, ω] >= Ω_rv[ω].ug[i] * zg[i, Ω_rv[ω].τ-1] )

        @constraint(model, [i in L, j in unique(Ω_rv[ω].Ilb[i])], νb[j, ω] >= Ω_rv[ω].ul[i] * zl[i, Ω_rv[ω].τ-1] )
        @constraint(model, [i in L, j in unique(Ω_rv[ω].Ilg[i])], νg[j, ω] >= Ω_rv[ω].ul[i] * zl[i, Ω_rv[ω].τ-1] )
        @constraint(model, [i in L, j in unique(Ω_rv[ω].Ill[i])], νl[j, ω] >= Ω_rv[ω].ul[i] * zl[i, Ω_rv[ω].τ-1] )
    end 
    
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


    return (OPT = JuMP.objective_value(model), 
            first_state_variable = first_state_variable, 
            second_state_variable = second_state_variable)
end














