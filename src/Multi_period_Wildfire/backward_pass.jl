#############################################################################################
###################################  function: backward pass ################################
#############################################################################################

"""
    This is the oracle in level set method, and it will return [F, ∇F]
    where ∇F is a Dict{Symbol, Vector}

"""

function backward_stage2_optimize!(indexSets::IndexSets, 
                                    paramDemand::ParamDemand, 
                                    paramOPF::ParamOPF, 
                                    ẑ::Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}},
                                    randomVariables::RandomVariables                          ## realization of the random time
                                    # π::Dict{Symbol, Vector{Float64}}
                                    )

    # (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, indexSets.Ω)
    # (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L) 

    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                "OutputFlag" => 0, 
                "Threads" => 0) 
                )

    @variable(Q, θ_angle[indexSets.B, 1:indexSets.T]) 
    @variable(Q, P[indexSets.L, 1:indexSets.T] >= 0) ## elements in L is Tuple (i, j)
    @variable(Q, s[indexSets.G, 1:indexSets.T] >= 0)
    @variable(Q, 0 <= x[indexSets.D, 1:indexSets.T] <= 1)

    @variable(Q, yb[indexSets.B], Bin)
    @variable(Q, yg[indexSets.G], Bin)
    @variable(Q, yl[indexSets.L], Bin)

    @variable(Q, νb[indexSets.B], Bin)
    @variable(Q, νg[indexSets.G], Bin)
    @variable(Q, νl[indexSets.L], Bin)

    @variable(Q, 0 <= zg[indexSets.G] <= 1)
    @variable(Q, 0 <= zb[indexSets.B] <= 1)
    @variable(Q, 0 <= zl[indexSets.L] <= 1)

    @variable(Q, slack_variable_b >= 0)
    @variable(Q, slack_variable_c <= 0)


    ## constraint 3e
    @constraint(Q, [i in indexSets.B, t in randomVariables.τ:indexSets.T], sum(s[g, t] for g in indexSets.Gᵢ[i]) + sum(P[(i, j), t] for j in indexSets.out_L[i] ) .== sum(paramDemand.demand[t][d] * x[d, t] for d in indexSets.Dᵢ[i]) )

    ## constraint k l m 
    @constraint(Q, [i in indexSets.B], yb[i] <= zb[i] )
    @constraint(Q, [g in indexSets.G], yg[g] <= zg[g] )
    @constraint(Q, [l in indexSets.L], yl[l] <= zl[l] )

    @constraint(Q, [i in indexSets.B], yb[i] <= 1- νb[i] )
    @constraint(Q, [g in indexSets.G], yg[g] <= 1- νg[g] )
    @constraint(Q, [l in indexSets.L], yl[l] <= 1- νl[l] )

    @constraint(Q, [i in indexSets.B], νb[i] >= randomVariables.vb[i] )
    @constraint(Q, [g in indexSets.G], νg[g] >= randomVariables.vg[g] )
    @constraint(Q, [l in indexSets.L], νl[l] >= randomVariables.vl[l] )

    
    for i in indexSets.B 
        ## constraint n
        @constraint(Q, [j in unique(randomVariables.Ibb[i])], νb[j] >= randomVariables.ub[i] * zb[i] )
        @constraint(Q, [j in unique(randomVariables.Ibg[i])], νg[j] >= randomVariables.ub[i] * zb[i] )
        @constraint(Q, [j in unique(randomVariables.Ibl[i])], νl[j] >= randomVariables.ub[i] * zb[i] )

         ## constraint g h i j
        @constraint(Q, [t in randomVariables.τ:indexSets.T, d in indexSets.Dᵢ[i]], yb[i] >= x[d, t] )
        @constraint(Q, [g in indexSets.Gᵢ[i]], yb[i] >= yg[g])
        @constraint(Q, [j in indexSets.out_L[i]], yb[i] >= yl[(i, j)] )
        @constraint(Q, [j in indexSets.in_L[i]], yb[i] >= yl[(j, i)] )
    end

    for g in indexSets.G
        ## constraint n 
        @constraint(Q, [j in unique(randomVariables.Igb[g])], νb[j] >= randomVariables.ug[g] * zg[g] )
        @constraint(Q, [j in unique(randomVariables.Igg[g])], νg[j] >= randomVariables.ug[g] * zg[g] )
        @constraint(Q, [j in unique(randomVariables.Igl[g])], νl[j] >= randomVariables.ug[g] * zg[g] )

        ## constraint 3f
        @constraint(Q, [g in indexSets.G, t in randomVariables.τ:indexSets.T], s[g, t] >= paramOPF.smin[g] * yg[g])
        @constraint(Q, [g in indexSets.G, t in randomVariables.τ:indexSets.T], s[g, t] <= paramOPF.smax[g] * yg[g])
    end

    for l in indexSets.L 
        ## constraint n
        @constraint(Q, [j in unique(randomVariables.Ilb[l])], νb[j] >= randomVariables.ul[l] * zl[l] )
        @constraint(Q, [j in unique(randomVariables.Ilg[l])], νg[j] >= randomVariables.ul[l] * zl[l] )
        @constraint(Q, [j in unique(randomVariables.Ill[l])], νl[j] >= randomVariables.ul[l] * zl[l] )

        ## constraint 3b 3c
        i = l[1]
        j = l[2]
        @constraint(Q, [t in randomVariables.τ:indexSets.T], P[l, t] <= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmax * (1 - yl[l] ) ) + slack_variable_b )
        @constraint(Q, [t in randomVariables.τ:indexSets.T], P[l, t] >= - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmin * (1 - yl[l] ) ) + slack_variable_c )

        ## constraint 3d
        @constraint(Q, [t in randomVariables.τ:indexSets.T], P[l, t] >= - paramOPF.W[l] * yl[l] )
        @constraint(Q, [t in randomVariables.τ:indexSets.T], P[l, t] <= paramOPF.W[l] * yl[l] )
    end

    ## objective function
    # @objective(Q, Min,  
    #         sum( sum(paramDemand.w[d] * paramDemand.demand[t][d] * (1 - x[d, t]) for d in indexSets.D ) for t in randomVariables.τ:indexSets.T) +
    #         sum(paramDemand.cb[i] * νb[i] for i in indexSets.B) + 
    #         sum(paramDemand.cg[g] * νg[g] for g in indexSets.G) + 
    #         sum(paramDemand.cl[l] * νl[l] for l in indexSets.L) +
    #         π[:zb]' * (ẑ[:zb] .- zb) + π[:zg]' * (ẑ[:zg] .- zg) + π[:zl]' * (ẑ[:zl] .- zl) 
    #         + paramDemand.penalty * slack_variable_b - paramDemand.penalty * slack_variable_c
    #         )


    # ####################################################### solve the model and display the result ###########################################################
    # optimize!(Q)
    # state_variable = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}(:zg => round.(JuMP.value.(zg)), 
    #                                                 :zb => round.(JuMP.value.(zb)), 
    #                                                 :zl => round.(JuMP.value.(zl))
    #                                                 )
    # F  = JuMP.objective_value(Q)
    # negative_∇F = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}(:zb => - ẑ[:zb] .+ state_variable[:zb],
    #                                                                     :zg => - ẑ[:zg] .+ state_variable[:zg],
    #                                                                     :zl => - ẑ[:zl] .+ state_variable[:zl]
    #                                                                     )

    backwardInfo = BackwardInfo(Q, x, νb, νg, νl, zg, zb, zl, slack_variable_b, slack_variable_c)
    return backwardInfo
end


#############################################################################################
###########################  auxiliary functions for level set method #######################
#############################################################################################

"""
    This function is to constraint the model for solving gap and alpha
"""

function Δ_model_formulation(function_info::FunctionInfo, f_star::Float64, iter::Int64; Output::Int64 = 0)
    
    model_alpha = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV),
            "OutputFlag" => Output, 
            "Threads" => 0)
            )

    @variable(model_alpha, z)
    @variable(model_alpha, 0 <= α <= 1)
    @constraint(model_alpha, con[j = 1:iter], z <=  α * ( function_info.f_his[j] - f_star) + (1 - α) * function_info.G_max_his[j] )
    
    # we first compute gap Δ
    @objective(model_alpha, Max, z)
    optimize!(model_alpha)
    st = termination_status(model_alpha)
    Δ = JuMP.objective_value(model_alpha)

    
    ## then we modify above model to compute alpha
    # alpha_min
    @constraint(model_alpha, z .>= 0)
    @objective(model_alpha, Min, α)
    optimize!(model_alpha)
    a_min = JuMP.value(α)

    # alpha_max
    @objective(model_alpha, Max, α)
    optimize!(model_alpha)
    a_max = JuMP.value(α)

    return Dict(1 => Δ, 2 => a_min, 3 => a_max)

end


"""
    This function is to add constraints for the model f_star and nxt pt.
"""

function add_constraint(function_info::FunctionInfo, model_info::ModelInfo, iter::Int64)
    m = length(function_info.G)

    xⱼ = function_info.x_his[iter]
    # add constraints     
    @constraint(model_info.model, model_info.z .>= function_info.f_his[iter] + function_info.df[:zb]' * (model_info.xb .- xⱼ[:zb]) 
                                                                    + function_info.df[:zg]' * (model_info.xg .- xⱼ[:zg] )
                                                                    + function_info.df[:zl]' * (model_info.xl .- xⱼ[:zl] )
                                                                    )

    @constraint(model_info.model, [k = 1:m], model_info.y .>= function_info.G[k] + sum(function_info.dG[k][:zb] .* (model_info.xb .- xⱼ[:zb]))
                                                                                 + sum(function_info.dG[k][:zg] .* (model_info.xg .- xⱼ[:zg]))
                                                                                 + sum(function_info.dG[k][:zl] .* (model_info.xl .- xⱼ[:zl]))  
                                                                                 ) 
end



#############################################################################################
#####################################  Main: Level Set Method ###############################
#############################################################################################


function LevelSetMethod_optimization!(  indexSets::IndexSets, 
                                        paramDemand::ParamDemand, 
                                        paramOPF::ParamOPF, 
                                        ẑ::Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}},
                                        randomVariables::RandomVariables;                          ## realization of the random time
                                        levelSetMethodParam::LevelSetMethodParam = levelSetMethodParam, 
                                        ϵ::Float64 = 1e-4, interior_value::Float64 = 0.5, Enhanced_Cut::Bool = true
                                        )

    ######################################################################################################################
    ###############################   auxiliary function for function information   ######################################
    ######################################################################################################################
    ##  μ larger is better
    (μ, λ, threshold, nxt_bound, max_iter, Output, Output_Gap, Adj) = (levelSetMethodParam.μ, levelSetMethodParam.λ, levelSetMethodParam.threshold, levelSetMethodParam.nxt_bound, levelSetMethodParam.max_iter, levelSetMethodParam.Output,levelSetMethodParam.Output_Gap, levelSetMethodParam.Adj)
    (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, indexSets.Ω)
    (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L) 

    x_interior= Dict{Symbol, Vector{Float64}}(:zb => [interior_value for i in 1:length(B)] , 
                                                :zg => [interior_value for i in 1:length(G)], 
                                                :zl => [interior_value for i in 1:length(L)]
                                                )
    # l_interior= [.8 for i in 1:n] - .5 * sum_generator
    # @info "$l_interior"
    # l_interior= [interior_value for i in 1:n]

    f_star = forward_stage2_optimize!(indexSets, 
                                        paramDemand, 
                                        paramOPF, 
                                        ẑ,
                                        randomVariables                       ## realization of the random time
                                        )
    f_star_value = f_star[2]

    backwardInfo = backward_stage2_optimize!(indexSets, 
                                    paramDemand, 
                                    paramOPF, 
                                    ẑ,
                                    randomVariables                          ## realization of the random time
                                    )


    ## collect the information from the objecive f, and constraints G
    function compute_f_G(   x₀::Dict{Symbol, Vector{Float64}}; 
                            Enhanced_Cut::Bool = true, f_star_value::Float64 = f_star_value, 
                            indexSets::IndexSets = indexSets, 
                            paramDemand::ParamDemand = paramDemand, 
                            paramOPF::ParamOPF = paramOPF, randomVariables::RandomVariables = randomVariables, 
                            ẑ::Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}} = ẑ, backwardInfo::BackwardInfo = backwardInfo
                            )

        # objective function
        @objective(backwardInfo.Q, Min,  
                sum( sum(paramDemand.w[d] * paramDemand.demand[t][d] * (1 - backwardInfo.x[d, t]) for d in indexSets.D ) for t in randomVariables.τ:indexSets.T) +
                sum(paramDemand.cb[i] * backwardInfo.νb[i] for i in indexSets.B) + 
                sum(paramDemand.cg[g] * backwardInfo.νg[g] for g in indexSets.G) + 
                sum(paramDemand.cl[l] * backwardInfo.νl[l] for l in indexSets.L) +
                x₀[:zb]' * (ẑ[:zb] .- backwardInfo.zb) + x₀[:zg]' * (ẑ[:zg] .- backwardInfo.zg) + x₀[:zl]' * (ẑ[:zl] .- backwardInfo.zl) 
                + paramDemand.penalty * backwardInfo.slack_variable_b - paramDemand.penalty * backwardInfo.slack_variable_c
                )


        ####################################################### solve the model and display the result ###########################################################
        optimize!(backwardInfo.Q)
        state_variable = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}(:zg => round.(JuMP.value.(backwardInfo.zg)), 
                                                        :zb => round.(JuMP.value.(backwardInfo.zb)), 
                                                        :zl => round.(JuMP.value.(backwardInfo.zl))
                                                        )
        F  = JuMP.objective_value(backwardInfo.Q)
        negative_∇F = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}(:zb => - ẑ[:zb] .+ state_variable[:zb],
                                                                            :zg => - ẑ[:zg] .+ state_variable[:zg],
                                                                            :zl => - ẑ[:zl] .+ state_variable[:zl]
                                                                            )
        F_solution = (F = F, negative_∇F = negative_∇F)
        if Enhanced_Cut
            function_value_info  = Dict(1 => - F_solution.F - 
                                                x₀[:zb]' * (x_interior[:zb] .- ẑ[:zb]) - 
                                                x₀[:zg]' * (x_interior[:zg] .- ẑ[:zg]) - 
                                                x₀[:zl]' * (x_interior[:zl] .- ẑ[:zl]),

                                        2 => Dict{Symbol,  Vector{Float64}}( :zb => F_solution.negative_∇F[:zb] .- (x_interior[:zb] .- ẑ[:zb]),
                                                                    :zg => F_solution.negative_∇F[:zg] .- (x_interior[:zg] .- ẑ[:zg]),
                                                                    :zl => F_solution.negative_∇F[:zl] .- (x_interior[:zl] .- ẑ[:zl])
                                                                    ),
                                        3 => Dict(1 => (1- ϵ) * f_star_value - F_solution.F),
                                        4 => Dict(1 => F_solution.negative_∇F),
                                        )

            # else
            #     function_value_info  = Dict(1 => - F_solution[1] - π' *  sum_generator,
            #                                 2 => - F_solution.negative_∇F - sum_generator,
            #                                 3 => Dict(1 => 0.0 ),
            #                                 4 => Dict(1 => F_solution.negative_∇F * 0),
            #                                 )
        end
        return function_value_info
        ## Com_f = function_value_info[1], Com_grad_f = function_value_info[2], 
        ## Com_G = function_value_info[3], Com_grad_G = function_value_info[4], Com_max_g = function_value_info[3]
    end

    ######################################################################################################################
    ##############################################   level set method   ##################################################
    ######################################################################################################################
    
    x₀ = Dict{Symbol, Vector{Float64}}(      :zb => [1. for i in 1:length(B)], 
                                    :zg => [1. for i in 1:length(G)], 
                                    :zl => [1. for i in 1:length(L)]
                                    )
    iter = 1
    α = 1/2

    ## trajectory
    function_value_info = compute_f_G(x₀; Enhanced_Cut = Enhanced_Cut)
    function_info = FunctionInfo(   Dict(1 => x₀), 
                                    function_value_info[3], 
                                    Dict(1 => function_value_info[1]), 
                                    function_value_info[2], 
                                    function_value_info[4], ## Dict{Int64, Dict{Symbol, Vector}} and index is the components of dG (i == 1)
                                    function_value_info[3]
                                    )

    ## model for oracle
    model_oracle = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => Output, 
            "Threads" => 0)
            )



    @variable(model_oracle, z)
    @variable(model_oracle, xb_oracle[B])
    @variable(model_oracle, xg_oracle[G])
    @variable(model_oracle, xl_oracle[L])
    @variable(model_oracle, y <= 0)

    # para_oracle_bound =  abs(α * function_info.f_his[1] + (1-α) * function_info.G_max_his[1] )
    # @variable(model_oracle, z >= - 10^(ceil(log10(-para_oracle_bound))))
    para_oracle_bound = abs(function_info.f_his[1])
    z_rhs = 15 * 10^(ceil(log10(para_oracle_bound)))
    @constraint(model_oracle, oracle_bound, z >= - z_rhs)

    @objective(model_oracle, Min, z)
    oracle_info = ModelInfo(model_oracle, xb_oracle, xg_oracle, xl_oracle, y, z)



    while true
        if Adj
            param_z_rhs = abs(function_info.f_his[iter])
            if z_rhs <  1.5 * param_z_rhs
                # @info "z level up $(z_rhs/param_z_rhs)"
                z_rhs = 2 * z_rhs
            end

            if z_rhs > 10 * param_z_rhs
                # @info "z level down $(z_rhs/param_z_rhs) "
                z_rhs = .1 * z_rhs
            end 
            set_normalized_rhs(oracle_bound, - z_rhs)  
        end

        add_constraint(function_info, oracle_info, iter)
        optimize!(model_oracle)

        st = termination_status(model_oracle)
        if st != MOI.OPTIMAL
            @info "oracle is infeasible"
            break
        end

        f_star = JuMP.objective_value(model_oracle)

        ## formulate alpha model

        result = Δ_model_formulation(function_info, f_star, iter, Output = Output)
        Δ, a_min, a_max = result[1], result[2], result[3]
        
        ## update α
        if μ/2 <= (α-a_min)/(a_max-a_min) .<= 1-μ/2
            α = α
        else
            α = (a_min+a_max)/2
        end

        ## update level
        w = α * f_star
        W = minimum( α * function_info.f_his[j] + (1-α) * function_info.G_max_his[j] for j in 1:iter) 
        level = w + λ * (W - w)

        if Output_Gap == true
            @info "Gap is $Δ, iter num is $iter, func_val is $( - function_value_info[1]), alpha is $α, w is $w, W is $W"
            @info "Constraint is $(function_info.G_max_his[iter])"
        end
        
        ######################################################################################################################
        #########################################     next iteration point   #################################################
        ######################################################################################################################

        ## obtain the next iteration point
        model_nxt = Model(
            optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => Output, 
            "Threads" => 0)
            )

        @variable(model_nxt, xb[B])
        @variable(model_nxt, xg[G])
        @variable(model_nxt, xl[L])
        @variable(model_nxt, z1 >= - nxt_bound)
        @variable(model_nxt, y1)

        @constraint(model_nxt, level_constraint, α * z1 + (1 - α) * y1 <= level)
        @constraint(model_nxt, z1 .>= function_info.f_his[iter] + 
                                        function_info.df[:zb]' * (xb .- function_info.x_his[iter][:zb]) +
                                        function_info.df[:zg]' * (xg .- function_info.x_his[iter][:zg]) +
                                        function_info.df[:zl]' * (xl .- function_info.x_his[iter][:zl]) 
                                        )
        @constraint(model_nxt, [k in keys(function_info.G)], y1 .>= function_info.G[k] + 
                                                                      sum(function_info.dG[k][:zb] .* (xb .- function_info.x_his[iter][:zb]))
                                                                    + sum(function_info.dG[k][:zg] .* (xg .- function_info.x_his[iter][:zg]))
                                                                    + sum(function_info.dG[k][:zl] .* (xl .- function_info.x_his[iter][:zl]))
                                                                    )

        @objective(model_nxt, Min, sum((xb .- function_info.x_his[iter][:zb]) .* (xb .- function_info.x_his[iter][:zb])) +
                                   sum((xg .- function_info.x_his[iter][:zg]) .* (xg .- function_info.x_his[iter][:zg])) +
                                   sum((xl .- function_info.x_his[iter][:zl]) .* (xl .- function_info.x_his[iter][:zl]))
                                    )
        optimize!(model_nxt)
        st = termination_status(model_nxt)
        if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
            x_nxt = Dict{Symbol, Vector{Float64}}(:zb => JuMP.value.(xb) , 
                                            :zg => JuMP.value.(xg), 
                                            :zl => JuMP.value.(xl)
                                            )
        elseif st == MOI.NUMERICAL_ERROR ## need to figure out why this case happened and fix it
            if Enhanced_Cut
                return [ - function_info.f_his[iter] - function_info.x_his[iter][:zb]' * x_interior[:zb] - 
                                                        function_info.x_his[iter][:zg]' * x_interior[:zg] - 
                                                        function_info.x_his[iter][:zl]' * x_interior[:zl],  function_info.x_his[iter]] 
            # else
            #     return [ - function_info.f_his[iter] - function_info.x_his[iter]' * sum_generator, 
            #                                     function_info.x_his[iter]]
            end
        else
            set_normalized_rhs( level_constraint, w + 1 * (W - w))
            optimize!(model_nxt)
            x_nxt = Dict{Symbol, Vector{Float64}}(:zb => JuMP.value.(xb) , 
                                            :zg => JuMP.value.(xg), 
                                            :zl => JuMP.value.(xl)
                                            )
            # break   
        end

        ## stop rule
        if Δ < threshold || iter > max_iter 
            if Enhanced_Cut
                return [ - function_info.f_his[iter] - function_info.x_his[iter][:zb]' * x_interior[:zb] - 
                                                        function_info.x_his[iter][:zg]' * x_interior[:zg] - 
                                                        function_info.x_his[iter][:zl]' * x_interior[:zl],  function_info.x_his[iter]] 
            # else
            #     return [ - function_info.f_his[iter] - function_info.x_his[iter]' * sum_generator, 
            #                                     function_info.x_his[iter]]
            end
        end
        ######################################################################################################################
        #####################################################    end   #######################################################
        ######################################################################################################################

        ## save the trajectory
        function_value_info = compute_f_G(x_nxt, Enhanced_Cut = Enhanced_Cut)
        iter = iter + 1
        function_info.x_his[iter]     = x_nxt
        function_info.G_max_his[iter] = function_value_info[3][1]
        function_info.f_his[iter]     = function_value_info[1]
        function_info.df              = function_value_info[2]
        function_info.dG              = function_value_info[4]
        function_info.G               = function_value_info[3]

    end

end

