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

function Δ_model_formulation(functionHistory::FunctionHistory, f_star::Float64, iter::Int64; Output::Int64 = 0)
    
    model_alpha = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV),
            "OutputFlag" => Output, 
            "Threads" => 0)
            )

    @variable(model_alpha, z)
    @variable(model_alpha, 0 <= α <= 1)
    @constraint(model_alpha, con[j = 1:iter], z <=  α * ( functionHistory.f_his[j] - f_star) + (1 - α) * functionHistory.G_max_his[j] )
    
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

function add_constraint(currentInfo::CurrentInfo, model_info::ModelInfo)
    m = length(currentInfo.G)

    xⱼ = currentInfo.x
    # add constraints     
    @constraint(model_info.model, model_info.z .>= currentInfo.f + currentInfo.df[:zb]' * (model_info.xb .- xⱼ[:zb]) 
                                                                    + currentInfo.df[:zg]' * (model_info.xg .- xⱼ[:zg] )
                                                                    + currentInfo.df[:zl]' * (model_info.xl .- xⱼ[:zl] )
                                                                    )

    @constraint(model_info.model, [k = 1:m], model_info.y .>= currentInfo.G[k] + sum(currentInfo.dG[k][:zb] .* (model_info.xb .- xⱼ[:zb]))
                                                                                 + sum(currentInfo.dG[k][:zg] .* (model_info.xg .- xⱼ[:zg]))
                                                                                 + sum(currentInfo.dG[k][:zl] .* (model_info.xl .- xⱼ[:zl]))  
                                                                                 ) 
end



#############################################################################################
#####################################  Main: Level Set Method ###############################
#############################################################################################


function LevelSetMethod_optimization!(  indexSets::IndexSets, 
                                        paramDemand::ParamDemand, 
                                        paramOPF::ParamOPF, 
                                        ẑ::Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}, f_star_value::Float64,
                                        randomVariables::RandomVariables;                          ## realization of the random time
                                        levelSetMethodParam::LevelSetMethodParam = levelSetMethodParam, 
                                        ϵ::Float64 = 1e-4, interior_value::Float64 = 0.5, Enhanced_Cut::Bool = true
                                        )

    ######################################################################################################################
    ###############################   auxiliary function for function information   ######################################
    ######################################################################################################################
    ##  μ larger is better
    (μ, λ, threshold, nxt_bound, max_iter, Output, Output_Gap) = (levelSetMethodParam.μ, levelSetMethodParam.λ, levelSetMethodParam.threshold, levelSetMethodParam.nxt_bound, levelSetMethodParam.max_iter, levelSetMethodParam.Output,levelSetMethodParam.Output_Gap)
    (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, indexSets.Ω)
    (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L) 

    x_interior= Dict{Symbol, Vector{Float64}}(:zb => [interior_value for i in 1:length(B)] , 
                                                :zg => [interior_value for i in 1:length(G)], 
                                                :zl => [interior_value for i in 1:length(L)]
                                                )
    # l_interior= [.8 for i in 1:n] - .5 * sum_generator
    # @info "$l_interior"
    # l_interior= [interior_value for i in 1:n]

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

        if Enhanced_Cut
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

            F  = JuMP.objective_value(backwardInfo.Q);
            negative_∇F = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}(:zb => - ẑ[:zb] .+ round.(JuMP.value.(backwardInfo.zb)),
                                                                                :zg => - ẑ[:zg] .+ round.(JuMP.value.(backwardInfo.zg)),
                                                                                :zl => - ẑ[:zl] .+ round.(JuMP.value.(backwardInfo.zl))
                                                                                );

            currentInfo = CurrentInfo( x₀,                                                                 ## current point
                                        - F - 
                                                x₀[:zb]' * (x_interior[:zb] .- ẑ[:zb]) - 
                                                x₀[:zg]' * (x_interior[:zg] .- ẑ[:zg]) - 
                                                x₀[:zl]' * (x_interior[:zl] .- ẑ[:zl]),                     ## obj function value
                                        Dict(1 => (1 - ϵ) * f_star_value - F),                              ## constraint value
                                        Dict{Symbol, Vector{Float64}}(  :zb => negative_∇F[:zb] .- (x_interior[:zb] .- ẑ[:zb]),
                                                                        :zg => negative_∇F[:zg] .- (x_interior[:zg] .- ẑ[:zg]),
                                                                        :zl => negative_∇F[:zl] .- (x_interior[:zl] .- ẑ[:zl])
                                                                        ),                                  ## obj gradient
                                        Dict(1 => negative_∇F )                                             ## constraint gradient
                                        );
        else
            # objective function
            @objective(backwardInfo.Q, Min,  
                    sum( sum(paramDemand.w[d] * paramDemand.demand[t][d] * (1 - backwardInfo.x[d, t]) for d in indexSets.D ) for t in randomVariables.τ:indexSets.T) +
                    sum(paramDemand.cb[i] * backwardInfo.νb[i] for i in indexSets.B) + 
                    sum(paramDemand.cg[g] * backwardInfo.νg[g] for g in indexSets.G) + 
                    sum(paramDemand.cl[l] * backwardInfo.νl[l] for l in indexSets.L) -
                    x₀[:zb]' * backwardInfo.zb - x₀[:zg]' * backwardInfo.zg - x₀[:zl]' * backwardInfo.zl
                    + paramDemand.penalty * backwardInfo.slack_variable_b - paramDemand.penalty * backwardInfo.slack_variable_c
                    );


            ####################################################### solve the model and display the result ###########################################################
            optimize!(backwardInfo.Q)
            F  = JuMP.objective_value(backwardInfo.Q)
            negative_∇F = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}(:zb => round.(JuMP.value.(backwardInfo.zb)),
                                                                                :zg => round.(JuMP.value.(backwardInfo.zg)),
                                                                                :zl => round.(JuMP.value.(backwardInfo.zl))
                                                                                );

            currentInfo = CurrentInfo( x₀,                                                                 ## current point
                                        - F - x₀[:zb]' * ẑ[:zb] - x₀[:zg]' * ẑ[:zg] - x₀[:zl]' * ẑ[:zl],    ## obj function value
                                        Dict(1 => 0.0 ), ## constraint value
                                        Dict{Symbol, Vector{Float64}}(  :zb => negative_∇F[:zb] .- ẑ[:zb],
                                                                        :zg => negative_∇F[:zg] .- ẑ[:zg],
                                                                        :zl => negative_∇F[:zl] .- ẑ[:zl]
                                                                        ),                                  ## obj gradient
                                        Dict{Int64, Dict{Symbol, Vector{Float64}}}(1 => 
                                        Dict{Symbol, Vector{Float64}}( :zb => zeros(size(B)),
                                                                                                :zg => zeros(size(G)),
                                                                                                :zl => zeros(size(L))
                                                                                                ))          ## constraint gradient
                                        )
        end
        return currentInfo
    end    
    ######################################################################################################################
    ##############################################   level set method   ##################################################
    ######################################################################################################################
    if Enhanced_Cut
        x₀ = Dict{Symbol, Vector{Float64}}(      :zb => ẑ[:zb] * 4 * f_star_value .- 2 * f_star_value, 
                                        :zg => ẑ[:zg] * 4 * f_star_value .- 2 * f_star_value, 
                                        :zl => ẑ[:zl] * 4 * f_star_value .- 2 * f_star_value
                                        )
    else
        x₀ = Dict{Symbol, Vector{Float64}}(      :zb => ẑ[:zb] * 0, 
                                        :zg => ẑ[:zg] * 0, 
                                        :zl => ẑ[:zl] * 0
                                        )
    end
    # x₀ = LevelSetMethod_InitialPoint()[2]
                                    
    iter = 1
    α = 1/2

    ## trajectory
    currentInfo = compute_f_G(x₀; Enhanced_Cut = Enhanced_Cut);

    functionHistory = FunctionHistory(  Dict(1 => currentInfo.f), 
                                        Dict(1 => maximum(currentInfo.G[k] for k in keys(currentInfo.G)) )
                                        );

    ## model for oracle
    model_oracle = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => Output, 
            "Threads" => 0)
            );



    @variable(model_oracle, z);
    @variable(model_oracle, xb_oracle[B]);
    @variable(model_oracle, xg_oracle[G]);
    @variable(model_oracle, xl_oracle[L]);
    @variable(model_oracle, y <= 0);

    para_oracle_bound = abs(currentInfo.f);
    z_rhs = 10 * 10^(ceil(log10(para_oracle_bound)));
    @constraint(model_oracle, oracle_bound, z >= - z_rhs);

    @objective(model_oracle, Min, z);
    oracle_info = ModelInfo(model_oracle, xb_oracle, xg_oracle, xl_oracle, y, z);


    gap_list = []
    iter_significance = Inf

    while true
        add_constraint(currentInfo, oracle_info)
        optimize!(model_oracle)

        st = termination_status(model_oracle)
        if st != MOI.OPTIMAL
            @info "Break --oracle is infeasible"
            break
        end

        f_star = JuMP.objective_value(model_oracle)

        ## formulate alpha model

        result = Δ_model_formulation(functionHistory, f_star, iter, Output = Output)
        Δ, a_min, a_max = result[1], result[2], result[3]

        push!(gap_list, Δ)
        
        ## update α
        if μ/2 <= (α-a_min)/(a_max-a_min) .<= 1-μ/2
            α = α
        else
            α = (a_min+a_max)/2
        end

        ## update level
        w = α * f_star
        W = minimum( α * functionHistory.f_his[j] + (1-α) * functionHistory.G_max_his[j] for j in 1:iter) 
        level = w + λ * (W - w)

        if Output_Gap
            @info "Gap is $Δ, iter num is $iter, func_val is $( - currentInfo.f), Constraint is $(currentInfo.G)"
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
        @constraint(model_nxt, z1 .>= currentInfo.f + 
                                        currentInfo.df[:zb]' * (xb .- currentInfo.x[:zb]) +
                                        currentInfo.df[:zg]' * (xg .- currentInfo.x[:zg]) +
                                        currentInfo.df[:zl]' * (xl .- currentInfo.x[:zl]) 
                                        )
        @constraint(model_nxt, [k in keys(currentInfo.G)], y1 .>= currentInfo.G[k] + 
                                                                      sum(currentInfo.dG[k][:zb] .* (xb .- currentInfo.x[:zb]))
                                                                    + sum(currentInfo.dG[k][:zg] .* (xg .- currentInfo.x[:zg]))
                                                                    + sum(currentInfo.dG[k][:zl] .* (xl .- currentInfo.x[:zl]))
                                                                    )

        @objective(model_nxt, Min, sum((xb .- currentInfo.x[:zb]) .* (xb .- currentInfo.x[:zb])) +
                                   sum((xg .- currentInfo.x[:zg]) .* (xg .- currentInfo.x[:zg])) +
                                   sum((xl .- currentInfo.x[:zl]) .* (xl .- currentInfo.x[:zl]))
                                    )
        optimize!(model_nxt)
        st = termination_status(model_nxt)
        if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
            x_nxt = Dict{Symbol, Vector{Float64}}(:zb => JuMP.value.(xb) , 
                                            :zg => JuMP.value.(xg), 
                                            :zl => JuMP.value.(xl)
                                            )
        elseif st == MOI.NUMERICAL_ERROR ## need to figure out why this case happened and fix it
            @info "Termination -- Numerical Error occures!"
            if Enhanced_Cut
                return [ - currentInfo.f - currentInfo.x[:zb]' * x_interior[:zb] - 
                                                        currentInfo.x[:zg]' * x_interior[:zg] - 
                                                        currentInfo.x[:zl]' * x_interior[:zl],  currentInfo.x] 
            else
                return [ - currentInfo.f - currentInfo.x[:zb]' * ẑ[:zb] - 
                                                currentInfo.x[:zg]' * ẑ[:zg] - 
                                                currentInfo.x[:zl]' * ẑ[:zl],  currentInfo.x] 
            end
        else
            @info "Re-compute Next Iteration Point -- change to a safe level!"
            set_normalized_rhs( level_constraint, w + 1 * (W - w))
            optimize!(model_nxt)
            x_nxt = Dict{Symbol, Vector{Float64}}(:zb => JuMP.value.(xb) , 
                                            :zg => JuMP.value.(xg), 
                                            :zl => JuMP.value.(xl)
                                            )
            # break   
        end

        if iter > 30
            iter_significance = abs(gap_list[iter] - sum(gap_list[iter-29:iter])/30)
            if iter_significance ≤ Δ * 1e-2 && currentInfo.G[1] ≤ 0
                @info "Termination -- Progress ($iter_significance) is too slow."
                if Enhanced_Cut
                    return [ - currentInfo.f - currentInfo.x[:zb]' * x_interior[:zb] - 
                                                            currentInfo.x[:zg]' * x_interior[:zg] - 
                                                            currentInfo.x[:zl]' * x_interior[:zl],  currentInfo.x] 
                else
                    return [ - currentInfo.f - currentInfo.x[:zb]' * ẑ[:zb] - 
                                                    currentInfo.x[:zg]' * ẑ[:zg] - 
                                                    currentInfo.x[:zl]' * ẑ[:zl],  currentInfo.x] 
                end
            end
        end

        ## stop rule
        if Δ < threshold * f_star_value || iter > max_iter
            if Enhanced_Cut
                return [ - currentInfo.f - currentInfo.x[:zb]' * x_interior[:zb] - 
                                                        currentInfo.x[:zg]' * x_interior[:zg] - 
                                                        currentInfo.x[:zl]' * x_interior[:zl],  currentInfo.x] 
            else
                return [ - currentInfo.f - currentInfo.x[:zb]' * ẑ[:zb] - 
                                                currentInfo.x[:zg]' * ẑ[:zg] - 
                                                currentInfo.x[:zl]' * ẑ[:zl],  currentInfo.x] 
            end
        end
        ######################################################################################################################
        #####################################################    end   #######################################################
        ######################################################################################################################

        ## save the trajectory
        currentInfo = compute_f_G(x_nxt, Enhanced_Cut = Enhanced_Cut)
        iter = iter + 1
        
        functionHistory.f_his[iter] = currentInfo.f
        functionHistory.G_max_his[iter] = maximum(currentInfo.G[k] for k in keys(currentInfo.G))
    end

end






## ============================================= Design for solve Lagrangian cut ================================================= ##

# function LevelSetMethod_InitialPoint(  ;indexSets::IndexSets = indexSets, 
#                                         paramDemand::ParamDemand = paramDemand, 
#                                         paramOPF::ParamOPF = paramOPF, 
#                                         ẑ::Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}} = ẑ, f_star_value::Float64 = f_star_value,
#                                         randomVariables::RandomVariables = randomVariables,                          ## realization of the random time
#                                         levelSetMethodParam::LevelSetMethodParam = levelSetMethodParam
#                                         )

#     ######################################################################################################################
#     ###############################   auxiliary function for function information   ######################################
#     ######################################################################################################################
#     ##  μ larger is better
#     (μ, λ, threshold, nxt_bound, max_iter, Output, Output_Gap, Adj) = (levelSetMethodParam.μ, levelSetMethodParam.λ, levelSetMethodParam.threshold, levelSetMethodParam.nxt_bound, levelSetMethodParam.max_iter, levelSetMethodParam.Output,levelSetMethodParam.Output_Gap, levelSetMethodParam.Adj)
#     (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, indexSets.Ω)
#     (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L) 

#     x_interior= Dict{Symbol, Vector{Float64}}(:zb => [interior_value for i in 1:length(B)] , 
#                                                 :zg => [interior_value for i in 1:length(G)], 
#                                                 :zl => [interior_value for i in 1:length(L)]
#                                                 )
#     # l_interior= [.8 for i in 1:n] - .5 * sum_generator
#     # @info "$l_interior"
#     # l_interior= [interior_value for i in 1:n]

#     backwardInfo = backward_stage2_optimize!(indexSets, 
#                                     paramDemand, 
#                                     paramOPF, 
#                                     ẑ,
#                                     randomVariables                          ## realization of the random time
#                                     )


#     ## collect the information from the objecive f, and constraints G
#     function compute_f_G(   x₀::Dict{Symbol, Vector{Float64}}; 
#                             indexSets::IndexSets = indexSets, 
#                             paramDemand::ParamDemand = paramDemand, 
#                             paramOPF::ParamOPF = paramOPF, randomVariables::RandomVariables = randomVariables, 
#                             ẑ::Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}} = ẑ, backwardInfo::BackwardInfo = backwardInfo
#                             )

#         @objective(backwardInfo.Q, Min,  
#                 sum( sum(paramDemand.w[d] * paramDemand.demand[t][d] * (1 - backwardInfo.x[d, t]) for d in indexSets.D ) for t in randomVariables.τ:indexSets.T) +
#                 sum(paramDemand.cb[i] * backwardInfo.νb[i] for i in indexSets.B) + 
#                 sum(paramDemand.cg[g] * backwardInfo.νg[g] for g in indexSets.G) + 
#                 sum(paramDemand.cl[l] * backwardInfo.νl[l] for l in indexSets.L) -
#                 x₀[:zb]' * backwardInfo.zb - x₀[:zg]' * backwardInfo.zg - x₀[:zl]' * backwardInfo.zl
#                 + paramDemand.penalty * backwardInfo.slack_variable_b - paramDemand.penalty * backwardInfo.slack_variable_c
#                 );


#         ####################################################### solve the model and display the result ###########################################################
#         optimize!(backwardInfo.Q)
#         F  = JuMP.objective_value(backwardInfo.Q)
#         negative_∇F = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}(:zb => round.(JuMP.value.(backwardInfo.zb)),
#                                                                             :zg => round.(JuMP.value.(backwardInfo.zg)),
#                                                                             :zl => round.(JuMP.value.(backwardInfo.zl))
#                                                                             );

#         currentInfo = CurrentInfo( x₀,                                                                 ## current point
#                                     - F - x₀[:zb]' * ẑ[:zb] - x₀[:zg]' * ẑ[:zg] - x₀[:zl]' * ẑ[:zl],    ## obj function value
#                                     Dict(1 => 0.0 ), ## constraint value
#                                     Dict{Symbol, Vector{Float64}}(  :zb => negative_∇F[:zb] .- ẑ[:zb],
#                                                                     :zg => negative_∇F[:zg] .- ẑ[:zg],
#                                                                     :zl => negative_∇F[:zl] .- ẑ[:zl]
#                                                                     ),                                  ## obj gradient
#                                     Dict{Int64, Dict{Symbol, Vector{Float64}}}(1 => 
#                                     Dict{Symbol, Vector{Float64}}( :zb => zeros(size(B)),
#                                                                                             :zg => zeros(size(G)),
#                                                                                             :zl => zeros(size(L))
#                                                                                             ))          ## constraint gradient
#                                     )
        
#         return currentInfo
#     end    
#     ######################################################################################################################
#     ##############################################   level set method   ##################################################
#     ######################################################################################################################

#     x₀ = Dict{Symbol, Vector{Float64}}(      :zb => ẑ[:zb] * 0, 
#                                     :zg => ẑ[:zg] * 0, 
#                                     :zl => ẑ[:zl] * 0
#                                     )
                                    
#     iter = 1
#     α = 1/2

#     ## trajectory
#     currentInfo = compute_f_G(x₀);

#     functionHistory = FunctionHistory(  Dict(1 => currentInfo.f), 
#                                         Dict(1 => maximum(currentInfo.G[k] for k in keys(currentInfo.G)) )
#                                         );

#     ## model for oracle
#     model_oracle = Model(
#         optimizer_with_attributes(
#             ()->Gurobi.Optimizer(GRB_ENV), 
#             "OutputFlag" => Output, 
#             "Threads" => 0)
#             );



#     @variable(model_oracle, z);
#     @variable(model_oracle, xb_oracle[B]);
#     @variable(model_oracle, xg_oracle[G]);
#     @variable(model_oracle, xl_oracle[L]);
#     @variable(model_oracle, y <= 0);

#     para_oracle_bound = abs(currentInfo.f);
#     z_rhs = 15 * 10^(ceil(log10(para_oracle_bound)));
#     @constraint(model_oracle, oracle_bound, z >= - z_rhs);

#     @objective(model_oracle, Min, z);
#     oracle_info = ModelInfo(model_oracle, xb_oracle, xg_oracle, xl_oracle, y, z);



#     while true
#         add_constraint(currentInfo, oracle_info)
#         optimize!(model_oracle)

#         st = termination_status(model_oracle)
#         if st != MOI.OPTIMAL
#             @info "oracle is infeasible"
#             break
#         end

#         f_star = JuMP.objective_value(model_oracle)

#         ## formulate alpha model

#         result = Δ_model_formulation(functionHistory, f_star, iter, Output = Output)
#         Δ, a_min, a_max = result[1], result[2], result[3]
        
#         ## update α
#         if μ/2 <= (α-a_min)/(a_max-a_min) .<= 1-μ/2
#             α = α
#         else
#             α = (a_min+a_max)/2
#         end

#         ## update level
#         w = α * f_star
#         W = minimum( α * functionHistory.f_his[j] + (1-α) * functionHistory.G_max_his[j] for j in 1:iter) 
#         level = w + λ * (W - w)

#         if Output_Gap == true
#             @info "Gap is $Δ, iter num is $iter, func_val is $( - currentInfo.f), Constraint is $(currentInfo.G)"
#         end
        
#         ######################################################################################################################
#         #########################################     next iteration point   #################################################
#         ######################################################################################################################

#         ## obtain the next iteration point
#         model_nxt = Model(
#             optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
#             "OutputFlag" => Output, 
#             "Threads" => 0)
#             )

#         @variable(model_nxt, xb[B])
#         @variable(model_nxt, xg[G])
#         @variable(model_nxt, xl[L])
#         @variable(model_nxt, z1 >= - nxt_bound)
#         @variable(model_nxt, y1)

#         @constraint(model_nxt, level_constraint, α * z1 + (1 - α) * y1 <= level)
#         @constraint(model_nxt, z1 .>= currentInfo.f + 
#                                         currentInfo.df[:zb]' * (xb .- currentInfo.x[:zb]) +
#                                         currentInfo.df[:zg]' * (xg .- currentInfo.x[:zg]) +
#                                         currentInfo.df[:zl]' * (xl .- currentInfo.x[:zl]) 
#                                         )
#         @constraint(model_nxt, [k in keys(currentInfo.G)], y1 .>= currentInfo.G[k] + 
#                                                                       sum(currentInfo.dG[k][:zb] .* (xb .- currentInfo.x[:zb]))
#                                                                     + sum(currentInfo.dG[k][:zg] .* (xg .- currentInfo.x[:zg]))
#                                                                     + sum(currentInfo.dG[k][:zl] .* (xl .- currentInfo.x[:zl]))
#                                                                     )

#         @objective(model_nxt, Min, sum((xb .- currentInfo.x[:zb]) .* (xb .- currentInfo.x[:zb])) +
#                                    sum((xg .- currentInfo.x[:zg]) .* (xg .- currentInfo.x[:zg])) +
#                                    sum((xl .- currentInfo.x[:zl]) .* (xl .- currentInfo.x[:zl]))
#                                     )
#         optimize!(model_nxt)
#         st = termination_status(model_nxt)
#         if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
#             x_nxt = Dict{Symbol, Vector{Float64}}(:zb => JuMP.value.(xb) , 
#                                             :zg => JuMP.value.(xg), 
#                                             :zl => JuMP.value.(xl)
#                                             )
#         elseif st == MOI.NUMERICAL_ERROR ## need to figure out why this case happened and fix it
#             return [ - currentInfo.f - currentInfo.x[:zb]' * ẑ[:zb] - 
#                                                 currentInfo.x[:zg]' * ẑ[:zg] - 
#                                                 currentInfo.x[:zl]' * ẑ[:zl],  currentInfo.x] 
#         else
#             set_normalized_rhs( level_constraint, w + 1 * (W - w))
#             optimize!(model_nxt)
#             x_nxt = Dict{Symbol, Vector{Float64}}(:zb => JuMP.value.(xb) , 
#                                             :zg => JuMP.value.(xg), 
#                                             :zl => JuMP.value.(xl)
#                                             )
#             break   
#         end

#         ## stop rule
#         if Δ < threshold * f_star_value || iter > max_iter 
#             return [ - currentInfo.f - currentInfo.x[:zb]' * ẑ[:zb] - 
#                                                 currentInfo.x[:zg]' * ẑ[:zg] - 
#                                                 currentInfo.x[:zl]' * ẑ[:zl],  currentInfo.x] 
#         end
#         ######################################################################################################################
#         #####################################################    end   #######################################################
#         ######################################################################################################################

#         ## save the trajectory
#         currentInfo = compute_f_G(x_nxt, Enhanced_Cut = Enhanced_Cut)
#         iter = iter + 1
        
#         functionHistory.f_his[iter] = currentInfo.f
#         functionHistory.G_max_his[iter] = maximum(currentInfo.G[k] for k in keys(currentInfo.G))
#     end

# end
