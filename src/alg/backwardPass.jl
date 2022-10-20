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
                                    randomVariables::RandomVariables;                          ## realization of the random time
                                    timelimit::Int64 = 5, 
                                    outputFlag::Int64 = 0,
                                    tightness::Bool = false
                                    )

    # (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, indexSets.Ω)
    # (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L) 

    Q = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                "OutputFlag" => outputFlag, 
                "Threads" => 0,
                "MIPGap" => 1e-4, 
                "TimeLimit" => timelimit ) 
                )
                

    @variable(Q, θ_angle[indexSets.B, 1:indexSets.T]) 
    @variable(Q, P[indexSets.L, 1:indexSets.T])                     ## elements in L is Tuple (i, j)
    @variable(Q, s[indexSets.G, 1:indexSets.T] ≥ 0)
    @variable(Q, 0 ≤ x[indexSets.D, 1:indexSets.T] ≤ 1)

    @variable(Q, yb[indexSets.B], Bin)
    @variable(Q, yg[indexSets.G], Bin)
    @variable(Q, yl[indexSets.L], Bin)

    @variable(Q, νb[indexSets.B], Bin)
    @variable(Q, νg[indexSets.G], Bin)
    @variable(Q, νl[indexSets.L], Bin)

    if tightness
        @variable(Q, zg[indexSets.G], Bin)
        @variable(Q, zb[indexSets.B], Bin)
        @variable(Q, zl[indexSets.L], Bin)
    else 
        @variable(Q, 0 ≤ zg[indexSets.G] ≤ 1)
        @variable(Q, 0 ≤ zb[indexSets.B] ≤ 1)
        @variable(Q, 0 ≤ zl[indexSets.L] ≤ 1)
    end


    ## constraint k l m 
    @constraint(Q, [i in indexSets.B], yb[i] ≤ zb[i] )
    @constraint(Q, [g in indexSets.G], yg[g] ≤ zg[g] )
    @constraint(Q, [l in indexSets.L], yl[l] ≤ zl[l] )

    @constraint(Q, [i in indexSets.B], yb[i] ≤ 1- νb[i] )
    @constraint(Q, [g in indexSets.G], yg[g] ≤ 1- νg[g] )
    @constraint(Q, [l in indexSets.L], yl[l] ≤ 1- νl[l] )

    @constraint(Q, [i in indexSets.B], νb[i] ≥ randomVariables.vb[i] )
    @constraint(Q, [g in indexSets.G], νg[g] ≥ randomVariables.vg[g] )
    @constraint(Q, [l in indexSets.L], νl[l] ≥ randomVariables.vl[l] )

    
    for i in indexSets.B 

        ## constraint 3e
        @constraint(Q, [t in randomVariables.τ:indexSets.T], 
                                                        sum(s[g, t] for g in indexSets.Gᵢ[i]) - 
                                                                sum(P[(i, j), t] for j in indexSets.out_L[i]) + 
                                                                    sum(P[(j, i), t] for j in indexSets.in_L[i]) 
                                                                        .== sum(paramDemand.demand[t][d] * x[d, t] for d in indexSets.Dᵢ[i]) )

        ## constraint n
        @constraint(Q, [j in unique(randomVariables.Ibb[i])], νb[j] ≥ randomVariables.ub[i] * zb[i] )
        @constraint(Q, [j in unique(randomVariables.Ibg[i])], νg[j] ≥ randomVariables.ub[i] * zb[i] )
        @constraint(Q, [j in unique(randomVariables.Ibl[i])], νl[j] ≥ randomVariables.ub[i] * zb[i] )

         ## constraint g h i j
        @constraint(Q, [t in randomVariables.τ:indexSets.T, d in indexSets.Dᵢ[i]], yb[i] ≥ x[d, t] )
        @constraint(Q, [g in indexSets.Gᵢ[i]], yb[i] ≥ yg[g])
        @constraint(Q, [j in indexSets.out_L[i]], yb[i] ≥ yl[(i, j)] )
        @constraint(Q, [j in indexSets.in_L[i]], yb[i] ≥ yl[(j, i)] )
    end

    for g in indexSets.G
        ## constraint n 
        @constraint(Q, [j in unique(randomVariables.Igb[g])], νb[j] ≥ randomVariables.ug[g] * zg[g] )
        @constraint(Q, [j in unique(randomVariables.Igg[g])], νg[j] ≥ randomVariables.ug[g] * zg[g] )
        @constraint(Q, [j in unique(randomVariables.Igl[g])], νl[j] ≥ randomVariables.ug[g] * zg[g] )

        ## constraint 3f
        @constraint(Q, [g in indexSets.G, t in randomVariables.τ:indexSets.T], s[g, t] ≥ paramOPF.smin[g] * yg[g])
        @constraint(Q, [g in indexSets.G, t in randomVariables.τ:indexSets.T], s[g, t] ≤ paramOPF.smax[g] * yg[g])
    end

    for l in indexSets.L 
        ## constraint n
        @constraint(Q, [j in unique(randomVariables.Ilb[l])], νb[j] ≥ randomVariables.ul[l] * zl[l] )
        @constraint(Q, [j in unique(randomVariables.Ilg[l])], νg[j] ≥ randomVariables.ul[l] * zl[l] )
        @constraint(Q, [j in unique(randomVariables.Ill[l])], νl[j] ≥ randomVariables.ul[l] * zl[l] )

        ## constraint 3b 3c
        i = l[1]
        j = l[2]
        @constraint(Q, [t in randomVariables.τ:indexSets.T], P[l, t] ≤ - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmax * (1 - yl[l] ) ) )
        @constraint(Q, [t in randomVariables.τ:indexSets.T], P[l, t] ≥ - paramOPF.b[l] * (θ_angle[i, t] - θ_angle[j, t] + paramOPF.θmin * (1 - yl[l] ) ) )

        ## constraint 3d
        @constraint(Q, [t in randomVariables.τ:indexSets.T], P[l, t] ≥ - paramOPF.W[l] * yl[l] )
        @constraint(Q, [t in randomVariables.τ:indexSets.T], P[l, t] ≤ paramOPF.W[l] * yl[l] )
    end

    ## objective function
    # @objective(Q, Min,  
    #         sum( sum(paramDemand.w[d] * (1 - x[d, t]) for d in indexSets.D ) for t in randomVariables.τ:indexSets.T) +
    #         sum(paramDemand.cb[i] * νb[i] for i in indexSets.B) + 
    #         sum(paramDemand.cg[g] * νg[g] for g in indexSets.G) + 
    #         sum(paramDemand.cl[l] * νl[l] for l in indexSets.L) +
    #         π[:zb]' * (ẑ[:zb] .- zb) + π[:zg]' * (ẑ[:zg] .- zg) + π[:zl]' * (ẑ[:zl] .- zl) 
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

    backwardInfo = BackwardInfo(Q, x, νb, νg, νl, zg, zb, zl)
    return backwardInfo
end


#############################################################################################
###########################  auxiliary functions for level set method #######################
#############################################################################################

"""
    This function is to constraint the model for solving gap and alpha
"""

function Δ_model_formulation(functionHistory::FunctionHistory, f_star::Float64, iter::Int64; Output::Int64 = 0)
    
    alphaModel = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV),
            "OutputFlag" => Output, 
            "Threads" => 0)
            )

    @variable(alphaModel, z)
    @variable(alphaModel, 0 ≤ α ≤ 1)
    @constraint(alphaModel, con[j = 1:iter], z ≤  α * ( functionHistory.f_his[j] - f_star) + (1 - α) * functionHistory.G_max_his[j] )
    
    # we first compute gap Δ
    @objective(alphaModel, Max, z)
    optimize!(alphaModel)
    st = termination_status(alphaModel)
    Δ = JuMP.objective_value(alphaModel)

    
    ## then we modify above model to compute alpha
    # alpha_min
    @constraint(alphaModel, z .≥ 0)
    @objective(alphaModel, Min, α)
    optimize!(alphaModel)
    a_min = JuMP.value(α)

    # alpha_max
    @objective(alphaModel, Max, α)
    optimize!(alphaModel)
    a_max = JuMP.value(α)

    return Dict(1 => Δ, 2 => a_min, 3 => a_max)

end


"""
    This function is to add constraints for the model f_star and nxt pt.
"""

function add_constraint(currentInfo::CurrentInfo, modelInfo::ModelInfo)
    m = length(currentInfo.G)

    xⱼ = currentInfo.x
    # add constraints     
    @constraint(modelInfo.model, modelInfo.z .≥ currentInfo.f + currentInfo.df[:zb]' * (modelInfo.xb .- xⱼ[:zb]) 
                                                                    + currentInfo.df[:zg]' * (modelInfo.xg .- xⱼ[:zg] )
                                                                    + currentInfo.df[:zl]' * (modelInfo.xl .- xⱼ[:zl] )
                                                                    )

    @constraint(modelInfo.model, [k = 1:m], modelInfo.y .≥ currentInfo.G[k] + sum(currentInfo.dG[k][:zb] .* (modelInfo.xb .- xⱼ[:zb]))
                                                                                 + sum(currentInfo.dG[k][:zg] .* (modelInfo.xg .- xⱼ[:zg]))
                                                                                 + sum(currentInfo.dG[k][:zl] .* (modelInfo.xl .- xⱼ[:zl]))  
                                                                                 ) 
end



#############################################################################################
#####################################  Main: Level Set Method ###############################
#############################################################################################
function LevelSetMethod_optimization!(ẑ::Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}, f_star_value::Float64, randomVariables::RandomVariables;
                                        indexSets::IndexSets = indexSets, paramDemand::ParamDemand = paramDemand, paramOPF::ParamOPF = paramOPF, ϵ::Float64 = 1e-4, 
                                        levelSetMethodParam::LevelSetMethodParam = levelSetMethodParam, 
                                        cutSelection::String = cutSelection,## "ELC", "LC", "ShrinkageLC" 
                                        x_interior::Union{Dict{Symbol, Vector{Float64}}, Nothing} = nothing, 
                                        x₀::Dict{Symbol, Vector{Float64}} = x₀, tightness::Bool = false
                                        )

    ## ==================================================== auxiliary function for function information ==================================================== ##
    # μ larger is better
    (μ, λ, threshold, nxt_bound, max_iter, Output, Output_Gap) = (levelSetMethodParam.μ, levelSetMethodParam.λ, levelSetMethodParam.threshold, levelSetMethodParam.nxt_bound, levelSetMethodParam.max_iter, levelSetMethodParam.Output,levelSetMethodParam.Output_Gap);
    (D, G, L, B) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B);

    backwardInfo = backward_stage2_optimize!(indexSets, 
                                    paramDemand, 
                                    paramOPF, 
                                    ẑ,
                                    randomVariables;                       ## realization of the random time
                                    outputFlag = 0,
                                    timelimit = 8,
                                    tightness = tightness
                                    );


    # collect the information from the objecive f, and constraints G
    function compute_f_G(   x₀::Dict{Symbol, Vector{Float64}}; 
                            cutSelection::String = "ELC", f_star_value::Float64 = f_star_value, 
                            indexSets::IndexSets = indexSets, 
                            paramDemand::ParamDemand = paramDemand, 
                            paramOPF::ParamOPF = paramOPF, randomVariables::RandomVariables = randomVariables, 
                            ẑ::Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}} = ẑ, backwardInfo::BackwardInfo = backwardInfo
                            )

        if cutSelection == "ELC"
            # objective function
            @objective(backwardInfo.model, Min,  
                    sum( sum(paramDemand.w[d] * (1 - backwardInfo.x[d, t]) for d in indexSets.D ) for t in randomVariables.τ:indexSets.T) +
                    sum(paramDemand.cb[i] * backwardInfo.νb[i] for i in indexSets.B) + 
                    sum(paramDemand.cg[g] * backwardInfo.νg[g] for g in indexSets.G) + 
                    sum(paramDemand.cl[l] * backwardInfo.νl[l] for l in indexSets.L) +
                    x₀[:zb]' * (ẑ[:zb] .- backwardInfo.zb) + x₀[:zg]' * (ẑ[:zg] .- backwardInfo.zg) + x₀[:zl]' * (ẑ[:zl] .- backwardInfo.zl) 
                    )


            ## ==================================================== solve the model and display the result ==================================================== ##
            optimize!(backwardInfo.model)

            F  = JuMP.objective_value(backwardInfo.model);
            negative_∇F = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}(:zb => - ẑ[:zb] .+ round.(JuMP.value.(backwardInfo.zb), digits = 5),
                                                                                :zg => - ẑ[:zg] .+ round.(JuMP.value.(backwardInfo.zg), digits = 5),
                                                                                :zl => - ẑ[:zl] .+ round.(JuMP.value.(backwardInfo.zl), digits = 5)
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
        elseif cutSelection == "LC"
            # objective function
            @objective(backwardInfo.model, Min,  
                    sum( sum(paramDemand.w[d] * (1 - backwardInfo.x[d, t]) for d in indexSets.D ) for t in randomVariables.τ:indexSets.T) +
                    sum(paramDemand.cb[i] * backwardInfo.νb[i] for i in indexSets.B) + 
                    sum(paramDemand.cg[g] * backwardInfo.νg[g] for g in indexSets.G) + 
                    sum(paramDemand.cl[l] * backwardInfo.νl[l] for l in indexSets.L) -
                    x₀[:zb]' * backwardInfo.zb - x₀[:zg]' * backwardInfo.zg - x₀[:zl]' * backwardInfo.zl
                    );
            ## ==================================================== solve the model and display the result ==================================================== ##
            optimize!(backwardInfo.model)
            F  = JuMP.objective_value(backwardInfo.model)
            negative_∇F = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}(:zb => round.(JuMP.value.(backwardInfo.zb), digits = 5),
                                                                                :zg => round.(JuMP.value.(backwardInfo.zg), digits = 5),
                                                                                :zl => round.(JuMP.value.(backwardInfo.zl), digits = 5)
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
        elseif cutSelection == "ShrinkageLC"
            @objective(backwardInfo.model, Min,  
                                        sum( sum(paramDemand.w[d] * (1 - backwardInfo.x[d, t]) for d in indexSets.D ) for t in randomVariables.τ:indexSets.T) +
                                        sum(paramDemand.cb[i] * backwardInfo.νb[i] for i in indexSets.B) + 
                                        sum(paramDemand.cg[g] * backwardInfo.νg[g] for g in indexSets.G) + 
                                        sum(paramDemand.cl[l] * backwardInfo.νl[l] for l in indexSets.L) +
                                        x₀[:zb]' * (ẑ[:zb] .- backwardInfo.zb) + x₀[:zg]' * (ẑ[:zg] .- backwardInfo.zg) + x₀[:zl]' * (ẑ[:zl] .- backwardInfo.zl) 
                                        )

            ## ==================================================== solve the model and display the result ==================================================== ##
            optimize!(backwardInfo.model)

            F  = JuMP.objective_value(backwardInfo.model);
            negative_∇F = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}(:zb => - ẑ[:zb] .+ round.(JuMP.value.(backwardInfo.zb), digits = 5),
                                                                                :zg => - ẑ[:zg] .+ round.(JuMP.value.(backwardInfo.zg), digits = 5),
                                                                                :zl => - ẑ[:zl] .+ round.(JuMP.value.(backwardInfo.zl), digits = 5)
                                                                                );

            currentInfo = CurrentInfo( x₀,                                                                                                ## current point
                                        1/2 * sum(x₀[:zb] .* x₀[:zb]) + 1/2 * sum(x₀[:zg] .* x₀[:zg]) + 1/2 * sum(x₀[:zl] .* x₀[:zl]) ,   ## obj function value
                                        Dict(1 => (1 - ϵ) * f_star_value - F),                                                            ## constraint value
                                        x₀,                                                                                               ## obj gradient
                                        Dict(1 => negative_∇F )                                                                           ## constraint gradient
                                        );
        end
        return (currentInfo = currentInfo, currentInfo_f = F)
    end    
    
    ## ==================================================== Levelset Method ============================================== ##    
    iter = 1;
    α = 1/2;

    # trajectory
    currentInfo, currentInfo_f = compute_f_G(x₀; cutSelection = cutSelection, f_star_value = f_star_value);

    functionHistory = FunctionHistory(  Dict(1 => currentInfo.f), 
                                        Dict(1 => maximum(currentInfo.G[k] for k in keys(currentInfo.G)) )
                                        );

    # model for oracle
    oracleModel = Model(
        optimizer_with_attributes(
            ()->Gurobi.Optimizer(GRB_ENV), 
            "OutputFlag" => Output, 
            "Threads" => 0)
            );

    ## ==================================================== Levelset Method ============================================== ##
    para_oracle_bound = abs(currentInfo.f);
    z_rhs = 5 * 10^(ceil(log10(para_oracle_bound)));
    @variable(oracleModel, z ≥ - z_rhs);
    @variable(oracleModel, xb_oracle[B]);
    @variable(oracleModel, xg_oracle[G]);
    @variable(oracleModel, xl_oracle[L]);
    @variable(oracleModel, y ≤ 0);

    @objective(oracleModel, Min, z);
    oracleInfo = ModelInfo(oracleModel, xb_oracle, xg_oracle, xl_oracle, y, z);


    nxtModel = Model(
        optimizer_with_attributes(
        ()->Gurobi.Optimizer(GRB_ENV), 
        "OutputFlag" => Output, 
        "Threads" => 0)
        )

    @variable(nxtModel, xb[B]);
    @variable(nxtModel, xg[G]);
    @variable(nxtModel, xl[L]);
    @variable(nxtModel, z1);
    @variable(nxtModel, y1);
    nxtInfo = ModelInfo(nxtModel, xb, xg, xl, y1, z1);

    Δ = Inf; τₖ = 1; τₘ = .5; μₖ = 1;

    if cutSelection == "ELC"
        cutInfo =  [ - currentInfo.f - currentInfo.x[:zb]' * x_interior[:zb] - 
                                                currentInfo.x[:zg]' * x_interior[:zg] - 
                                                currentInfo.x[:zl]' * x_interior[:zl],  currentInfo.x] 
    elseif cutSelection == "LC"
        cutInfo = [ - currentInfo.f - currentInfo.x[:zb]' * ẑ[:zb] - 
                                        currentInfo.x[:zg]' * ẑ[:zg] - 
                                        currentInfo.x[:zl]' * ẑ[:zl],  currentInfo.x] 
    elseif cutSelection == "ShrinkageLC"
        cutInfo = [ currentInfo_f - currentInfo.x[:zb]' * ẑ[:zb] - 
                                        currentInfo.x[:zg]' * ẑ[:zg] - 
                                        currentInfo.x[:zl]' * ẑ[:zl],  
                                                                        currentInfo.x] 
    end 

    while true
        add_constraint(currentInfo, oracleInfo);
        optimize!(oracleModel);
        st = termination_status(oracleModel);
        if st == MOI.OPTIMAL
            f_star = JuMP.objective_value(oracleModel);
        else 
            return cutInfo
        end

        # formulate alpha model
        result = Δ_model_formulation(functionHistory, f_star, iter, Output = Output);
        previousΔ = copy.(Δ);
        Δ, a_min, a_max = result[1], result[2], result[3];

        if Output_Gap # && (iter % 30 == 0)
            if iter == 1
                println("------------------------------------ Iteration Info --------------------------------------")
                println("Iter |   Gap                              Objective                             Constraint")
            end
            @printf("%3d  |   %5.3g                         %5.3g                              %5.3g\n", iter, Δ, - currentInfo.f, currentInfo.G[1])
        end

        # push!(gap_list, Δ);
        x₀ = currentInfo.x;
        if round(previousΔ) > round(Δ)
            # x₀ = currentInfo.x;
            τₖ = μₖ * τₖ;
            if cutSelection == "ELC"
                cutInfo =  [ - currentInfo.f - currentInfo.x[:zb]' * x_interior[:zb] - 
                                                        currentInfo.x[:zg]' * x_interior[:zg] - 
                                                        currentInfo.x[:zl]' * x_interior[:zl],  currentInfo.x];
            elseif cutSelection == "LC"
                cutInfo = [ - currentInfo.f - currentInfo.x[:zb]' * ẑ[:zb] - 
                                                currentInfo.x[:zg]' * ẑ[:zg] - 
                                                currentInfo.x[:zl]' * ẑ[:zl],  currentInfo.x];
            elseif cutSelection == "ShrinkageLC"
                cutInfo = [ currentInfo_f - currentInfo.x[:zb]' * ẑ[:zb] - 
                                                currentInfo.x[:zg]' * ẑ[:zg] - 
                                                currentInfo.x[:zl]' * ẑ[:zl],  
                                                                                currentInfo.x] 
            end 
        else
            τₖ = (τₖ + τₘ) / 2;
        end

        # update α
        if μ/2 ≤ (α-a_min)/(a_max-a_min) .≤ 1-μ/2
            α = α;
        else
            α = round.((a_min+a_max)/2, digits = 6);
        end

        # update level
        w = α * f_star;
        W = minimum( α * functionHistory.f_his[j] + (1-α) * functionHistory.G_max_his[j] for j in 1:iter);

<<<<<<< HEAD
        λ = iter ≤ 10 ? 0.05 : 0.15;
        λ = iter ≥ 20 ? 0.25 : λ;
        λ = iter ≥ 30 ? 0.35 : λ;
        λ = iter ≥ 40 ? 0.45 : λ;
        λ = iter ≥ 50 ? 0.55 : λ;
        λ = iter ≥ 60 ? 0.65 : λ;
        λ = iter ≥ 70 ? 0.75 : λ;
        λ = iter ≥ 85 ? 0.85 : λ;
        λ = iter ≥ 90 ? 0.95 : λ;
=======
        λ = iter ≤ 10 ? 0.05 : 0.1;
        λ = iter ≥ 20 ? 0.15 : λ;
        λ = iter ≥ 30 ? 0.25 : λ;
        λ = iter ≥ 40 ? 0.35 : λ;
        λ = iter ≥ 50 ? 0.45 : λ;
        λ = iter ≥ 60 ? 0.55 : λ;
        λ = iter ≥ 70 ? 0.6 : λ;
        λ = iter ≥ 85 ? 0.7 : λ;
        λ = iter ≥ 90 ? 0.8 : λ;
>>>>>>> b0ba41fadf36400d3b50b9271b856ef7b4c90548
        
        level = round.(w + λ * (W - w), digits = 5)
        
        ## ==================================================== next iteration point ============================================== ##
        # obtain the next iteration point
        if iter == 1
            @constraint(nxtModel, levelConstraint, α * z1 + (1 - α) * y1 ≤ level);
        else 
            delete(nxtModel, nxtModel[:levelConstraint]);
            unregister(nxtModel, :levelConstraint);
            @constraint(nxtModel, levelConstraint, α * z1 + (1 - α) * y1 ≤ level);
        end
        add_constraint(currentInfo, nxtInfo);
        @objective(nxtModel, Min, sum((xb .- x₀[:zb]) .* (xb .- x₀[:zb])) +
                                   sum((xg .- x₀[:zg]) .* (xg .- x₀[:zg])) +
                                   sum((xl .- x₀[:zl]) .* (xl .- x₀[:zl])) #+ 2 * (α * z1 + (1 - α) * y1) * τₖ 
                                    );
        optimize!(nxtModel);
        st = termination_status(nxtModel)
        if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
            x_nxt = Dict{Symbol, Vector{Float64}}(:zb => round.(JuMP.value.(xb), digits = 5) , 
                                            :zg => round.(JuMP.value.(xg), digits = 5), 
                                            :zl => round.(JuMP.value.(xl), digits = 5)
                                            );
            λₖ = abs(dual(levelConstraint)); μₖ = λₖ + 1; 
        elseif st == MOI.NUMERICAL_ERROR ## need to figure out why this case happened and fix it
            # @info "Numerical Error occures! -- Build a new nxtModel"
            nxtModel = Model(
                optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                "OutputFlag" => Output, 
                "Threads" => 0)
                );
        
            @variable(nxtModel, xb[B]);
            @variable(nxtModel, xg[G]);
            @variable(nxtModel, xl[L]);
            @variable(nxtModel, z1);
            @variable(nxtModel, y1);

            nxtInfo = ModelInfo(nxtModel, xb, xg, xl, y1, z1);
            @constraint(nxtModel, levelConstraint, α * z1 + (1 - α) * y1 ≤ level);
            add_constraint(currentInfo, nxtInfo);
            @objective(nxtModel, Min, sum((xb .- x₀[:zb]) .* (xb .- x₀[:zb])) +
                                   sum((xg .- x₀[:zg]) .* (xg .- x₀[:zg])) +
                                   sum((xl .- x₀[:zl]) .* (xl .- x₀[:zl])) #+ 2 * (α * z1 + (1 - α) * y1) * τₖ 
                                    );
            optimize!(nxtModel);
            x_nxt = Dict{Symbol, Vector{Float64}}(:zb => round.(JuMP.value.(xb), digits = 5) , 
                                            :zg => round.(JuMP.value.(xg), digits = 5), 
                                            :zl => round.(JuMP.value.(xl), digits = 5)
                                            );
            λₖ = abs(dual(levelConstraint)); μₖ = λₖ + 1; 
        else
            # @info "Re-compute Next Iteration Point -- change to a safe level!"
            set_normalized_rhs( levelConstraint, w + .99 * (W - w))
            optimize!(nxtModel)
            x_nxt = Dict{Symbol, Vector{Float64}}(:zb => round.(JuMP.value.(xb), digits = 5) , 
                                            :zg => round.(JuMP.value.(xg), digits = 5), 
                                            :zl => round.(JuMP.value.(xl), digits = 5)
                                            );
        end

        ## stop rule: gap ≤ .07 * function-value && constraint ≤ 0.05 * LagrangianFunction
        if ( Δ ≤ threshold * 10 && currentInfo.G[1] ≤ threshold ) || (iter > max_iter) || (currentInfo.G[1] ≤ 0.0)
            return cutInfo
        end
        
        ## ==================================================== end ============================================== ##
        ## save the trajectory
        currentInfo, currentInfo_f = compute_f_G(x_nxt, cutSelection = cutSelection, f_star_value = f_star_value);
        iter = iter + 1;
        
        functionHistory.f_his[iter] = currentInfo.f;
        functionHistory.G_max_his[iter] = maximum(currentInfo.G[k] for k in keys(currentInfo.G));
    end

end




# function LevelSetMethod_Shrinkage!(  indexSets::IndexSets, 
#                                         paramDemand::ParamDemand, 
#                                         paramOPF::ParamOPF, 
#                                         ẑ::Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}, f_star_value::Float64,
#                                         randomVariables::RandomVariables;                          ## realization of the random time
#                                         levelSetMethodParam::LevelSetMethodParam = levelSetMethodParam, 
#                                         ϵ::Float64 = 1e-4
#                                         )

#     ## ==================================================== auxiliary function for function information ==================================================== ##
#     # μ larger is better
#     (μ, λ, threshold, nxt_bound, max_iter, Output, Output_Gap) = (levelSetMethodParam.μ, levelSetMethodParam.λ, levelSetMethodParam.threshold, levelSetMethodParam.nxt_bound, levelSetMethodParam.max_iter, levelSetMethodParam.Output,levelSetMethodParam.Output_Gap);
#     (D, G, L, B, T, Ω) = (indexSets.D, indexSets.G, indexSets.L, indexSets.B, indexSets.T, indexSets.Ω);
#     (Dᵢ, Gᵢ, in_L, out_L) = (indexSets.Dᵢ, indexSets.Gᵢ, indexSets.in_L, indexSets.out_L);

#     backwardInfo = backward_stage2_optimize!(indexSets, 
#                                     paramDemand, 
#                                     paramOPF, 
#                                     ẑ,
#                                     randomVariables                          ## realization of the random time
#                                     );


#     # collect the information from the objecive f, and constraints G
#     function compute_f_G(   x₀::Dict{Symbol, Vector{Float64}}; 
#                             f_star_value::Float64 = f_star_value, 
#                             indexSets::IndexSets = indexSets, 
#                             paramDemand::ParamDemand = paramDemand, 
#                             paramOPF::ParamOPF = paramOPF, randomVariables::RandomVariables = randomVariables, 
#                             ẑ::Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}} = ẑ, backwardInfo::BackwardInfo = backwardInfo
#                             )

#             @objective(backwardInfo.model, Min,  
#                                         sum( sum(paramDemand.w[d] * (1 - backwardInfo.x[d, t]) for d in indexSets.D ) for t in randomVariables.τ:indexSets.T) +
#                                         sum(paramDemand.cb[i] * backwardInfo.νb[i] for i in indexSets.B) + 
#                                         sum(paramDemand.cg[g] * backwardInfo.νg[g] for g in indexSets.G) + 
#                                         sum(paramDemand.cl[l] * backwardInfo.νl[l] for l in indexSets.L) +
#                                         x₀[:zb]' * (ẑ[:zb] .- backwardInfo.zb) + x₀[:zg]' * (ẑ[:zg] .- backwardInfo.zg) + x₀[:zl]' * (ẑ[:zl] .- backwardInfo.zl) 
#                                         )

#             ## ==================================================== solve the model and display the result ==================================================== ##
#             optimize!(backwardInfo.model)

#             F  = JuMP.objective_value(backwardInfo.Q);
#             negative_∇F = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 1}}(:zb => - ẑ[:zb] .+ round.(JuMP.value.(backwardInfo.zb)),
#                                                                                 :zg => - ẑ[:zg] .+ round.(JuMP.value.(backwardInfo.zg)),
#                                                                                 :zl => - ẑ[:zl] .+ round.(JuMP.value.(backwardInfo.zl))
#                                                                                 );

#             currentInfo = CurrentInfo( x₀,                                                                                                ## current point
#                                         1/2 * sum(x₀[:zb] .* x₀[:zb]) + 1/2 * sum(x₀[:zg] .* x₀[:zg]) + 1/2 * sum(x₀[:zl] .* x₀[:zl]) ,   ## obj function value
#                                         Dict(1 => (1 - ϵ) * f_star_value - F),                                                            ## constraint value
#                                         x₀,                                                                                               ## obj gradient
#                                         Dict(1 => negative_∇F )                                                                           ## constraint gradient
#                                         );
#         return (currentInfo = currentInfo, currentInfo_f = F)
#     end    
    
#     ## ==================================================== Levelset Method ============================================== ##
#     x₀ = Dict{Symbol, Vector{Float64}}(     :zb => ẑ[:zb] * 0, 
#                                             :zg => ẑ[:zg] * 0, 
#                                             :zl => ẑ[:zl] * 0
#                                             );
    
#     iter = 1
#     α = 1/2

#     # trajectory
#     currentInfo, currentInfo_f = compute_f_G(x₀);

#     functionHistory = FunctionHistory(  Dict(1 => currentInfo.f), 
#                                         Dict(1 => maximum(currentInfo.G[k] for k in keys(currentInfo.G)) )
#                                         );

#     # model for oracle
#     oracleModel = Model(
#         optimizer_with_attributes(
#             ()->Gurobi.Optimizer(GRB_ENV), 
#             "OutputFlag" => Output, 
#             "Threads" => 0)
#             );

#     ## ==================================================== Levelset Method ============================================== ##
#     para_oracle_bound = abs(currentInfo.f);
#     z_rhs = 2 * 10^(ceil(log10(para_oracle_bound)));
#     @variable(oracleModel, z ≥ - nxt_bound);
#     @variable(oracleModel, xb_oracle[B]);
#     @variable(oracleModel, xg_oracle[G]);
#     @variable(oracleModel, xl_oracle[L]);
#     @variable(oracleModel, y ≤ 0);

#     @objective(oracleModel, Min, z);
#     oracleInfo = ModelInfo(oracleModel, xb_oracle, xg_oracle, xl_oracle, y, z);


#     nxtModel = Model(
#         optimizer_with_attributes(
#         ()->Gurobi.Optimizer(GRB_ENV), 
#         "OutputFlag" => Output, 
#         "Threads" => 0)
#         )

#     @variable(nxtModel, xb[B]);
#     @variable(nxtModel, xg[G]);
#     @variable(nxtModel, xl[L]);
#     @variable(nxtModel, z1);
#     @variable(nxtModel, y1);
#     nxtInfo = ModelInfo(nxtModel, xb, xg, xl, y1, z1);

#     Δ = Inf; τₖ = 1; τₘ = .5; μₖ = 1;

#     cutInfo = [ currentInfo_f - currentInfo.x[:zb]' * ẑ[:zb] - 
#                                         currentInfo.x[:zg]' * ẑ[:zg] - 
#                                         currentInfo.x[:zl]' * ẑ[:zl],  
#                                                                         currentInfo.x] 

#     while true
#         add_constraint(currentInfo, oracleInfo);
#         optimize!(oracleModel);
#         f_star = JuMP.objective_value(oracleModel);

#         # formulate alpha model
#         result = Δ_model_formulation(functionHistory, f_star, iter, Output = Output);
#         previousΔ = Δ
#         Δ, a_min, a_max = result[1], result[2], result[3];

#         if Output_Gap # && (iter % 30 == 0)
#             if iter == 1
#                 println("------------------------------------ Iteration Info --------------------------------------")
#                 println("Iter |   Gap                              Objective                             Constraint")
#             end
#             @printf("%3d  |   %5.3g                         %5.3g                              %5.3g\n", iter, Δ,  currentInfo.f ,currentInfo.G[1])
#         end

#         # push!(gap_list, Δ);

#         if round(previousΔ) > round(Δ)
#             x₀ = currentInfo.x; τₖ = μₖ * τₖ;
#             cutInfo = [ currentInfo_f - currentInfo.x[:zb]' * ẑ[:zb] - 
#                                                 currentInfo.x[:zg]' * ẑ[:zg] - 
#                                                 currentInfo.x[:zl]' * ẑ[:zl],  currentInfo.x];
#         else
#             τₖ = (τₖ + τₘ) / 2;
#         end

#         # update α
#         if μ/2 ≤ (α-a_min)/(a_max-a_min) .≤ 1-μ/2
#             α = α;
#         else
#             α = (a_min+a_max)/2;
#         end

#         # update level
#         w = α * f_star;
#         W = minimum( α * functionHistory.f_his[j] + (1-α) * functionHistory.G_max_his[j] for j in 1:iter);

#         λ = iter ≤ 10 ? 0.05 : 0.15
#         λ = iter ≥ 20 ? 0.25 : λ
#         λ = iter ≥ 30 ? 0.4 : λ
#         λ = iter ≥ 40 ? 0.6 : λ
#         λ = iter ≥ 50 ? 0.7 : λ
#         λ = iter ≥ 55 ? 0.8 : λ
        
#         level = w + λ * (W - w);
        
#         ## ==================================================== next iteration point ============================================== ##
#         # obtain the next iteration point
#         if iter == 1
#             @constraint(nxtModel, levelConstraint, α * z1 + (1 - α) * y1 ≤ level);
#         else 
#             delete(nxtModel, nxtModel[:levelConstraint]);
#             unregister(nxtModel, :levelConstraint);
#             @constraint(nxtModel, levelConstraint, α * z1 + (1 - α) * y1 ≤ level);
#         end
#         add_constraint(currentInfo, nxtInfo);
#         @objective(nxtModel, Min, sum((xb .- x₀[:zb]) .* (xb .- x₀[:zb])) +
#                                    sum((xg .- x₀[:zg]) .* (xg .- x₀[:zg])) +
#                                    sum((xl .- x₀[:zl]) .* (xl .- x₀[:zl])) #+ 2 * (α * z1 + (1 - α) * y1) * τₖ 
#                                     );
#         optimize!(nxtModel);
#         st = termination_status(nxtModel)
#         if st == MOI.OPTIMAL || st == MOI.LOCALLY_SOLVED   ## local solution
#             x_nxt = Dict{Symbol, Vector{Float64}}(:zb => JuMP.value.(xb) , 
#                                             :zg => JuMP.value.(xg), 
#                                             :zl => JuMP.value.(xl)
#                                             );
#             λₖ = abs(dual(levelConstraint)); μₖ = λₖ + 1; 
#         elseif st == MOI.NUMERICAL_ERROR ## need to figure out why this case happened and fix it
#             @info "Numerical Error occures! -- Build a new nxtModel"

#             nxtModel = Model(
#                 optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
#                 "OutputFlag" => Output, 
#                 "Threads" => 0)
#                 )
        
#             @variable(nxtModel, xb[B]);
#             @variable(nxtModel, xg[G]);
#             @variable(nxtModel, xl[L]);
#             @variable(nxtModel, z1);
#             @variable(nxtModel, y1);

#             nxtInfo = ModelInfo(nxtModel, xb, xg, xl, y1, z1);
#             @constraint(nxtModel, levelConstraint, α * z1 + (1 - α) * y1 ≤ level);
#             add_constraint(currentInfo, nxtInfo);
#              @objective(nxtModel, Min, sum((xb .- x₀[:zb]) .* (xb .- x₀[:zb])) +
#                                    sum((xg .- x₀[:zg]) .* (xg .- x₀[:zg])) +
#                                    sum((xl .- x₀[:zl]) .* (xl .- x₀[:zl])) #+ 2 * (α * z1 + (1 - α) * y1) * τₖ 
#                                     );
#             optimize!(nxtModel);
#             x_nxt = Dict{Symbol, Vector{Float64}}(:zb => JuMP.value.(xb) , 
#                                                 :zg => JuMP.value.(xg), 
#                                                 :zl => JuMP.value.(xl)
#                                                 );
#         else
#             # @info "Re-compute Next Iteration Point -- change to a safe level!"
#             set_normalized_rhs( levelConstraint, w + .99 * (W - w))
#             optimize!(nxtModel)
#             x_nxt = Dict{Symbol, Vector{Float64}}(:zb => JuMP.value.(xb) , 
#                                             :zg => JuMP.value.(xg), 
#                                             :zl => JuMP.value.(xl)
#                                             );   
#         end

#         ## stop rule
#         if ( Δ < threshold && currentInfo.G[1] ≤ threshold / 2 ) || iter > max_iter
#             return cutInfo
#         end
        
#         ## ==================================================== end ============================================== ##
#         ## save the trajectory
#         currentInfo, currentInfo_f = compute_f_G(x_nxt);
#         iter = iter + 1;
        
#         functionHistory.f_his[iter] = currentInfo.f;
#         functionHistory.G_max_his[iter] = maximum(currentInfo.G[k] for k in keys(currentInfo.G));
#     end

# end
