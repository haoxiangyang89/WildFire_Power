
# ----------------------------------------------------------------------------------------------------------------------- #
# --------------------------------------------- Decomposition Algorithm ------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------------------- #
@everywhere begin 
    function backwardParalllel(f_star_value::Float64, randomVariables::RandomVariables, ẑ;
                                indexSets::IndexSets = indexSets, 
                                paramDemand::ParamDemand = paramDemand, 
                                paramOPF::ParamOPF = paramOPF)
        cutSelection = "ShrinkageLC"
        (x_interior, levelSetMethodParam, x₀) = setupLevelSetMethod(ẑ, f_star_value; cutSelection = cutSelection, 
                                                                                Output_Gap = true, 
                                                                                    ℓ1 = 0., # or 0   ## adjust x0
                                                                                        ℓ2 = .5, ## adjust interior points
                                                                                        λ = nothing);  
        coef = LevelSetMethod_optimization!(ẑ, f_star_value, randomVariables;
                                                levelSetMethodParam = levelSetMethodParam,
                                                    cutSelection = cutSelection,            ## "ELC", "LC", "ShrinkageLC" 
                                                        x_interior = x_interior, 
                                                            x₀ = x₀, tightness = true
                                            )
        return coef
    end
end

@everywhere function setupLevelSetMethod(ẑ, f_star_value::Float64; 
                                    cutSelection::String = cutSelection, 
                                        Output_Gap::Bool = false,  
                                            ℓ1::Real = 2, ℓ2::Real = .8, λ::Union{Real, Nothing} = .1  
                            )
    if cutSelection == "ShrinkageLC"
        λ_value = λ; Output = 0; threshold = 1e-3 * f_star_value; 
        levelSetMethodParam = LevelSetMethodParam(0.9, λ_value, threshold, 1e13, 50, Output, Output_Gap);
        x_interior = nothing;
    elseif cutSelection == "ELC"
        λ_value = λ; Output = 0; threshold = 5e-3 * f_star_value; 
        levelSetMethodParam = LevelSetMethodParam(0.9, λ_value, threshold, 1e13, 100, Output, Output_Gap);
        x_interior = Dict{Symbol, Vector{Float64}}(:zg => ẑ[:zg] .* ℓ2  .+ (1 - ℓ2)/2, 
                                :zb => ẑ[:zb] .* ℓ2  .+ (1 - ℓ2)/2, 
                                        :zl => ẑ[:zl] .* ℓ2  .+ (1 - ℓ2)/2
                                        );
    elseif cutSelection == "LC"
        λ_value = λ; Output = 0; threshold = 1e-5 * f_star_value; 
        levelSetMethodParam = LevelSetMethodParam(0.95, λ_value, threshold, 1e15, 30, Output, Output_Gap);
        x_interior = nothing;
    end

    x₀ = Dict{Symbol, Vector{Float64}}(  :zb => round.(ẑ[:zb] * ℓ1 * f_star_value .- f_star_value * (ℓ1 / 2), digits = 2)  , 
                                                :zg =>  round.(ẑ[:zg] * ℓ1 * f_star_value .- f_star_value * (ℓ1 / 2), digits = 2)  , 
                                                    :zl => round.(ẑ[:zl] * ℓ1 * f_star_value .- f_star_value * (ℓ1 / 2), digits = 2)  
                                                );

    return (x_interior = x_interior, 
            levelSetMethodParam = levelSetMethodParam,
            x₀ = x₀
            )
end



initial = now(); i = 1; LB = - Inf; UB = Inf;
iter_time = 0; total_Time = 0; t0 = 0.0;

col_names = [:iter, :LB, :OPT, :UB, :gap, :time, :Time]; # needs to be a vector Symbols
col_types = [Int64, Float64, Union{Float64,Nothing}, Float64, String, Float64, Float64];
named_tuple = (; zip(col_names, type[] for type in col_types )...);
sddipResult = DataFrame(named_tuple); # 0×7 DataFrame
gapList = [];

# forward Master problem 
forwardInfo = forward_stage1_model!(indexSets, 
                                            paramDemand, 
                                            paramOPF, 
                                            Ω_rv,
                                            prob;
                                            θ_bound = 0.0);

while true
    t0 = now();                                     ## since we will enumerate all of realizations, hence, we only need to set M = 1
    Stage1_collection = Dict();                     # to store every iteration results
    Stage2_collection = [];     # to store every iteration results

    ####################################################### Forwad Steps ###########################################################
    # stage 1
    optimize!(forwardInfo.model);
    state_variable = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 2}}(:zg => round.(JuMP.value.(forwardInfo.zg)), 
                                                                                :zb => round.(JuMP.value.(forwardInfo.zb)), 
                                                                                :zl => round.(JuMP.value.(forwardInfo.zl)));
    state_value    = JuMP.objective_value(forwardInfo.model) - sum(prob[ω] * JuMP.value(forwardInfo.θ[ω]) for ω in indexSets.Ω);       ## 1a first term

    Stage1_collection[1] = (state_variable = state_variable, 
                            state_value = state_value, 
                            obj_value = JuMP.objective_value(forwardInfo.model));  ## returen [state_variable, first_stage value, objective_value(Q)]
    LB = round(Stage1_collection[1].obj_value, digits = 5);

    ẑList = Dict(ω => Dict( :zg => Stage1_collection[1][1][:zg][:, Ω_rv[ω].τ - 1], 
                            :zb => Stage1_collection[1][1][:zb][:, Ω_rv[ω].τ - 1], 
                            :zl => Stage1_collection[1][1][:zl][:, Ω_rv[ω].τ - 1]) 
                                for ω in indexSets.Ω
                        );
    ## stage 2
    Stage2_collection = pmap(forward_stage2_optimize!, [ẑList[ω] for ω in indexSets.Ω], [Ω_rv[ω] for ω in indexSets.Ω])

    ## compute the upper bound
    UB = minimum([Stage1_collection[1].state_value + sum(prob[ω] * Stage2_collection[ω] for ω in indexSets.Ω), UB]); 
    gap = round((UB-LB)/UB * 100 ,digits = 2); gapString = string(gap,"%");
    push!(sddipResult, [i, LB, OPT, UB, gapString, iter_time, total_Time]); push!(gapList, gap);
    if i == 1
        println("---------------------------------- Iteration Info ------------------------------------")
        println("Iter |   LB                              UB                             gap")
    end
    @printf("%3d  |   %5.3g                         %5.3g                              %1.3f%s\n", i, LB, UB, gap, "%")
    if UB-LB ≤ 1e-2 * UB || i > max_iter
        # Stage1_collection[1].state_variable[:zl] == gurobiResult.first_state_variable[:zl]
        return Dict(:solHistory => sddipResult, 
                        :solution => Stage1_collection[1], 
                        :gapHistory => gapList
                        ) 
    end

    ####################################################### Backward Steps ###########################################################
    cutCollection = pmap(backwardParalllel, 
                                [Stage2_collection[ω] for ω in indexSets.Ω], 
                                [Ω_rv[ω] for ω in indexSets.Ω], 
                                [ẑList[ω] for ω in indexSets.Ω]
                                );
    # add adds 
    for ω in indexSets.Ω
        τ = Ω_rv[ω].τ;
        @constraint(forwardInfo.model, forwardInfo.θ[ω] .≥ cutCollection[ω][1] + 
                                                            cutCollection[ω][2][:zb]' * forwardInfo.zb[:, τ-1] +
                                                            cutCollection[ω][2][:zg]' * forwardInfo.zg[:, τ-1] +
                                                            cutCollection[ω][2][:zl]' * forwardInfo.zl[:, τ-1]
                                                    );
    end

    t1 = now(); iter_time = (t1 - t0).value/1000; total_Time = (t1 - initial).value/1000; i = i + 1;

end
