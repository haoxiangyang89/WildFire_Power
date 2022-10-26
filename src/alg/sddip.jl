
## ---------------------------------------------   Auxiliary Function   --------------------------------------------- ##
# setup the parameters for different cuts
# setup the interiors point and the initial points

function setupLevelSetMethod(ẑ, f_star_value::Float64; 
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

## ---------------------------------------------   Main Function   --------------------------------------------- ##

"""
    SDDiP_algorithm( ;
                            Ω_rv::Dict{Int64, RandomVariables} = Ω_rv, 
                            prob::Dict{Int64,Float64} = prob, 
                            indexSets::IndexSets = indexSets, 
                            paramDemand::ParamDemand = paramDemand, 
                            paramOPF::ParamOPF = paramOPF, max_iter::Int64 = 100)

TBW
"""
function SDDiP_algorithm( ; Ω_rv::Dict{Int64, RandomVariables} = Ω_rv, 
                            prob::Dict{Int64,Float64} = prob, 
                            indexSets::IndexSets = indexSets, 
                            paramDemand::ParamDemand = paramDemand, 
                            paramOPF::ParamOPF = paramOPF, max_iter::Int64 = 100)
    ## d: x dim
    ## M: num of scenarios when doing one iteration
    initial = now(); T = 2; i = 1; LB = - Inf; UB = Inf;
    iter_time = 0; total_Time = 0; t0 = 0.0;
    cut_collection = Dict{Int64, CutCoefficient}();  # here, the index is ω

    for ω in indexSets.Ω

        cut_collection[ω] = CutCoefficient(
                                            Dict(1=> Dict(1=> 0.0), 2 => Dict()),   ## v
                                            Dict(1=> Dict(1=> zeros(Float64, length(indexSets.B))), 2 => Dict()),  ## πb
                                            Dict(1=> Dict(1=> zeros(Float64, length(indexSets.G))), 2 => Dict()),  ## πg
                                            Dict(1=> Dict(1=> zeros(Float64, length(indexSets.L))), 2 => Dict()),  ## πl
                                          );
    end

    col_names = [:iter, :LB, :OPT, :UB, :gap, :time, :Time]; # needs to be a vector Symbols
    col_types = [Int64, Float64, Union{Float64,Nothing}, Float64, String, Float64, Float64];
    named_tuple = (; zip(col_names, type[] for type in col_types )...);
    sddipResult = DataFrame(named_tuple); # 0×7 DataFrame
    gapList = [];
    @time gurobiResult = gurobiOptimize!(indexSets, 
                                    paramDemand, 
                                    paramOPF, 
                                    Ω_rv,
                                    prob);  
    OPT = gurobiResult.OPT;
    # OPT = nothing;
    forwardInfo = forward_stage1_model!(indexSets, 
                                                    paramDemand, 
                                                    paramOPF, 
                                                    Ω_rv,
                                                    prob,
                                                    cut_collection;  ## the index is ω
                                                    θ_bound = 0.0);
    forward2Info_List = Dict{Int64, Forward2Info}();
    for ω in indexSets.Ω
        forward2Info_List[ω] = forward_stage2_model!(indexSets, 
                                                    paramDemand,
                                                    paramOPF,
                                                    Ω_rv[ω]                        ## realization of the random time
                                                    );
    end 

    while true
        t0 = now();
        M = 1;  ## since we will enumerate all of realizations, hence, we only need to set M = 1
        Stage1_collection = Dict();  # to store every iteration results
        Stage2_collection = Dict{Int64, Float64}();  # to store every iteration results
        u = Vector{Float64}(undef, M);  # to compute upper bound

        ####################################################### Forwad Steps ###########################################################
        for k in 1:M
            ## stage 1
            optimize!(forwardInfo.model);
            state_variable = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 2}}(:zg => round.(JuMP.value.(forwardInfo.zg)), 
                                                                                        :zb => round.(JuMP.value.(forwardInfo.zb)), 
                                                                                        :zl => round.(JuMP.value.(forwardInfo.zl)));
            state_value    = JuMP.objective_value(forwardInfo.model) - sum(prob[ω] * JuMP.value(forwardInfo.θ[ω]) for ω in indexSets.Ω);       ## 1a first term
 
            Stage1_collection[k] = (state_variable = state_variable, 
                                    state_value = state_value, 
                                    obj_value = JuMP.objective_value(forwardInfo.model));  ## returen [state_variable, first_stage value, objective_value(Q)]
            LB = round(Stage1_collection[k].obj_value, digits = 5);

            ## stage 2
            # first_stage_decision = Stage1_collection[k].state_variable
            c = 0.0
            for ω in indexSets.Ω
                randomVariables = Ω_rv[ω];

                ẑ = Dict(   :zg => Stage1_collection[k].state_variable[:zg][:, randomVariables.τ - 1], 
                            :zb => Stage1_collection[k].state_variable[:zb][:, randomVariables.τ - 1], 
                            :zl => Stage1_collection[k].state_variable[:zl][:, randomVariables.τ - 1]
                            );
                ## modify the constraints according to the first stage state variables
                forward_stage2_modify_constraints!(indexSets, 
                                                    forward2Info_List[ω],
                                                    i,
                                                    ẑ,
                                                    randomVariables
                                                    );

                ####################################################### solve the model and display the result ###########################################################
                optimize!(forward2Info_List[ω].model);
                state_obj_value    = JuMP.objective_value(forward2Info_List[ω].model);
                Stage2_collection[ω] = state_obj_value;
                c = c + prob[ω] * state_obj_value;
            end
            u[k] = Stage1_collection[k].state_value + c;
        end

        ## compute the upper bound
        # UB = minimum([mean(u), UB]); 
        UB = mean(u);
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
                            :gapHistory => gapList, 
                            :cutHistory => cut_collection) 
        end

        # sddipResult
        ####################################################### Backward Steps ###########################################################
        for k in 1:M 
            for ω in keys(Ω_rv)
                # @info "$i $ω"
                randomVariables = Ω_rv[ω];
                ẑ = Dict(   :zg => Stage1_collection[k].state_variable[:zg][:, randomVariables.τ - 1], 
                            :zb => Stage1_collection[k].state_variable[:zb][:, randomVariables.τ - 1], 
                            :zl => Stage1_collection[k].state_variable[:zl][:, randomVariables.τ - 1]
                            );
                f_star_value = Stage2_collection[ω];

                cutSelection = "ShrinkageLC";                                                               ## "ELC", "LC", "ShrinkageLC" 
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
                
                # add cut
                if i ≥ 3 
                    cut_collection[ω].v[i] = Dict{Int64, Float64}();
                    cut_collection[ω].πb[i] = Dict{Int64, Vector{Float64}}();
                    cut_collection[ω].πg[i] = Dict{Int64, Vector{Float64}}();
                    cut_collection[ω].πl[i] = Dict{Int64, Vector{Float64}}();
                end
                cut_collection[ω].v[i][k] = coef[1];
                cut_collection[ω].πb[i][k] = coef[2][:zb];
                cut_collection[ω].πg[i][k] = coef[2][:zg];
                cut_collection[ω].πl[i][k] = coef[2][:zl];
            end
        end
        
        ########################################################### add cut ###############################################################
        # add cut for value functions 
        for ω in indexSets.Ω
            cutCoefficient = cut_collection[ω];
            ωk = length(keys(cutCoefficient.v[1]));  ## scenario num
            τ = Ω_rv[ω].τ;
            @constraint(forwardInfo.model, [m in 1:ωk], forwardInfo.θ[ω] .≥ cutCoefficient.v[i][m] + 
                                                        cutCoefficient.πb[i][m]' * forwardInfo.zb[:, τ-1] +
                                                        cutCoefficient.πg[i][m]' * forwardInfo.zg[:, τ-1] +
                                                        cutCoefficient.πl[i][m]' * forwardInfo.zl[:, τ-1]
                                                        );
        end 
    
        t1 = now(); iter_time = (t1 - t0).value/1000; total_Time = (t1 - initial).value/1000; i = i + 1;

    end

end
