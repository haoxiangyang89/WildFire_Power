using JuMP, Gurobi, PowerModels
using Statistics, StatsBase, Random, Dates
using Distributed, ParallelDataTransfer
using Agents
using InteractiveDynamics
using CairoMakie
using Distributions
using Geodesy
using CSV, DataFrames


const GRB_ENV = Gurobi.Env()


include("src/Multi_period_Wildfire/data_struct.jl")
include("src/Multi_period_Wildfire/backward_pass.jl")
include("src/Multi_period_Wildfire/forward_pass.jl")
include("src/Multi_period_Wildfire/gurobiTest.jl")

include("src/Multi_period_Wildfire/runtests_RTS_GMLC.jl")
# include("src/Multi_period_Wildfire/runtests_case30.jl") 

#############################################################################################
####################################    main function   #####################################
#############################################################################################

max_iter = 200; ϵ = 1e-2; Enhanced_Cut = true;

λ_value = .8; Output = 0; Output_Gap = false; Adj = false; Enhanced_Cut = true; threshold = 1e2; 
levelSetMethodParam = LevelSetMethodParam(0.95, λ_value, threshold, 1e14, 3e3, Output, Output_Gap, Adj)




function SDDiP_algorithm(Ω_rv::Dict{Int64, RandomVariables}, 
                            prob::Dict{Int64,Float64}, 
                            indexSets::IndexSets, 
                            paramDemand::ParamDemand, 
                            paramOPF::ParamOPF; 
                            levelSetMethodParam::LevelSetMethodParam = levelSetMethodParam,
                            ϵ::Float64 = 0.001, M::Int64 = 30, max_iter::Int64 = 200, 
                            Enhanced_Cut::Bool = true)
    ## d: x dim
    ## M: num of scenarios when doing one iteration
    initial = now();
    T = 2;
    
    i = 1;
    LB = - Inf;
    UB = Inf;
    
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
    col_types = [Int64, Float64, Float64, Float64, String, Float64, Float64];
    named_tuple = (; zip(col_names, type[] for type in col_types )...);
    sddipResult = DataFrame(named_tuple); # 0×7 DataFrame
    gapList = [];
    gurobiResult = gurobiOptimize!(indexSets, 
                                    paramDemand, 
                                    paramOPF, 
                                    Ω_rv,
                                    prob);  
    OPT = gurobiResult.OPT;
    
    forwardInfo = forward_stage1_model!(indexSets, 
                                                    paramDemand, 
                                                    paramOPF, 
                                                    Ω_rv,
                                                    prob,
                                                    cut_collection;  ## the index is ω
                                                    θ_bound = 0.0);
    forward2Info_List = Dict{Int64, Forward2Info}()
    for ω in indexSets.Ω
        forward2Info_List[ω] = forward_stage2_model!(indexSets, 
                                                    paramDemand,
                                                    paramOPF,
                                                    Ω_rv[ω]                        ## realization of the random time
                                                    )
    end 
    println("---------------- print out iteration information -------------------")
    while true
        t0 = now()
        M = 1  ## since we will enumerate all of realizations, hence, we only need to set M = 1
        Stage1_collection = Dict();  # to store every iteration results
        Stage2_collection = Dict{Int64, Float64}();  # to store every iteration results
        u = Vector{Float64}(undef, M);  # to compute upper bound

        ####################################################### Forwad Steps ###########################################################
        for k in 1:M
            ## stage 1
            optimize!(forwardInfo.model)
            state_variable = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 2}}(:zg => round.(JuMP.value.(forwardInfo.zg)), 
                                                                                        :zb => round.(JuMP.value.(forwardInfo.zb)), 
                                                                                        :zl => round.(JuMP.value.(forwardInfo.zl)))
            state_value    = JuMP.objective_value(forwardInfo.model) - sum(prob[ω] * JuMP.value(forwardInfo.θ[ω]) for ω in indexSets.Ω)       ## 1a first term
 
            Stage1_collection[k] = (state_variable = state_variable, 
                                    state_value = state_value, 
                                    obj_value = JuMP.objective_value(forwardInfo.model))  ## returen [state_variable, first_stage value, objective_value(Q)]
            LB = Stage1_collection[k].obj_value;
            if i > 1
                gap = round((OPT-LB)/OPT * 100 ,digits = 2);
                gapString = string(gap,"%");
                push!(sddipResult, [i, LB, OPT, UB, gapString, iter_time, total_Time]); push!(gapList, gap);

                @info "iter num is $(i-1), LB is $LB, OPT is $OPT UB is $UB"
                if OPT-LB <= ϵ * OPT || i > max_iter
                    # println(Stage1_collection[k].state_variable[:zg])
                    # println(gurobiResult.first_state_variable[:zg])
                    return Dict(:solHistory => sddipResult, :solution => Stage1_collection[k], :gapHistory => gapList) 
                end
            end


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
                                                    )

                ####################################################### solve the model and display the result ###########################################################
                optimize!(forward2Info_List[ω].model)
                state_obj_value    = JuMP.objective_value(forward2Info_List[ω].model)
                Stage2_collection[ω] = state_obj_value
                c = c + prob[ω] * state_obj_value;
            end
            u[k] = Stage1_collection[k].state_value + c;
        end

        ## compute the upper bound
        UB = mean(u)

        ####################################################### Backward Steps ###########################################################

        for k in 1:M 
            for ω in keys(Ω_rv)
                # @info "$t $k $j"
                randomVariables = Ω_rv[ω]
                ẑ = Dict(   :zg => Stage1_collection[k][1][:zg][:, randomVariables.τ - 1], 
                            :zb => Stage1_collection[k][1][:zb][:, randomVariables.τ - 1], 
                            :zl => Stage1_collection[k][1][:zl][:, randomVariables.τ - 1]
                            )
                coef = LevelSetMethod_optimization!(indexSets, paramDemand, paramOPF, 
                                                                    ẑ,  
                                                                    Stage2_collection[ω], randomVariables,                 
                                                                    levelSetMethodParam = levelSetMethodParam, 
                                                                    ϵ = 1e-3, 
                                                                    interior_value = 0.5, 
                                                                    Enhanced_Cut = true
                                                                    )
                # add cut
                if i ≥ 3 
                    cut_collection[ω].v[i] = Dict{Int64, Float64}()
                    cut_collection[ω].πb[i] = Dict{Int64, Vector{Float64}}()
                    cut_collection[ω].πg[i] = Dict{Int64, Vector{Float64}}()
                    cut_collection[ω].πl[i] = Dict{Int64, Vector{Float64}}()
                end
                cut_collection[ω].v[i][k] = coef[1]
                cut_collection[ω].πb[i][k] = coef[2][:zb]
                cut_collection[ω].πg[i][k] = coef[2][:zg]
                cut_collection[ω].πl[i][k] = coef[2][:zl]
            end
        end
        
        ########################################################################################################
        # add cut for value functions
        for ω in indexSets.Ω
            cut_coefficient = cut_collection[ω]
            ωk = length(keys(cut_coefficient.v[1]))  ## scenario num
            τ = Ω_rv[ω].τ
            @constraint(forwardInfo.model, [m in 1:ωk], forwardInfo.θ[ω] .>= cut_coefficient.v[i][m] + 
                                                        cut_coefficient.πb[i][m]' * forwardInfo.zb[:, τ-1] +
                                                        cut_coefficient.πg[i][m]' * forwardInfo.zg[:, τ-1] +
                                                        cut_coefficient.πl[i][m]' * forwardInfo.zl[:, τ-1]
                                                        )
        end

        t1 = now();
        iter_time = (t1 - t0).value/1000;
        total_Time = (t1 - initial).value/1000;
        
        i = i + 1;

    end

end


λ_value = .5; Output = 0; Output_Gap = true; Adj = false; Enhanced_Cut = true; threshold = 1e2; 
levelSetMethodParam = LevelSetMethodParam(0.95, λ_value, threshold, 1e14, 3e3, Output, Output_Gap, Adj)
LevelSetMethod_optimization!(indexSets, paramDemand, paramOPF, 
                                                                    ẑ,  
                                                                    Stage2_collection[ω], randomVariables,                 
                                                                    levelSetMethodParam = levelSetMethodParam, 
                                                                    ϵ = 1e-3, 
                                                                    interior_value = 0.5, 
                                                                    Enhanced_Cut = true
                                                                    )