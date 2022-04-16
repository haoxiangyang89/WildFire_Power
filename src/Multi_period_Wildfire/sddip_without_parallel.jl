using JuMP, Test, Statistics, StatsBase, Gurobi, Distributed, Distributions, ParallelDataTransfer, Random, DataFrames, Dates, PowerModels

const GRB_ENV = Gurobi.Env()


include("data_struct.jl")
include("backward_pass.jl")
include("forward_pass.jl")
include("gurobiTest.jl")
include("runtests_small3.jl")  ## M = 4

#############################################################################################
####################################    main function   #####################################
#############################################################################################
max_iter = 200; ϵ = 1e-2; Enhanced_Cut = true


function SDDiP_algorithm(Ω_rv::Dict{Int64,Dict{Int64,RandomVariables}}, 
                            prob::Dict{Int64,Float64}, 
                            indexSets::IndexSets; 
                            scenario_sequence::Dict{Int64, Dict{Int64, Any}} = scenario_sequence, 
                            ϵ::Float64 = 0.001, M::Int64 = 30, max_iter::Int64 = 200, 
                            Enhanced_Cut::Bool = true, binaryInfo::BinaryInfo = binaryInfo)
    ## d: x dim
    ## M: num of scenarios when doing one iteration
    initial = now()
    T = 2
    
    i = 1
    LB = - Inf 
    UB = Inf
    
    cut_collection = Dict{Int64, CutCoefficient}()  # here, the index is ω

    for ω in indexSets.Ω

        cut_collection[ω] = CutCoefficient(
                                            Dict(1=> Dict(1=> 0.0), 2 => Dict()),   ## v
                                            Dict(1=> Dict(1=> zeros(Float64, length(indexSets.B))), 2 => Dict()),  ## πb
                                            Dict(1=> Dict(1=> zeros(Float64, length(indexSets.G))), 2 => Dict()),  ## πg
                                            Dict(1=> Dict(1=> zeros(Float64, length(indexSets.L))), 2 => Dict()),  ## πl
                                          )
    end

    col_names = [:iter, :LB, :OPT, :UB, :gap, :time, :Time] # needs to be a vector Symbols
    col_types = [Int64, Float64, Float64, Float64, String, Float64, Float64]
    named_tuple = (; zip(col_names, type[] for type in col_types )...)
    sddipResult = DataFrame(named_tuple) # 0×7 DataFrame
    gapList = []
    # gurobiResult = gurobiOptimize!(Ω, prob, StageCoefficient,
    #                                 binaryInfo = binaryInfo)
    # # OPT = round!(gurobiResult[1])[3]
    # OPT = gurobiResult[1]
    println("---------------- print out iteration information -------------------")
    while true
        t0 = now()
        M = 1  ## since we will enumerate all of realizations, hence, we only need to set M = 1
        Stage1_collection = Dict();  # to store every iteration results
        Stage2_collection = Dict();  # to store every iteration results
        u = Vector{Float64}(undef, M);  # to compute upper bound
        
        ## Forward Step
        for k in 1:M
            ## stage 1
            Stage1_collection[k] = forward_stage1_optimize!(indexSets, 
                                                            paramDemand, 
                                                            paramOPF, 
                                                            Ω_rv,
                                                            prob,
                                                            cut_collection;  ## the index is ω
                                                            θ_bound = 0.0)
            
            ## stage 2
            first_stage_decision = Stage1_collection[k][1]
            c = 0.0
            for ω in indexSets.Ω
                randomVariables = Ω_rv[ω]

                ẑ = Dict(   :zg => first_stage_decision[:zg][:, randomVariables.τ - 1], 
                            :zb => first_stage_decision[:zb][:, randomVariables.τ - 1], 
                            :zl => first_stage_decision[:zl][:, randomVariables.τ - 1]
                            )

                Stage2_collection[ω, k] = forward_stage2_optimize!(indexSets, 
                                                                paramDemand,
                                                                paramOPF,
                                                                ẑ,
                                                                randomVariables                        ## realization of the random time
                                                                )[2]
                c = c + prob[ω] * Stage2_collection[ω, k]
            end
            u[k] = Stage1_collection[k][2] + c
        end

        ## compute the upper bound
        μ̄ = mean(u)
        # σ̂² = Distributions.var(u)
        # UB = μ̄ + 1.96 * sqrt(σ̂²/M)
        # UB = round!(UB)[3]
        ##################################### Parallel Computation for backward step ###########################

        for k in 1:M 
            for ω in Ω
                # @info "$t $k $j"
                ϵ_value = 1e-5 # 1e-5
                λ_value = .1; Output = 0; Output_Gap = false; Adj = false; Enhanced_Cut = true; threshold = 1e-5; 
                levelSetMethodParam = LevelSetMethodParam(0.95, λ_value, threshold, 1e14, 3e3, Output, Output_Gap, Adj)
                randomVariables = Ω_rv[ω]
                ẑ = Dict(   :zg => first_stage_decision[:zg][:, randomVariables.τ - 1], 
                            :zb => first_stage_decision[:zb][:, randomVariables.τ - 1], 
                            :zl => first_stage_decision[:zl][:, randomVariables.τ - 1]
                            )
                coef = LevelSetMethod_optimization!(indexSets, paramDemand, paramOPF, 
                                                                    ẑ, randomVariables,                 
                                                                    levelSetMethodParam = levelSetMethodParam, 
                                                                    ϵ = 1e-4, 
                                                                    interior_value = 0.5, 
                                                                    Enhanced_Cut = true
                                                                    )
                # add cut
                cut_collection[ω].v[i][k] = coef[1]
                cut_collection[ω].πb[i][k] = coef[2][:zb]
                cut_collection[ω].πg[i][k] = coef[2][:zg]
                cut_collection[ω].πl[i][k] = coef[2][:zl]
            end
        end
        
        ########################################################################################################


        ## compute the lower bound
        _LB = forward_stage1_optimize!(indexSets, 
                                        paramDemand, 
                                        paramOPF, 
                                        Ω_rv,
                                        prob,
                                        cut_collection;  ## the index is ω
                                        θ_bound = 0.0
                                        )
        # LB = round!(_LB[3] + _LB[4])[3]
        LB = _LB[3]
        t1 = now()
        iter_time = (t1 - t0).value/1000
        total_Time = (t1 - initial).value/1000
        gap = round((OPT-LB)/OPT * 100 ,digits = 2)
        gapString = string(gap,"%")
        push!(sddipResult, [i, LB, OPT, UB, gapString, iter_time, total_Time]); push!(gapList, 100-gap);
        
        
        i = i + 1
        
        @info "iter num is $(i-1), LB is $LB, UB is $UB, first stage decision is $(A*_LB[1])"
        if OPT-LB <= ϵ * OPT || i > max_iter
            return Dict(:solHistory => sddipResult, :solution => _LB, :gapHistory => gapList) 
        end

    end

end



using JLD2, FileIO, DataFrames
result_enhanced = copy(Dict(:solHistory => sddipResult, :solution => _LB, :gapHistory => gapList) )
cut_enhanced = copy(cut_collection)

@save "runtests_small2_enhanced.jld2" result_enhanced cut_enhanced
# # @load "runtests_small2_enhanced.jld2" result_enhanced cut_enhanced

# result_LC = copy(Dict(:solHistory => sddipResult, :solution => _LB, :gapHistory => gapList) )
# cut_LC = copy(cut_collection)

# @save "runtests_small2_LC.jld2" result_LC cut_LC
# # @load "runtests_small2_LC.jld2" result_LC cut_LC



# using DataFrames
# using Latexify
# df = DataFrame(A = 'x':'z', B = ["M", "F", "F"])
# latexify(df; env=:table, latex=false)




