#############################################################################################
#######################################   Parallel   #######################################
#############################################################################################
using Distributed
addprocs(4)

@everywhere begin
    using JuMP, Gurobi
    using Statistics, StatsBase, Random, Dates
    using Distributed, ParallelDataTransfer
    using Distributions
    using DataFrames, Printf

    # using Agents, Distributions
    # using CSV, Geodesy
    # using InteractiveDynamics
    # using CairoMakie
    # using PowerModels


    const GRB_ENV = Gurobi.Env()
    
    include("src/alg/def.jl")
    include("src/alg/backwardPass.jl")
    include("src/alg/forwardPass.jl")
    include("src/alg/extFormGurobi.jl")
end

# include("src/alg/CellularAutomaton.jl")
# include("src/alg/readin.jl")
# # include("src/alg/runtests_case30.jl")  
# include("src/alg/runtests_RTS_GMLC.jl") 
include("src/alg/sddipParallel.jl") 

using JLD2, FileIO
indexSets = load("src/testData/RTS_24_4/indexSets.jld2")["indexSets"]
paramOPF = load("src/testData/RTS_24_4/paramOPF.jld2")["paramOPF"]
paramDemand = load("src/testData/RTS_24_4/paramDemand.jld2")["paramDemand"]
Ω_rv = load("src/testData/RTS_24_4/Ω_rv.jld2")["Ω_rv"]
prob = load("src/testData/RTS_24_4/prob.jld2")["prob"]
@passobj 1 workers() indexSets
@passobj 1 workers() paramOPF
@passobj 1 workers() paramDemand
@passobj 1 workers() Ω_rv
@passobj 1 workers() prob


sddipResultParallel = SDDiP_algorithm(; ϵ = 1e-4, max_iter = 100) 
@save "sddipResultParallel.jld2" sddipResultParallel
@load "sddipResultParallel.jld2" sddipResultParallel



@save "cut_collection.jld2" cut_collection
@load "cut_collection.jld2" cut_collection


# paramOPF = load("src/testData/RTS_24_200/paramOPF.jld2")["paramOPF"]
# for l in indexSets.L 
#     paramOPF.b[l] = paramOPF.b[l] * 10
# end
# save("src/testData/RTS_24_4/paramOPF.jld2", "paramOPF", paramOPF)