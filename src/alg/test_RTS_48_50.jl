using JuMP, Gurobi, PowerModels
using Statistics, StatsBase, Random, Dates, Distributions
using Distributed, ParallelDataTransfer
using CSV, DataFrames, Printf
using JLD2, FileIO

const GRB_ENV = Gurobi.Env()


include("src/alg/def.jl")
include("src/alg/backwardPass.jl")
include("src/alg/forwardPass.jl")
include("src/alg/extFormGurobi.jl")
include("src/alg/sddip.jl")



indexSets = load("src/test/RTS_48_50/indexSets.jld2")["indexSets"]
paramOPF = load("src/test/RTS_48_50/paramOPF.jld2")["paramOPF"]
paramDemand = load("src/test/RTS_48_50/paramDemand.jld2")["paramDemand"]
Ω_rv = load("src/test/RTS_48_50/Ω_rv.jld2")["Ω_rv"]
prob = load("src/test/RTS_48_50/prob.jld2")["prob"]

sddipResult_Enhanced = SDDiP_algorithm(; ϵ = 1e-4, max_iter = 100) 
@save "src/test/RTS_48_50/enhancedRTS.jld2" sddipResult_Enhanced
@load "src/test/RTS_48_50/enhancedRTS.jld2" sddipResult_Enhanced


@save "src/test/RTS_48_50/enhancedRTS.jld2" sddipResult
@load "src/test/RTS_48_50/enhancedRTS.jld2" sddipResult

# using Latexify
# latexify(df; env=:table, latex=false)
