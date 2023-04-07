#############################################################################################
####################################    Without Parallel   ##################################
#############################################################################################
using JuMP, Gurobi, PowerModels
using Statistics, StatsBase, Random, Dates, Distributions
using Distributed, ParallelDataTransfer
# using Agents
# using InteractiveDynamics
# using CairoMakie
# using Geodesy
using CSV, DataFrames, Printf
using JLD2, FileIO

const GRB_ENV = Gurobi.Env()


include("src/alg/def.jl")
include("src/alg/backwardPass.jl")
include("src/alg/forwardPass.jl")
include("src/alg/extFormGurobi.jl")
include("src/alg/sddip.jl")

# ## generate data, simulation
# include("src/alg/CellularAutomaton.jl")
# include("src/alg/readin.jl")
# include("src/alg/runtests_RTS_GMLC.jl")

## RTS test
indexSets = load("data/testData_RTS/indexSets.jld2")["indexSets"]
paramOPF = load("data/testData_RTS/paramOPF.jld2")["paramOPF"]
paramDemand = load("data/testData_RTS/paramDemand.jld2")["paramDemand"]
Ω_rv = load("data/testData_RTS/Ω_rv.jld2")["Ω_rv"]
prob = load("data/testData_RTS/prob.jld2")["prob"]

sddipResult = SDDiP_algorithm(; ϵ = 1e-4, max_iter = 100) 
@save "testData_RTS/sampleSize2/enhancedRTS.jld2" sddipResult
@load "testData_RTS/sampleSize2/enhancedRTS.jld2" sddipResult

# using Latexify
# latexify(df; env=:table, latex=false)

# for i in [1,2,3,4,6,7]
# sddipResult[:,i] = round.(sddipResult[:,i], digits = 1)
# end
# df = sddipResult[:,[2,4,5,7]]
# latexify(df; env=:table, latex=false)