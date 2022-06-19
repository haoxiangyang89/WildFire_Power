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



indexSets = load("indexSets.jld2")["indexSets"]
paramOPF = load("paramOPF.jld2")["paramOPF"]
paramDemand = load("paramDemand.jld2")["paramDemand"]
Ω_rv = load("Ω_rv.jld2")["Ω_rv"]
prob = load("prob.jld2")["prob"]

sddipResult_Enhanced = SDDiP_algorithm(; ϵ = 1e-4, max_iter = 100) 
@save "enhancedRTS.jld2" sddipResult_Enhanced
@load "enhancedRTS.jld2" sddipResult_Enhanced


@save "enhancedRTS.jld2" sddipResult
@load "enhancedRTS.jld2" sddipResult



