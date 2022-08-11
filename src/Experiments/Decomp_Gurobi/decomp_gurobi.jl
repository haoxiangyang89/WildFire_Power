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

## -------------------------------------- Scenario = 20

indexSets = load("src/testData/RTS_24_20/indexSets.jld2")["indexSets"]
paramOPF = load("src/testData/RTS_24_20/paramOPF.jld2")["paramOPF"]
paramDemand = load("src/testData/RTS_24_20/paramDemand.jld2")["paramDemand"]
Ω_rv = load("src/testData/RTS_24_20/Ω_rv.jld2")["Ω_rv"]
prob = load("src/testData/RTS_24_20/prob.jld2")["prob"]

sddipResult_Enhanced = SDDiP_algorithm(; ϵ = 1e-4, max_iter = 100) 
@time gurobiResult = gurobiOptimize!(indexSets, 
                                    paramDemand, 
                                    paramOPF, 
                                    Ω_rv,
                                    prob; timelimit = 14400); 
@save "src/testData/RTS_24_20/gurobiResult.jld2" gurobiResult
@save "src/testData/RTS_24_20/enhancedRTS.jld2" sddipResult_Enhanced
@load "src/testData/RTS_24_20/enhancedRTS.jld2" sddipResult_Enhanced

# using Latexify
# latexify(df; env=:table, latex=false)

## -------------------------------------- Scenario = 50


indexSets = load("src/testData/RTS_24_50/indexSets.jld2")["indexSets"]
paramOPF = load("src/testData/RTS_24_50/paramOPF.jld2")["paramOPF"]
paramDemand = load("src/testData/RTS_24_50/paramDemand.jld2")["paramDemand"]
Ω_rv = load("src/testData/RTS_24_50/Ω_rv.jld2")["Ω_rv"]
prob = load("src/testData/RTS_24_50/prob.jld2")["prob"]

sddipResult_Enhanced = SDDiP_algorithm(; ϵ = 1e-4, max_iter = 100) 
@time gurobiResult = gurobiOptimize!(indexSets, 
                                    paramDemand, 
                                    paramOPF, 
                                    Ω_rv,
                                    prob; timelimit = 14400); 
@save "src/testData/RTS_24_50/gurobiResult.jld2" gurobiResult
@save "src/testData/RTS_24_50/enhancedRTS.jld2" sddipResult_Enhanced
@load "src/testData/RTS_24_50/enhancedRTS.jld2" sddipResult_Enhanced


## -------------------------------------- Scenario = 100

indexSets = load("src/testData/RTS_24_100/indexSets.jld2")["indexSets"]
paramOPF = load("src/testData/RTS_24_100/paramOPF.jld2")["paramOPF"]
paramDemand = load("src/testData/RTS_24_100/paramDemand.jld2")["paramDemand"]
Ω_rv = load("src/testData/RTS_24_100/Ω_rv.jld2")["Ω_rv"]
prob = load("src/testData/RTS_24_100/prob.jld2")["prob"]

sddipResult_Enhanced = SDDiP_algorithm(; ϵ = 1e-4, max_iter = 100) 
@time gurobiResult = gurobiOptimize!(indexSets, 
                                    paramDemand, 
                                    paramOPF, 
                                    Ω_rv,
                                    prob; timelimit = 14400); 
@save "src/testData/RTS_24_100/gurobiResult.jld2" gurobiResult
@save "src/testData/RTS_24_100/enhancedRTS.jld2" sddipResult_Enhanced
@load "src/testData/RTS_24_100/enhancedRTS.jld2" sddipResult_Enhanced


## -------------------------------------- Scenario = 200

indexSets = load("src/testData/RTS_24_200/indexSets.jld2")["indexSets"]
paramOPF = load("src/testData/RTS_24_200/paramOPF.jld2")["paramOPF"]
paramDemand = load("src/testData/RTS_24_200/paramDemand.jld2")["paramDemand"]
Ω_rv = load("src/testData/RTS_24_200/Ω_rv.jld2")["Ω_rv"]
prob = load("src/testData/RTS_24_200/prob.jld2")["prob"]

sddipResult_Enhanced = SDDiP_algorithm(; ϵ = 1e-4, max_iter = 100) 
@time gurobiResult = gurobiOptimize!(indexSets, 
                                    paramDemand, 
                                    paramOPF, 
                                    Ω_rv,
                                    prob; timelimit = 14400);  
@save "src/testData/RTS_24_200/gurobiResult.jld2" gurobiResult
@save "src/testData/RTS_24_200/enhancedRTS.jld2" sddipResult_Enhanced
@load "src/testData/RTS_24_200/enhancedRTS.jld2" sddipResult_Enhanced



## -------------------------------------- Scenario = 500

indexSets = load("src/testData/RTS_24_500/indexSets.jld2")["indexSets"]
paramOPF = load("src/testData/RTS_24_500/paramOPF.jld2")["paramOPF"]
paramDemand = load("src/testData/RTS_24_500/paramDemand.jld2")["paramDemand"]
Ω_rv = load("src/testData/RTS_24_500/Ω_rv.jld2")["Ω_rv"]
prob = load("src/testData/RTS_24_500/prob.jld2")["prob"]

sddipResult_Enhanced = SDDiP_algorithm(; ϵ = 1e-4, max_iter = 100) 
@time gurobiResult = gurobiOptimize!(indexSets, 
                                    paramDemand, 
                                    paramOPF, 
                                    Ω_rv,
                                    prob; timelimit = 14400);  
@save "src/testData/RTS_24_500/gurobiResult.jld2" gurobiResult
@save "src/testData/RTS_24_500/enhancedRTS.jld2" sddipResult_Enhanced
@load "src/testData/RTS_24_500/enhancedRTS.jld2" sddipResult_Enhanced





## ------------------------------------- Comparsion --------------------------------- ##
