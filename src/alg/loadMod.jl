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


const GRB_ENV = Gurobi.Env()


include("src/alg/def.jl")
include("src/alg/backwardPass.jl")
include("src/alg/forwardPass.jl")
include("src/alg/extFormGurobi.jl")

include("src/alg/CellularAutomaton.jl")
include("src/alg/readin.jl")

include("src/alg/runtests_RTS_GMLC.jl")
include("src/alg/sddip.jl")
# include("src/alg/runtests_case30.jl") 