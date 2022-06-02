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
using CSV, DataFrames


const GRB_ENV = Gurobi.Env()


include("src/MultiPeriod_v2/def.jl")
include("src/MultiPeriod_v2/backwardPass.jl")
include("src/MultiPeriod_v2/forwardPass.jl")
include("src/MultiPeriod_v2/extFormGurobi.jl")

include("src/MultiPeriod_v2/CellularAutomaton.jl")
include("src/MultiPeriod_v2/readin.jl")

include("src/MultiPeriod_v2/runtests_RTS_GMLC.jl")
include("src/MultiPeriod_v2/sddip.jl")
# include("src/MultiPeriod_v2/runtests_case30.jl") 


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
    using DataFrames

    # using Agents, Distributions
    # using CSV, Geodesy
    # using InteractiveDynamics
    # using CairoMakie
    # using PowerModels


    const GRB_ENV = Gurobi.Env()
    
    include("src/MultiPeriod_v2/def.jl")
    include("src/MultiPeriod_v2/backwardPass.jl")
    include("src/MultiPeriod_v2/forwardPass.jl")
    include("src/MultiPeriod_v2/extFormGurobi.jl")
end

include("src/MultiPeriod_v2/CellularAutomaton.jl")
include("src/MultiPeriod_v2/readin.jl")

# include("src/MultiPeriod_v2/runtests_case30.jl")  
include("src/MultiPeriod_v2/runtests_RTS_GMLC.jl") 
include("src/MultiPeriod_v2/sddipParallel.jl") 
