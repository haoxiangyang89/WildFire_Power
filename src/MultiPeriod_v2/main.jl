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


include("src/MultiPeriod_v2/data_struct.jl")
include("src/MultiPeriod_v2/backward_pass.jl")
include("src/MultiPeriod_v2/forward_pass.jl")
include("src/MultiPeriod_v2/gurobiTest.jl")

include("src/MultiPeriod_v2/CellularAutomaton.jl")
include("src/MultiPeriod_v2/readin.jl")

include("src/MultiPeriod_v2/runtests_RTS_GMLC.jl")
include("src/MultiPeriod_v2/SDDiP.jl")
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
    
    include("src/MultiPeriod_v2/data_struct.jl")
    include("src/MultiPeriod_v2/backward_pass.jl")
    include("src/MultiPeriod_v2/forward_pass.jl")
    include("src/MultiPeriod_v2/gurobiTest.jl")
end

include("src/MultiPeriod_v2/CellularAutomaton.jl")
include("src/MultiPeriod_v2/readin.jl")

# include("src/MultiPeriod_v2/runtests_case30.jl")  
include("src/MultiPeriod_v2/runtests_RTS_GMLC.jl") 
include("src/MultiPeriod_v2/SDDiP_parallel.jl") 
