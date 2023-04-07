## generate data, simulation
using JuMP, Gurobi, PowerModels;
using Statistics, StatsBase, Random, Dates, Distributions;
using CSV, DataFrames, Printf;
using JLD2, FileIO;


include("src/alg/CellularAutomaton.jl");
include("src/alg/def.jl");
include("src/alg/readin.jl");

network_data = PowerModels.parse_file("data/RTS_GMLC/case_RTS_GMLC.m");
businfo = CSV.read("data/RTS_GMLC/bus.csv", DataFrame);
branchInfo = CSV.read("data/RTS_GMLC/branch.csv", DataFrame);
WFPI_file = CSV.read("data/RTS_GMLC/RTS_GMLC_WFPI_with_Mean.csv", DataFrame);
WFPI_Info = WFPI_file[:, [:From_Bus, :To_Bus, :Mean]];



T = 24;
Ω = 20; ## Int
(indexSets, paramOPF, paramDemand, multiLines) = prepareIndexSets(network_data, T, Ω);



(x_grid_num, y_grid_num, line_location_id, 
                            line_id_location, 
                            line_id_bus,
                            bus_id_location, 
                            bus_location_id) = prepareSimulation(businfo, branchInfo, WFPI_Info; n = 2, grid_length = 3000);

 

Ω_rv = prepareScenarios( ;period_span = 1, 
                                    T = T, 
                                    Ω = Ω, 
                                    indexSets = indexSets, line_id_bus = line_id_bus);


# @load "emptyScenario.jld2" emptyScenario
# Ω_rv[9] = emptyScenario

prob = Dict{Int64, Float64}();
for ω in indexSets.Ω 
    prob[ω] = 1/Ω;
end


# save("data/testData_RTS/sampleSize_20/indexSets.jld2", "indexSets", indexSets)
# save("data/testData_RTS/sampleSize_20/paramOPF.jld2", "paramOPF", paramOPF)
# save("data/testData_RTS/sampleSize_20/paramDemand.jld2", "paramDemand", paramDemand)
save("data/testData_RTS/sampleSize_20/Ω_rv.jld2", "Ω_rv", Ω_rv)
save("data/testData_RTS/sampleSize_20/prob.jld2", "prob", prob)



# T = 24;
# Ω = 5000; ## Int
# (x_grid_num, y_grid_num, line_location_id, 
#                             line_id_location, 
#                             line_id_bus,
#                             bus_id_location, 
#                             bus_location_id) = prepareSimulation(businfo, branchInfo, WFPI_Info; n = 2, grid_length = 3000);

# Ω_rv = prepareScenarios( ;period_span = 1, 
#                                     T = T, 
#                                     Ω = Ω, 
#                                     indexSets = indexSets, line_id_bus = line_id_bus);
# save("data/testData_RTS/Ω_rv5000.jld2", "Ω_rv", Ω_rv)


# Ω_rvList = Dict()
# for num in [20, 50, 100, 200, 500] 
#     for i in 1:20 
#         T = 24;
#         Ω = num; ## Int

#         Ω_rvList[num, i] = prepareScenarios( ;period_span = 1, 
#                                             T = T, 
#                                             Ω = Ω, 
#                                             indexSets = indexSets, line_id_bus = line_id_bus);
#     end
    # save("data/testData_RTS/Ω_rvList.jld2", "Ω_rvList", Ω_rvList)
# end


# probList = Dict{Any, Dict{Int64, Float64}}()
# prob = Dict{Int64, Float64}()
# Ω_rv_New_List = Dict()
# for num in [20, 50, 100, 200, 500] 
#     for i in 1:20 
#         Ω_rv = Ω_rvList[num, i]
#         Ω_rv_new = Dict{Int64, RandomVariables}()
#         prob = Dict{Int64, Float64}()
#         num_empty = 0 
#         num_inex = 0
#         for ω in keys(Ω_rv)
#             if Ω_rv[ω].τ == 24 
#                 num_empty = num_empty + 1
#             else
#                 num_inex = num_inex + 1
#                 Ω_rv_new[num_inex] = Ω_rv[ω]
#                 prob[num_inex] = 1/num
#             end
#         end
#         Ω_rv_new[num_inex+1] = emptyScenario
#         prob[num_inex+1] = round(1 - num_inex/num, digits = 6)
#         probList[(num, i)] = prob
#         Ω_rv_New_List[num, i] = Ω_rv_new
#     end
# end
# save("data/testData_RTS/prob5000_reduced.jld2", "prob", prob)
# save("data/testData_RTS/Ω_rv5000_reduced.jld2", "Ω_rv", Ω_rv_new)

