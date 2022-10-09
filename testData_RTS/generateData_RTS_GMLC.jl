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
Ω = 500; ## Int
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

# @load "src/testData/emptyScenario.jld2" emptyScenario
Ω_rv[9] = emptyScenario

prob = Dict{Int64, Float64}();
for ω in indexSets.Ω 
    prob[ω] = 1/Ω;
end


save("testData_RTS/indexSets.jld2", "indexSets", indexSets)
save("testData_RTS/paramOPF.jld2", "paramOPF", paramOPF)
save("testData_RTS/paramDemand.jld2", "paramDemand", paramDemand)
# save("testData_RTS/Ω_rv.jld2", "Ω_rv", Ω_rv)
save("testData_RTS/prob.jld2", "prob", prob)





T = 24;
Ω = 5000; ## Int
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
save("testData_RTS/Ω_rv5000.jld2", "Ω_rv", Ω_rv)


T = 24;
Ω = 10000; ## Int
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
save("testData_RTS/Ω_rv10000.jld2", "Ω_rv", Ω_rv)

# Ω_rvList = Dict()
# for num in [20, 50, 100, 200, 500] 
#     for i in 1:20 
#         T = 24;
#         Ω = num; ## Int
#         (indexSets, paramOPF, paramDemand, multiLines) = prepareIndexSets(network_data, T, Ω);
#         (x_grid_num, y_grid_num, line_location_id, 
#                                     line_id_location, 
#                                     line_id_bus,
#                                     bus_id_location, 
#                                     bus_location_id) = prepareSimulation(businfo, branchInfo, WFPI_Info; n = 2, grid_length = 3000);

        

#         Ω_rvList[num, i] = prepareScenarios( ;period_span = 1, 
#                                             T = T, 
#                                             Ω = Ω, 
#                                             indexSets = indexSets, line_id_bus = line_id_bus);
#     end
#     save("testData_RTS/Ω_rvList.jld2", "Ω_rvList", Ω_rvList)
# end
