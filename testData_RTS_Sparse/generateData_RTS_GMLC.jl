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
Ω = 101; ## Int
(indexSets, paramOPF, paramDemand, multiLines) = prepareIndexSets(network_data, T, Ω);



(x_grid_num, y_grid_num, line_location_id, 
                            line_id_location, 
                            line_id_bus,
                            bus_id_location, 
                            bus_location_id) = prepareSimulation(businfo, branchInfo, WFPI_Info; n = 2, grid_length = 3000);

 

Ω_rv = prepareScenarios( ;period_span = 1, 
                                    T = T, 
                                    Ω = Ω, 
                                    indexSets = indexSets, line_id_bus = line_id_bus)

@load "emptyScenario.jld2" emptyScenario
Ω_rv[Ω] = emptyScenario

prob = Dict{Int64, Float64}();
for ω in 1:(Ω - 1) 
    prob[ω] = round(0.1/(Ω - 1), digits = 4);
end
prob[Ω] = .9

save("testData_RTS_Sparse/indexSets.jld2", "indexSets", indexSets)
save("testData_RTS_Sparse/paramOPF.jld2", "paramOPF", paramOPF)
save("testData_RTS_Sparse/paramDemand.jld2", "paramDemand", paramDemand)
save("testData_RTS_Sparse/Ω_rv.jld2", "Ω_rv", Ω_rv)
save("testData_RTS_Sparse/prob.jld2", "prob", prob)
# save("testData_RTS_Sparse/wholeSpace.jld2", "wholeSpace", wholeSpace)
# @load "testData_RTS_Sparse/wholeSpace.jld2" wholeSpace



Ω_rvList = Dict()
for N in [50, 100, 200, 500] 
    for i in 1:20 
        SAAscenario = rand!(collect(1:N),collect(1:1000));
        saa = findall(x -> x <= 100, SAAscenario);
        Ω_rv = Dict{Int64, RandomVariables}();
        for i in 1:N 
            Ω_rv[i] = wholeSpace[101]
        end 
        for i in saa 
            Ω_rv[i] = wholeSpace[SAAscenario[i]]
        end

        Ω_rvList[N, i] = Ω_rv
    end
    save("testData_RTS_Sparse/Ω_rvList.jld2", "Ω_rvList", Ω_rvList)
end




ignitionList = Dict()
for N in [50] 
    for i in 1:20 
        SAAscenario = rand!(collect(1:N),collect(1:100));
        Ω_rv = Dict{Int64, RandomVariables}();
        for i in 1:N 
            Ω_rv[i] = wholeSpace[SAAscenario[i]]
        end 

        Ω_rv[N] = wholeSpace[101]

        ignitionList[N, i] = Ω_rv
    end
    save("testData_RTS_Sparse/ignitionList.jld2", "ignitionList", ignitionList)
end




# T = 24;
# Ω = 5000; ## Int
# (indexSets, paramOPF, paramDemand, multiLines) = prepareIndexSets(network_data, T, Ω);
# (x_grid_num, y_grid_num, line_location_id, 
#                             line_id_location, 
#                             line_id_bus,
#                             bus_id_location, 
#                             bus_location_id) = prepareSimulation(businfo, branchInfo, WFPI_Info; n = 2, grid_length = 3000);

 

# Ω_rv = prepareScenarios( ;period_span = 1, 
#                                     T = T, 
#                                     Ω = Ω, 
#                                     indexSets = indexSets, line_id_bus = line_id_bus);
# save("testData_RTS_Sparse/Ω_rv5000.jld2", "Ω_rv", Ω_rv)


# T = 24;
# Ω = 10000; ## Int
# (indexSets, paramOPF, paramDemand, multiLines) = prepareIndexSets(network_data, T, Ω);
# (x_grid_num, y_grid_num, line_location_id, 
#                             line_id_location, 
#                             line_id_bus,
#                             bus_id_location, 
#                             bus_location_id) = prepareSimulation(businfo, branchInfo, WFPI_Info; n = 2, grid_length = 3000);

 

# Ω_rv = prepareScenarios( ;period_span = 1, 
#                                     T = T, 
#                                     Ω = Ω, 
#                                     indexSets = indexSets, line_id_bus = line_id_bus);
# save("testData_RTS_Sparse/Ω_rv10000.jld2", "Ω_rv", Ω_rv)


# num = 100; i = 9, 15, 14

# for num in [20, 50, 100, 200, 500] 
#     for i in [9, 15, 14]
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
#     save("testData_RTS_Sparse/Ω_rvList.jld2", "Ω_rvList", Ω_rvList)
# end
