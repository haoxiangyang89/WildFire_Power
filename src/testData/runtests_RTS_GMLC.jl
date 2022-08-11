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
(indexSets, paramOPF, paramDemand) = prepareIndexSets(network_data, T, Ω);




(x_grid_num, y_grid_num, line_location_id, 
                            line_id_location, 
                            line_id_bus,
                            bus_id_location, 
                            bus_location_id) = prepareSimulation(businfo, branchInfo, WFPI_Info; n = 10, grid_length = 3000);

 



Ω_rv = prepareScenarios( ;period_span = 4, 
                                    T = T, 
                                    Ω = Ω, 
                                    indexSets = indexSets, line_id_bus = line_id_bus);



prob = Dict{Int64, Float64}();
for ω in indexSets.Ω 
    prob[ω] = 1/Ω;
end

save("src/testData/RTS_24_500/indexSets.jld2", "indexSets", indexSets)
save("src/testData/RTS_24_500/paramOPF.jld2", "paramOPF", paramOPF)
save("src/testData/RTS_24_500/paramDemand.jld2", "paramDemand", paramDemand)
save("src/testData/RTS_24_500/Ω_rv.jld2", "Ω_rv", Ω_rv)
save("src/testData/RTS_24_500/prob.jld2", "prob", prob)