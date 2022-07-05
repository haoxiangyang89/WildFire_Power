
using JLD2, FileIO

network_data = PowerModels.parse_file("data/RTS_GMLC/case_RTS_GMLC.m");
businfo = CSV.read("data/RTS_GMLC/bus.csv", DataFrame);
branchInfo = CSV.read("data/RTS_GMLC/branch.csv", DataFrame);
WFPI_file = CSV.read("data/RTS_GMLC/RTS_GMLC_WFPI_with_Mean.csv", DataFrame);
WFPI_Info = WFPI_file[:, [:From_Bus, :To_Bus, :Mean]];

T = 24;
Ω = 20; ## Int
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

save("indexSets.jld2", "indexSets", indexSets)
save("paramOPF.jld2", "paramOPF", paramOPF)
save("paramDemand.jld2", "paramDemand", paramDemand)
save("Ω_rv.jld2", "Ω_rv", Ω_rv)
save("prob.jld2", "prob", prob)