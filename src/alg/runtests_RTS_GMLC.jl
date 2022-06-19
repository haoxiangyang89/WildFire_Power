
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

# save("indexSets.jld2", "indexSets", indexSets)
# save("paramOPF.jld2", "paramOPF", paramOPF)
# save("paramDemand.jld2", "paramDemand", paramDemand)
# save("Ω_rv.jld2", "Ω_rv", Ω_rv)
# save("prob.jld2", "prob", prob)



# @passobj 1 workers() indexSets
# @passobj 1 workers() paramOPF
# @passobj 1 workers() paramDemand
# @passobj 1 workers() Ω_rv
# @passobj 1 workers() prob




# #############################################################################################################
# @everywhere begin
#     max_iter = 200; ϵ = 1e-3; Enhanced_Cut = true;

#     λ_value = .1; Output = 0; Output_Gap = false; Adj = false; Enhanced_Cut = true; threshold = 1e2; 
#     levelSetMethodParam = LevelSetMethodParam(0.95, λ_value, threshold, 1e14, 3e3, Output, Output_Gap, Adj)
# end



# resultDict = SDDiP_algorithm(Ω_rv, prob, 
#                     indexSets, 
#                     paramDemand, 
#                     paramOPF; 
#                     levelSetMethodParam = levelSetMethodParam,
#                     ϵ = 0.001, M = 1, max_iter = 200, 
#                     Enhanced_Cut = true)





















