
using JLD2, FileIO

network_data = PowerModels.parse_file("data/RTS_GMLC/case_RTS_GMLC.m");
businfo = CSV.read("data/RTS_GMLC/bus.csv", DataFrame);
branchInfo = CSV.read("data/RTS_GMLC/branch.csv", DataFrame);
WFPI_file = CSV.read("data/RTS_GMLC/RTS_GMLC_WFPI_with_Mean.csv", DataFrame);
WFPI_Info = WFPI_file[:, [:From_Bus, :To_Bus, :Mean]];



(x_grid_num, y_grid_num, line_location_id, 
line_id_location, 
line_id_bus,
bus_id_location, 
bus_location_id) = prepareSimulation(businfo, branchInfo, WFPI_Info; n = 10, grid_length = 3000);

iter = 1; Ω_List = [20, 50, 100]; gurobiResultList = Dict{Int64, Any}();
while iter <= 30
    T = 24; Ω = 20;
    (indexSets, paramOPF, paramDemand) = prepareIndexSets(network_data, T, Ω);
    Ω_rv = prepareScenarios( ;period_span = 4, 
                                        T = T, 
                                        Ω = Ω, 
                                        indexSets = indexSets, line_id_bus = line_id_bus);

    prob = Dict{Int64, Float64}();
    for ω in indexSets.Ω 
        prob[ω] = 1/Ω;
    end

    @time gurobiResult = gurobiOptimize!(indexSets, 
                                    paramDemand, 
                                    paramOPF, 
                                    Ω_rv,
                                    prob; timelimit = 7200); 
    gurobiResultList[(iter, Ω)] = gurobiResult;
    iter = iter + 1;
end


iter = 1; 
while iter <= 30
    T = 24; Ω = 50;
    (indexSets, paramOPF, paramDemand) = prepareIndexSets(network_data, T, Ω);
    Ω_rv = prepareScenarios( ;period_span = 4, 
                                        T = T, 
                                        Ω = Ω, 
                                        indexSets = indexSets, line_id_bus = line_id_bus);



    prob = Dict{Int64, Float64}();
    for ω in indexSets.Ω 
        prob[ω] = 1/Ω;
    end

    @time gurobiResult = gurobiOptimize!(indexSets, 
                                    paramDemand, 
                                    paramOPF, 
                                    Ω_rv,
                                    prob; timelimit = 7200); 
    gurobiResultList[(iter, Ω)] = gurobiResult;
    iter = iter + 1;
end



iter = 1; 
while iter <= 30
    T = 24; Ω = 100;
    (indexSets, paramOPF, paramDemand) = prepareIndexSets(network_data, T, Ω);
    Ω_rv = prepareScenarios( ;period_span = 4, 
                                        T = T, 
                                        Ω = Ω, 
                                        indexSets = indexSets, line_id_bus = line_id_bus);



    prob = Dict{Int64, Float64}();
    for ω in indexSets.Ω 
        prob[ω] = 1/Ω;
    end

    @time gurobiResult = gurobiOptimize!(indexSets, 
                                    paramDemand, 
                                    paramOPF, 
                                    Ω_rv,
                                    prob; timelimit = 14400); 
    gurobiResultList[(iter, Ω)] = gurobiResult;
    iter = iter + 1;
end

save("src/Experiments/ScenariosSizeTest/totalCost.jld2", "totalCost", totalCost)
