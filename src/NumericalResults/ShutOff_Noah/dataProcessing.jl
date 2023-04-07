
# load packages 
import JSON
using JuMP, Gurobi, PowerModels;
using CSV, DataFrames, Printf;
using JLD2, FileIO, JSON;


include("src/alg/def.jl")

## ----------------------------------------------- Noah Solution ------------------------------------------------ ##

solution = load("testData_RTS_Sparse/Solution.jld2")["Solution"].first_state_variable;
network_data = PowerModels.parse_file("data/RTS_GMLC/case_RTS_GMLC.m");

NoahSolutionList = Dict()
heuristicSolutionList = Dict()
for alpha in [10, 20, 30, 40, 50, 60, 70, 80, 90]
    solution[:zb] .= 1; 
    solution[:zg] .= 1; 
    solution[:zl] .= 1; 

    Shutoff = JSON.parsefile("testData_RTS_Sparse/Noah/alpha_$alpha.json")

    for t in 1:24 
        shutoff_period_t = Shutoff["$t"]
        for component in keys(shutoff_period_t) 
            shutOffSet = shutoff_period_t[component]
            if !isempty(shutOffSet)
                if component == "gen"
                    for gen in shutOffSet 
                        g = network_data[component][gen]["index"]
                        solution[:zg][g, t:24] =  [0. for i in t:24]
                    end
                elseif component == "branch"
                    for branch in shutOffSet 
                        l = (network_data[component][branch]["f_bus"], network_data[component][branch]["t_bus"])
                        solution[:zl][l, t:24] =  [0. for i in t:24]
                    end
                elseif component == "bus"
                    for bus in shutOffSet 
                        b = network_data[component][bus]["index"]
                        solution[:zb][b, t:24] =  [0. for i in t:24]
                    end
                end
            end
        end
    end

    heuristicSolutionList[alpha] = deepcopy(solution)
end
save("src/Experiments/ShutOff_Noah/heuristicSolutionList.jld2", "heuristicSolutionList", heuristicSolutionList)




for alpha in [10, 20, 30, 40, 50, 60, 70, 80, 90]
    solution[:zb] .= 1; 
    solution[:zg] .= 1; 
    solution[:zl] .= 1; 

    Shutoff = JSON.parsefile("testData_RTS_Sparse/Noah/alpha_$alpha.json")

    for t in 1:24 
        shutoff_period_t = Shutoff["$t"]
        for component in keys(shutoff_period_t) 
            shutOffSet = shutoff_period_t[component]
            if !isempty(shutOffSet)
                if component == "gen"
                    for gen in shutOffSet 
                        g = network_data[component][gen]["index"]
                        solution[:zg][g, t] =  0.0 # [0. for i in t:24]
                    end
                elseif component == "branch"
                    for branch in shutOffSet 
                        l = (network_data[component][branch]["f_bus"], network_data[component][branch]["t_bus"])
                        solution[:zl][l, t] =  0.0 # [0. for i in t:24]
                    end
                elseif component == "bus"
                    for bus in shutOffSet 
                        b = network_data[component][bus]["index"]
                        solution[:zb][b, t] =  0.0 # [0. for i in t:24]
                    end
                end
            end
        end
    end

    NoahSolutionList[alpha] = deepcopy(solution)
end
save("src/Experiments/ShutOff_Noah/NoahSolutionList.jld2", "NoahSolutionList", NoahSolutionList)

## ----------------------------------------------- Wildfire Risk Value ------------------------------------------------ ##

riskValue = JSON.parsefile("src/Experiments/ShutOff_Noah/Noah/wildfire_risk.json")
riskValueDict = load("testData_RTS_Sparse/Solution.jld2")["Solution"].first_state_variable;

for t in 1:24 
    riskValue_t = riskValue["$t"]
    for b in keys(riskValue_t["bus"])
        riskValueDict[:zb][parse(Int64, b), t] = riskValue_t["bus"][b]
    end

    for g in keys(riskValue_t["gen"])
        riskValueDict[:zg][parse(Int64, g), t] = riskValue_t["gen"][g]
    end

    for l in keys(riskValue_t["line"])
        riskValueDict[:zl][(parse(Int, l[2:4]), parse(Int, l[7:9])), t] = riskValue_t["line"][l]
    end
end

save("src/Experiments/ShutOff_Noah/riskValueDict.jld2", "riskValueDict", riskValueDict)


