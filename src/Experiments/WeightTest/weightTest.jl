## This test is to evaluate the influence caused by different weight of no wildfire
## RTS case, 24 periods, 20 scenarios
## 19 wildfire ignition + 1 no wildfire
## Ω_rv[9] is the one without ignition

indexSets = load("src/testData/RTS_24_20/indexSets.jld2")["indexSets"]
paramOPF = load("src/testData/RTS_24_20/paramOPF.jld2")["paramOPF"]
paramDemand = load("src/testData/RTS_24_20/paramDemand.jld2")["paramDemand"]
@load "src/Experiments/WeightTest/Ω_rv_no_wildfire.jld2" Ω_rv


Ω = 20;
probList = Dict{Float64, Any}()
for p in [0.0, .3, .5, .7, .9]
    prob = Dict{Int64, Float64}();
    for ω in 1:Ω
        prob[ω] = (1 - p)/Ω;
    end
    prob[9] = prob[9] + p  
    probList[p] = prob
end


gurobiResultList = Dict{Float64, Any}();
for p in [0.0, .3, .5, .7, .9]
    prob = probList[p] 
    @time gurobiResult = gurobiOptimize!(indexSets, 
                                        paramDemand, 
                                        paramOPF, 
                                        Ω_rv,
                                        prob);  
    gurobiResultList[p] = gurobiResult;
end

save("src/Experiments/WeightTest/gurobiResultList.jld2", "gurobiResultList", gurobiResultList)
