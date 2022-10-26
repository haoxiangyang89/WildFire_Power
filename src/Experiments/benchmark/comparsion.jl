include("src/Experiments/benchmark//optimalShutoff.jl")
include("src/Experiments/benchmark/benchmark1_trivial.jl")
include("src/Experiments/benchmark//benchmark2_oracle.jl")

using Plots, Unitful, UnitfulRecipes

indexSets = load("testData_RTS/indexSets.jld2")["indexSets"]
paramOPF = load("testData_RTS/paramOPF.jld2")["paramOPF"]
paramDemand = load("testData_RTS/paramDemand.jld2")["paramDemand"]
Ω_rvList = load("testData_RTS/Ω_rvList.jld2")["Ω_rvList"]
gurobiResultList = load("src/Experiments/ScenariosSizeTest/gurobiResultList.jld2")["gurobiResultList"]
# (13, 200)

zstar = gurobiResultList[200,13].first_state_variable

costShutOffDict = Dict()
costTrivialDict = Dict()
costOracleDict = Dict()

for Ω in [20, 50, 100, 200]
    for i in 1:20
        @info "$Ω, $i"
        Ω_rv = Ω_rvList[(Ω, i)]; indexSets.Ω = [1:Ω...];
        prob = Dict{Int64, Float64}();
        for ω in indexSets.Ω 
            prob[ω] = 1/Ω;
        end
        
        costShutOff = optimalShutOff!(zstar; Ω_rv = Ω_rv);
        costTrivial = benchmarkTrivial!(; Ω_rv = Ω_rv);
        costOracle = benchmarkOracle!(; Ω_rv = Ω_rv);
        costShutOffDict[Ω, i] = costShutOff; costTrivialDict[Ω, i] = costTrivial; costOracleDict[Ω, i] = costOracle; 
        save("src/Experiments/benchmark/costShutOffDict.jld2", "costShutOffDict", costShutOffDict)
        save("src/Experiments/benchmark/costTrivialDict.jld2", "costTrivialDict", costTrivialDict)
        save("src/Experiments/benchmark/costOracleDict.jld2", "costOracleDict", costOracleDict)

        mean(values(costTrivial)) - mean(values(costShutOff))
        mean(values(costTrivial))
        (mean(values(costShutOff)) - mean(values(costOracle)))/mean(values(costShutOff)) * 100
        (mean(values(costTrivial)) - mean(values(costOracle)))/mean(values(costTrivial)) * 100
        (mean(values(costTrivial)) - mean(values(costShutOff)))/mean(values(costTrivial)) * 100

        col_names = [:scenario, :Oracle, :ShutOff, :gapShutOff, :Trivial, :gapTrivial]; # needs to be a vector Symbols
        col_types = [Float64, Float64, Float64, Float64, Float64, Float64];
        named_tuple = (; zip(col_names, type[] for type in col_types )...);
        comparsionResult = DataFrame(named_tuple); # 0×6 DataFrame

        for ω in 1:length(indexSets.Ω)
            push!(comparsionResult, [ω, costOracle[ω], costShutOff[ω], round(costShutOff[ω] - costOracle[ω], digits = 3), costTrivial[ω], round(costTrivial[ω] - costOracle[ω], digits = 3)]);
        end
        comparsionResult[:,:diff] = comparsionResult[:,:gapTrivial] - comparsionResult[:,:gapShutOff]
        comparsionResult[:,:gap] = round.(comparsionResult[:,:diff] ./ comparsionResult[:,:Trivial] * 100, digits = 5)


        ## ============================================ plots part ================================================ ##
        x = comparsionResult[:,1]; oracle = comparsionResult[:,2]; shutoff = comparsionResult[:,3]; trivial = comparsionResult[:,5]; diff = comparsionResult[:,7]; # These are the plotting data
        y = hcat(oracle, shutoff, trivial, diff)

        # Plots.plot(x, diff,  
        #                 label = ["diff"], 
        #                 ls=[ :solid], 
        #                 shape=[:x], 
        #                 lw =2, bg="white", 
        #                 xlab = "Scenario", 
        #                 ylab = "diff"
        #                 )

        Plots.plot(x, y[:,1:3],  
                        label = ["Wait & See" "Shut Off" "No ShutOff"], 
                        ls=[:dot :solid :dot ], 
                        shape=[:star5 :circle :star8], 
                        lw =2, bg="white", # ylims=(minimum(y[:,1:3]), maximum(y[:,1:3])),
                        xlab = "Index of scenarios", 
                        ylab = "Cost"
                        )
        hline!(mean(y[:,1:3], dims=1), label = ["meanWS" "meanSO" "meanNoSO"], line=(2, [:solid :solid :solid], 0.5, [:blue :orange :green]))
    end 
end

Ω = 200
costShutOff_sum = Dict(i => 0 for i in 1:20)
costTrivial_sum = Dict(i => 0 for i in 1:20)
costOracle_sum = Dict(i => 0 for i in 1:20)
for i in 1:20 
    costShutOff_sum =  merge(+, costShutOff_sum, costShutOffDict[Ω, i])
    costTrivial_sum =  merge(+, costTrivial_sum, costTrivialDict[Ω, i])
    costOracle_sum =  merge(+, costOracle_sum, costOracleDict[Ω, i])
end 
(mean(values(costShutOff_sum)) - mean(values(costOracle_sum)))/mean(values(costShutOff_sum)) * 100
(mean(values(costTrivial_sum)) - mean(values(costOracle_sum)))/mean(values(costTrivial_sum)) * 100
(mean(values(costTrivial_sum)) - mean(values(costShutOff_sum)))/mean(values(costTrivial_sum)) * 100




costShutOffDict = load("src/Experiments/benchmark/costShutOffDict.jld2")["costShutOffDict"];
costTrivialDict = load("src/Experiments/benchmark/costTrivialDict.jld2")["costTrivialDict"];
costOracleDict = load("src/Experiments/benchmark/costOracleDict.jld2")["costOracleDict"];

Ω = 20; i = 10;
costShutOff = costShutOffDict[Ω, i];
costTrivial = costTrivialDict[Ω, i];
costOracle = costOracleDict[Ω, i];

mean(values(costShutOff))
mean(values(costTrivial))
mean(values(costOracle))
mean(values(costTrivial)) - mean(values(costShutOff))
(mean(values(costShutOff)) - mean(values(costOracle)))/mean(values(costShutOff)) * 100
(mean(values(costTrivial)) - mean(values(costOracle)))/mean(values(costTrivial)) * 100
(mean(values(costTrivial)) - mean(values(costShutOff)))/mean(values(costTrivial)) * 100

col_names = [:scenario, :Oracle, :ShutOff, :gapShutOff, :Trivial, :gapTrivial]; # needs to be a vector Symbols
col_types = [Float64, Float64, Float64, Float64, Float64, Float64];
named_tuple = (; zip(col_names, type[] for type in col_types )...);
comparsionResult = DataFrame(named_tuple); # 0×6 DataFrame

for ω in 1:Ω
    push!(comparsionResult, [ω, costOracle[ω], costShutOff[ω], round(costShutOff[ω] - costOracle[ω], digits = 3), costTrivial[ω], round(costTrivial[ω] - costOracle[ω], digits = 3)]);
end
comparsionResult[:,:diff] = comparsionResult[:,:gapTrivial] - comparsionResult[:,:gapShutOff]
comparsionResult[:,:gap] = round.(comparsionResult[:,:diff] ./ comparsionResult[:,:Trivial] * 100, digits = 5)


## ============================================ plots part ================================================ ##
x = comparsionResult[:,1]; oracle = comparsionResult[:,2]; shutoff = comparsionResult[:,3]; trivial = comparsionResult[:,5]; diff = comparsionResult[:,7]; # These are the plotting data
y = hcat(oracle, shutoff, trivial, diff)

# Plots.plot(x, diff,  
#                 label = ["diff"], 
#                 ls=[ :solid], 
#                 shape=[:x], 
#                 lw =2, bg="white", 
#                 xlab = "Scenario", 
#                 ylab = "diff"
#                 )

Plots.plot(x, y[:,1:3],  
            label = ["Wait & See" "Shut Off" "No ShutOff"], 
            ls=[:dot :solid :dot ], 
            shape=[:star5 :circle :star8], 
            lw =2, bg="white", # ylims=(0, 1e4),
            xlab = "Index of scenarios", 
            ylab = "Cost"
            )
hline!(mean(y[:,1:3], dims=1), label = ["meanWS" "meanSO" "meanNoSO"], line=(2, [:solid :solid :solid], 0.5, [:blue :orange :green]))
# https://docs.juliahub.com/UnitfulRecipes/KPSlU/1.0.0/examples/2_Plots/#




# p1 = Plots.plot(x, y,  
#                 label = ["Wait & See" "Shut Off" "Trivial" "diff"], 
#                 ls=[:dash :solid :dash :solid], 
#                 shape=[:star5 :circle :star8 :x], 
#                 lw =2, bg="white", ylims=(6e7, 1.12e8),
#                 xlab = "Scenario", 
#                 ylab = "Cost"
#                 )
# p2 = Plots.plot(x, diff,  
#                 shape=[:x], 
#                 lw =2, bg="white", 
#                 xlab = "Scenario", 
#                 ylab = "Progress",
#                 linecolor = [:black]
#                 )


# Plots.plot(p1, p2, layout = (1, 2), legend = false)


