include("src/benchmark/optimalShutoff.jl")
include("src/benchmark/benchmark1_trivial.jl")
include("src/benchmark/benchmark2_oracle.jl")

indexSets = load("testData_RTS/indexSets.jld2")["indexSets"]
paramOPF = load("testData_RTS/paramOPF.jld2")["paramOPF"]
paramDemand = load("testData_RTS/paramDemand.jld2")["paramDemand"]
Ω_rvList = load("testData_RTS/Ω_rvList.jld2")["Ω_rvList"]


costShutOffDict = Dict()
costTrivialDict = Dict()
costOracleDict = Dict()

for Ω in [20, 50, 100, 200]
    for i in 1:1
        Ω_rv = Ω_rvList[(Ω, i)]; indexSets.Ω = [1:Ω...];
        costShutOff = optimalShutOff!();
        costTrivial = benchmarkTrivial!();
        costOracle = benchmarkOracle!();
        costShutOffDict[Ω] = costShutOff; costTrivialDict[Ω] = costTrivial; costOracleDict[Ω] = costOracle; 
        save("src/Experiments/benchmark/costShutOffDict.jld2", "costShutOffDict", costShutOffDict)
        save("src/Experiments/benchmark/costTrivialDict.jld2", "costTrivialDict", costTrivialDict)
        save("src/Experiments/benchmark/costOracleDict.jld2", "costOracleDict", costOracleDict)

        mean(values(costTrivial)) - mean(values(costShutOff))
        mean(values(costTrivial))
        (mean(values(costShutOff)) - mean(values(costOracle)))/mean(values(costOracle)) * 100
        (mean(values(costTrivial)) - mean(values(costOracle)))/mean(values(costOracle)) * 100
        (mean(values(costTrivial)) - mean(values(costShutOff)))/mean(values(costOracle)) * 100

        col_names = [:scenario, :Oracle, :ShutOff, :gapShutOff, :Trivial, :gapTrivial]; # needs to be a vector Symbols
        col_types = [Float64, Float64, Float64, Float64, Float64, Float64];
        named_tuple = (; zip(col_names, type[] for type in col_types )...);
        comparsionResult = DataFrame(named_tuple); # 0×6 DataFrame

        for ω in 1:length(indexSets.Ω)
            push!(comparsionResult, [ω, costOracle[ω], costShutOff[ω], round(costShutOff[ω] - costOracle[ω], digits = 3), costTrivial[ω], round(costTrivial[ω] - costOracle[ω], digits = 3)]);
        end
        comparsionResult[:,:diff] = comparsionResult[:,:gapTrivial] - comparsionResult[:,:gapShutOff]
        comparsionResult[:,:gap] = round.(comparsionResult[:,:diff] ./ comparsionResult[:,:Oracle] * 100, digits = 5)


        ## ============================================ plots part ================================================ ##
        using Plots, Unitful, UnitfulRecipes
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
                        lw =2, bg="white", ylims=(miniumu(y[:,1:3]), maximum(y[:,1:3])),
                        xlab = "Index of scenarios", 
                        ylab = "Cost"
                        )
        hline!(mean(y[:,1:3], dims=1), label = ["meanWS" "meanSO" "meanNoSO"], line=(2, [:solid :solid :solid], 0.5, [:blue :orange :green]))
    end 
end

# costShutOff = optimalShutOff!();
# costTrivial = benchmarkTrivial!();
# costOracle = benchmarkOracle!();

# mean(values(costTrivial)) - mean(values(costShutOff))
# mean(values(costTrivial))
# (mean(values(costShutOff)) - mean(values(costOracle)))/mean(values(costOracle)) * 100
# (mean(values(costTrivial)) - mean(values(costOracle)))/mean(values(costOracle)) * 100
# (mean(values(costTrivial)) - mean(values(costShutOff)))/mean(values(costOracle)) * 100

# col_names = [:scenario, :Oracle, :ShutOff, :gapShutOff, :Trivial, :gapTrivial]; # needs to be a vector Symbols
# col_types = [Float64, Float64, Float64, Float64, Float64, Float64];
# named_tuple = (; zip(col_names, type[] for type in col_types )...);
# comparsionResult = DataFrame(named_tuple); # 0×6 DataFrame

# for ω in 1:length(indexSets.Ω)
#     push!(comparsionResult, [ω, costOracle[ω], costShutOff[ω], round(costShutOff[ω] - costOracle[ω], digits = 3), costTrivial[ω], round(costTrivial[ω] - costOracle[ω], digits = 3)]);
# end
# comparsionResult[:,:diff] = comparsionResult[:,:gapTrivial] - comparsionResult[:,:gapShutOff]
# comparsionResult[:,:gap] = round.(comparsionResult[:,:diff] ./ comparsionResult[:,:Oracle] * 100, digits = 5)


# ## ============================================ plots part ================================================ ##
# using Plots, Unitful, UnitfulRecipes
# x = comparsionResult[:,1]; oracle = comparsionResult[:,2]; shutoff = comparsionResult[:,3]; trivial = comparsionResult[:,5]; diff = comparsionResult[:,7]; # These are the plotting data
# y = hcat(oracle, shutoff, trivial, diff)

# Plots.plot(x, diff,  
#                 label = ["diff"], 
#                 ls=[ :solid], 
#                 shape=[:x], 
#                 lw =2, bg="white", 
#                 xlab = "Scenario", 
#                 ylab = "diff"
#                 )

# Plots.plot(x, y[:,1:3],  
#                 label = ["Wait & See" "Shut Off" "No ShutOff"], 
#                 ls=[:dot :solid :dot ], 
#                 shape=[:star5 :circle :star8], 
#                 lw =2, bg="white", ylims=(2e5, 3.1e5),
#                 xlab = "Index of scenarios", 
#                 ylab = "Cost"
#                 )
# hline!(mean(y[:,1:3], dims=1), label = ["meanWS" "meanSO" "meanNoSO"], line=(2, [:solid :solid :solid], 0.5, [:blue :orange :green]))
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


