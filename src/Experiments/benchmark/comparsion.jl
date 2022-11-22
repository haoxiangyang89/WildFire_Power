include("src/Experiments/benchmark//optimalShutoff.jl")
include("src/Experiments/benchmark/benchmark1_trivial.jl")
include("src/Experiments/benchmark//benchmark2_oracle.jl")

using Plots, Unitful, UnitfulRecipes

indexSets = load("testData_RTS_Sparse/indexSets.jld2")["indexSets"]
paramOPF = load("testData_RTS_Sparse/paramOPF.jld2")["paramOPF"]
paramDemand = load("testData_RTS_Sparse/paramDemand.jld2")["paramDemand"]
Ω_rv = load("testData_RTS_Sparse/Ω_rv.jld2")["Ω_rv"]
Solution = load("testData_RTS_Sparse/Solution.jld2")["Solution"]
zstar = Solution.first_state_variable

prob = Dict{Int64, Float64}(); Ω = length(indexSets.Ω);
for ω in 1:(Ω - 1) 
    prob[ω] = round(0.1/(Ω - 1), digits = 4);
end
prob[Ω] = .9

costShutOff = optimalShutOff!(zstar; Ω_rv = Ω_rv);
costTrivial = benchmarkTrivial!(; Ω_rv = Ω_rv);
costOracle = benchmarkOracle!(; Ω_rv = Ω_rv);

# costOracleDict = load("src/Experiments/benchmark/costOracleDict.jld2")["costOracleDict"]
# costShutOffDict = load("src/Experiments/benchmark/costShutOffDict.jld2")["costShutOffDict"]
# costTrivialDict = load("src/Experiments/benchmark/costTrivialDict.jld2")["costTrivialDict"]

# n = 20;m = 11;
# costOracle = costOracleDict[n,m]
# costShutOff = costShutOffDict[n,m]
# costTrivial = costTrivialDict[n,m]

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

for ω in 1:length(indexSets.Ω)
    push!(comparsionResult, [ω, costOracle[ω], costShutOff[ω], round(costShutOff[ω] - costOracle[ω], digits = 3), costTrivial[ω], round(costTrivial[ω] - costOracle[ω], digits = 3)]);
end
comparsionResult[:,:diff] = comparsionResult[:,:gapTrivial] - comparsionResult[:,:gapShutOff]
comparsionResult[:,:gap] = round.(comparsionResult[:,:diff] ./ comparsionResult[:,:Trivial] * 100, digits = 5)


## ============================================ plots part ================================================ ##
x = comparsionResult[:,1]; oracle = comparsionResult[:,2]; shutoff = comparsionResult[:,3]; trivial = comparsionResult[:,5]; diff = comparsionResult[:,7]; # These are the plotting data
y = hcat(oracle, shutoff, trivial, diff)


p = Plots.plot(x[5:15], y[5:15,1:3],  
                label = ["Wait & See" "Shut Off" "Det"], 
                ls=[:dot :solid :dot ], 
                shape=[:star5 :circle :star8], 
                lw =2, bg="white", # ylims=(minimum(y[:,1:3]), maximum(y[:,1:3])),
                xlab = "Index of scenarios", 
                ylab = "Total Cost",
                ylim = [-500,13000]
                )
# hline!(mean(y[:,1:3], dims=1), label = ["meanWS" "meanSO" "meanDET"], line=(2, [:solid :solid :solid], 0.5, [:blue :orange :green]))
# https://docs.juliahub.com/UnitfulRecipes/KPSlU/1.0.0/examples/2_Plots/#
savefig(p, "/Users/aaron/WildFire_Power/src/Experiments/benchmark/benchmark.pdf")


