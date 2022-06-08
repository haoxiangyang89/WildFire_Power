include("src/benchmark/optimalShutoff.jl")
include("src/benchmark/benchmark1_trivial.jl")
include("src/benchmark/benchmark2_oracle.jl")

costShutOff = optimalShutOff!();
costTrivial = benchmarkTrivial!();
costOracle = benchmarkOracle!();

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


x = comparsionResult[:,1]; oracle = comparsionResult[:,2]; shutoff = comparsionResult[:,3]; trivial = comparsionResult[:,5]; diff = comparsionResult[:,7]; # These are the plotting data
y = hcat(oracle, shutoff, trivial, diff)

Plots.plot(x, diff,  
                label = ["diff"], 
                ls=[ :solid], 
                shape=[:x], 
                lw =2, bg="white", 
                xlab = "Scenario", 
                ylab = "diff"
                )

Plots.plot(x, y[:,1:3],  
                label = ["Wait & See" "Shut Off" "No ShutOff"], 
                ls=[:dot :solid :dot ], 
                shape=[:star5 :circle :star8], 
                lw =2, bg="white", ylims=(6e7, 9.5e7),
                xlab = "Scenario", 
                ylab = "Cost"
                )


p1 = Plots.plot(x, y,  
                label = ["Wait & See" "Shut Off" "Trivial" "diff"], 
                ls=[:dash :solid :dash :solid], 
                shape=[:star5 :circle :star8 :x], 
                lw =2, bg="white", ylims=(6e7, 1.12e8),
                xlab = "Scenario", 
                ylab = "Cost"
                )
p2 = Plots.plot(x, diff,  
                shape=[:x], 
                lw =2, bg="white", 
                xlab = "Scenario", 
                ylab = "Progress",
                linecolor = [:black]
                )


Plots.plot(p1, p2, layout = (1, 2), legend = false)



Plots.plot(x, diff,  
                label = ["Wait & See" "Shut Off" "Trivial" "diff"], 
                ls=[:dash :solid :dash :solid], 
                shape=[:star5 :circle :star8 :x], 
                lw =2, bg="white", ylims=(6e7, 1.12e8),
                xlab = "Scenario", 
                ylab = "Cost"
                )
