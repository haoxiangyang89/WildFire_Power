using CSV, DataFrames, Printf, Gurobi, JuMP
using JLD2, FileIO

totalCost = load("src/Experiments/ScenariosSizeTest/totalCost.jld2")["totalCost"]
gurobiResultList = load("src/Experiments/ScenariosSizeTest/gurobiResultList.jld2")["gurobiResultList"]

UBDF = Dict()
for (i, j) in keys(gurobiResultList)
  UBDF[j, i] = gurobiResultList[i,j].OPT
end

# for i in 1:20 
#   totalCost[i, 500] = totalCost[i, 500] - 250 
#   UBDF[i, 500] = UBDF[i, 500] - 50
# end




## ------------------------------ Transfer into two dataframes ----------------------------- ##
# first one
## boxplot
using DataFrames, Random, StatsPlots
StatsPlots.theme(:ggplot2)
default(dpi=600)

iter_num = 20; sampleList = []; solClass = []; costList = [];
for (i,j) in keys(totalCost) 
  push!(sampleList, i)
  push!(solClass, j)
  push!(costList, totalCost[(i,j)])
end
df = DataFrame(solClass = solClass, 
                Cost = costList
                )
nam0 = first.(keys(DataFrames.groupby(df, :solClass)))
cdic = Dict(nam0 .=> palette(:default)[1:length(nam0)])
col0 = [cdic[x] for x in df.solClass]

# Table-like data structures, including DataFrames, IndexedTables, DataStreams, etc... (see here for an exhaustive list), 
# are supported thanks to the macro @df which allows passing columns as symbols.
# https://github.com/JuliaPlots/StatsPlots.jl
p1 = @df df StatsPlots.plot(:solClass ./10, :Cost, marker=:circle, ms=3, group=:solClass, legend=:outertopright)
p2 = @df df StatsPlots.boxplot(:solClass ./10, :Cost, ms=3, group=:solClass, legend=:outertopright)
plotd = StatsPlots.plot(p1, p2, layout = (1, 2), legend = true, ylab = "Cost")
# StatsPlots.savefig(plotd,"src/Experiments/ScenariosSizeTest/file4.png")
# @df df Plots.violin(string.(:solClass), :Cost, linewidth=0)


## second one for UB
costMatrix = zeros(20,5)
for (i,j) in keys(totalCost)
  index = 0
  for k in [20, 50, 100, 200, 500]
    if k <= j 
      index = index + 1
    end
  end
  costMatrix[i,index] = totalCost[(i,j)]
end
costUBDF = DataFrame(costMatrix, :auto)
costUBDF = DataFrames.select(costUBDF, :x1 => :sz20, :x2 => :sz50, :x3 => :sz100, :x4 => :sz200, :x5 => :sz500)



## second one for LB
costMatrix = zeros(20,5)
for (i,j) in keys(UBDF)
  index = 0
  for k in [20, 50, 100, 200, 500]
    if k <= j 
      index = index + 1
    end
  end
  costMatrix[i,index] = UBDF[(i,j)]
end
costLBDF = DataFrame(costMatrix, :auto)
costLBDF = DataFrames.select(costLBDF, :x1 => :sz20, :x2 => :sz50, :x3 => :sz100, :x4 => :sz200, :x5 => :sz500)


using StatsPlots, DataFrames
X = 0:2.5:10 # range(1,10, step = .5)
n = length(X)
Y = [costUBDF[:, i ] for i in 1:5]
X2 = [fill(x,length(y)) for (x,y) in zip(X,Y)]
df = DataFrame(X = X2, Y = Y)
@df df StatsPlots.boxplot(:X, :Y, legend=false)

X = 0:2.5:10 # range(1,10, step = .5)
n = length(X)
Y = [costLBDF[:, i ] for i in 1:5]
X2 = [fill(x,length(y)) for (x,y) in zip(X,Y)]
df = DataFrame(X = X2, Y = Y)
@df df StatsPlots.boxplot(:X, :Y, legend=false)

# https://discourse.julialang.org/t/how-to-plot-mean-deviation-in-plots-jl/1975/15
using Statistics
using Plots; pyplot()
xs = [20, 50, 100, 200, 500]
μs1, σs1, max1, min1 = mean.(eachcol(costUBDF)),   std.(eachcol(costUBDF)), maximum.(eachcol(costUBDF)), minimum.(eachcol(costUBDF))
μs2, σs2, max2, min2 = mean.(eachcol(costLBDF)),   std.(eachcol(costLBDF)), maximum.(eachcol(costLBDF)), minimum.(eachcol(costLBDF))
# plot ribbon
p = Plots.plot( xs, μs1, color=:lightblue, ribbon= (μs1 .- min1, max1 .- μs1),label=false)
Plots.plot!( xs, μs2, color=:pink, ribbon= (μs2 .- min2, max2 .- μs2),label=false)
# plot mean point
Plots.plot!(xs, μs1, color=:blue, marker=(:circle, 8, 1.), label="Upper Bound")
Plots.plot!(xs, μs2, color=:red, marker=(:circle, 8, 1.), label="Lower Bound", xlab = "Sample size", ylab = "Value of Bounds")
# plot CI
# https://discourse.julialang.org/t/asymmetric-error-bars-and-box-plots-in-julia/73647/3
n = 5; ## the number of columns
YLB = [costLBDF[:, i] for i in 1:n]
Ym = mean.(YLB)
ϵ⁻ = 1.96 .* σs2
ϵ⁺ = 1.96 .* σs2
scatter!(xs, Ym, ms=6, yerror=(ϵ⁻, ϵ⁺), label= false)

YUB = [costUBDF[:, i] for i in 1:n]
Ym = mean.(YUB)
ϵ⁻ = 1.96 .* σs1
ϵ⁺ = 1.96 .* σs1
scatter!(xs, Ym, ms=6, yerror=(ϵ⁻, ϵ⁺), label= false, title = "Confidence intervals and point estimates of bounds")
savefig(p, "/Users/aaron/WildFire_Power/src/Experiments/ScenariosSizeTest/ConfidenceInterval.pdf")



## ------------------------------- Data Processing ------------------------------ #
## compute CI
# using Distributions
# function t_test(x; conf_level=0.95)
#     alpha = (1 - conf_level)
#     tstar = quantile(TDist(length(x)-1), 1 - alpha/2)
#     SE = std(x)/sqrt(length(x))

#     lo, hi = mean(x) .+ [-1, 1] * tstar * SE
#     "($lo, $hi)"
# end

# t_test(costUBDF[:,1])
# t_test(costUBDF[:,2])
# t_test(costUBDF[:,3])
# t_test(costUBDF[:,4])
# t_test(costUBDF[:,5])



## ------------------------------- Data Processing ------------------------------ #
# using DataFrames
# using Statistics

# DF = DataFrame(ID = 1:10, Col1 = rand(15:0.01:40,10),
#                           Col2 = rand(30:0.01:55,10),
#                           Col3 = rand(20:0.01:65,10))

# round3(x) = round(x, digits=3)

# transform(DF, AsTable(2:4) .=>
#               ByRow.([minimum, maximum, round3∘var, round3∘std]) .=> [:Min, :Max, :Var, :StdDev])

# transform(DF, AsTable(:) => ByRow(argmax) => :prediction)
# var.(eachcol(DF))




# using DataFrames, Random, StatsPlots
# Plots.theme(:ggplot2)
# default(dpi=600)

# Random.seed!(123)
# df = DataFrame(name = repeat(["A","B","C","D","E","F"], inner=4), 
#                 time=repeat([0,1,3,6], outer=6), value = rand(24)
#                 )
# #
# dx, dy = extrema.((df.time, df.value))
# dx, dy = 0.025 .* (dx[2] - dx[1], dy[2] - dy[1])
# nam0 = first.(keys(DataFrames.groupby(df, :name)))
# cdic = Dict(nam0 .=> palette(:default)[1:length(nam0)])
# col0 = [cdic[x] for x in df.name]

# @df df Plots.plot(:time, :value, marker=:circle, ms=3, group=:name, legend=:outertopright, legendtitle="name")
# for (x,y,nm,c) in zip(df.time, df.value, df.name, col0)
#     annotate!(x + rand((-dx,dx)), y + rand((-dy,dy)), (nm, 7, c))
# end
# Plots.current()
