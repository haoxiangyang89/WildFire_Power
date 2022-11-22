using CSV, DataFrames, Printf, Gurobi, JuMP
using JLD2, FileIO

totalCost = load("testData_RTS_New/Experiments/ProbTest/totalCost.jld2")["totalCost"]

## ------------------------------ Transfer into two dataframes ----------------------------- ##
# first one
## boxplot
using DataFrames, Random, StatsPlots
StatsPlots.theme(:ggplot2)
default(dpi=100)

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
p1 = @df df StatsPlots.plot(:solClass, :Cost, marker=:circle, ms=3, group=:solClass, legend=:outertopright)
p2 = @df df StatsPlots.boxplot(:solClass, :Cost, ms=3, group=:solClass, legend=:outertopright)
plotd = StatsPlots.plot(p1, p2, layout = (1, 2), xlim = (-1, 10),legend = true, ylab = "Cost", xlab = "Prob. Without Ignition")
# StatsPlots.savefig(plotd,"src/Experiments/WeightTest/WeightTestBoxplot.png")
# @df df Plots.violin(string.(:solClass), :Cost, linewidth=0)
@df df StatsPlots.boxplot(:solClass .* 10, :Cost, ms=3, group=:solClass, legend=true, 
                                ylab = "Cost", title = "Robustness of the first-stage decision", 
                                     framestyle = :box)

## second one for UB
costMatrix = zeros(20, 8)
for (i,j) in keys(totalCost)
  index = 0
  for k in [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.9, 1.0]
    if k <= j 
      index = index + 1
    end
  end
  costMatrix[i,index] = totalCost[(i,j)]
end
costUBDF = DataFrame(costMatrix, :auto)
# costUBDF = DataFrames.select(costUBDF, :x1 => :prob0, :x2 => :prob3, :x3 => :prob6)
costUBDF = DataFrames.select(costUBDF, :x1 => :prob01, :x2 => :prob05, :x3 => :prob1, :x4 => :prob2, 
                                          :x5 => :prob3, :x6 => :prob5, :x7 => :prob9, :x8 => :prob10)



using StatsPlots, DataFrames
X = 0:2.5:17.5 # range(1,10, step = 2.5)
n = length(X)
Y = [costUBDF[:, i ] for i in 1:7]
X2 = [fill(x,length(y)) for (x,y) in zip(X,Y)]
df = DataFrame(X = X2, Y = Y)
# labels = reshape(["Ignition Prob 1.0", "Ignition Prob 0.5", "Ignition Prob 0.1", "Ignition Prob 0.05", "Ignition Prob 0.01"], 1, 5)
labels = reshape([0.0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5], 1, 7)
@df df StatsPlots.boxplot(:X, :Y, legend=:outertopright, label = labels, xticks = false, 
                    ylab = "Total Cost", title = "Influence of Distributions", ylim = [4700, 7300])

savefig("testData_RTS_New/Experiments/ProbTest/WeightTest.pdf")







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
