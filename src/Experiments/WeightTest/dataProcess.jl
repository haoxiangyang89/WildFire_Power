using CSV, DataFrames, Printf, Gurobi, JuMP
using JLD2, FileIO

totalCost = load("src/Experiments/WeightTest/totalCostNew.jld2")["totalCost"]
gurobiResultList = load("src/Experiments/WeightTest/gurobiResultListNew.jld2")["gurobiResultList"]

LBDF = Dict()
for (i, j) in keys(gurobiResultList)
  LBDF[j, i] = gurobiResultList[i,j].OPT
end


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
  for k in [0.0, 0.01, 0.05, .5, .9, .95, .99, 1.0]
    if k <= j 
      index = index + 1
    end
  end
  costMatrix[i,index] = totalCost[(i,j)]
end
costUBDF = DataFrame(costMatrix, :auto)
# costUBDF = DataFrames.select(costUBDF, :x1 => :prob0, :x2 => :prob3, :x3 => :prob6)
costUBDF = DataFrames.select(costUBDF, :x1 => :prob0, :x2 => :prob01, :x3 => :prob05, :x4 => :prob5, :x5 => :prob9, :x6 => :prob95, :x7 => :prob99, :x8 => :prob100)



## second one for LB
costMatrix = zeros(20, 8)
for (i,j) in keys(LBDF)
  index = 0
  for k in [0.0, 0.01, 0.05, .5, .9, .95, .99, 1.0]
    if k <= j 
      index = index + 1
    end
  end
  costMatrix[i,index] = LBDF[(i,j)]
end
costLBDF = DataFrame(costMatrix, :auto)
# costLBDF = DataFrames.select(costLBDF, :x1 => :prob0, :x2 => :prob3, :x3 => :prob6)
costLBDF = DataFrames.select(costLBDF, :x1 => :prob0, :x2 => :prob01, :x3 => :prob05, :x4 => :prob5, :x5 => :prob9, :x6 => :prob95, :x7 => :prob99, :x8 => :prob100)


using StatsPlots, DataFrames
X = 0:2.5:17.5 # range(1,10, step = 2.5)
n = length(X)
Y = [costUBDF[:, i ] for i in 1:8]
X2 = [fill(x,length(y)) for (x,y) in zip(X,Y)]
df = DataFrame(X = X2, Y = Y)
# labels = reshape(["Ignition Prob 1.0", "Ignition Prob 0.5", "Ignition Prob 0.1", "Ignition Prob 0.05", "Ignition Prob 0.01"], 1, 5)
labels = reshape([1.0, 0.99, 0.95, 0.5, 0.1, 0.05, 0.01, 0.0], 1, 8)
@df df StatsPlots.boxplot(:X, :Y, legend=:topleft, label = labels, xticks = false, ylab = "Total Cost", title = "Robustness of the first-stage decision")

X = 0:2.5:17.5 # range(1,10, step = .5)
n = length(X)
Y = [costLBDF[:, i ] for i in 1:8]
X2 = [fill(x,length(y)) for (x,y) in zip(X,Y)]
df = DataFrame(X = X2, Y = Y)
@df df StatsPlots.boxplot(:X, :Y, legend=false)













# https://discourse.julialang.org/t/how-to-plot-mean-deviation-in-plots-jl/1975/15
using Statistics
using Plots; pyplot()
xs =  [1.0, .7, .4, 0.0]
μs1, σs1 = mean.(eachcol(costUBDF)),   std.(eachcol(costUBDF))
μs2, σs2 = mean.(eachcol(costLBDF)),   std.(eachcol(costLBDF))
# plot ribbon
Plots.plot( xs, μs1, color=:lightblue, ribbon=σs1, label=false)
Plots.plot!( xs, μs2, color=:pink, ribbon=σs1, label=false)
# plot mean point
Plots.plot!(xs, μs1, color=:blue, marker=(:circle, 8, 1.), label="Upper Bound")
Plots.plot!(xs, μs2, color=:red, marker=(:circle, 8, 1.), label="LB")
# plot CI
n = 4; 
YLB = [costLBDF[:, i] for i in 1:n]
Ym = mean.(YLB)
ϵ⁻ = Ym .- quantile.(YLB, fill(0.25, n))
ϵ⁺ = quantile.(YLB, fill(0.75, n)) .- Ym
scatter!(xs, Ym, ms=6, yerror=(ϵ⁻, ϵ⁺), label=false)

YUB = [costUBDF[:, i] for i in 1:n]
Ym = mean.(YUB)
ϵ⁻ = Ym .- quantile.(YUB, fill(0.25, n))
ϵ⁺ = quantile.(YUB, fill(0.75, n)) .- Ym
scatter!(xs, Ym, ms=6, yerror=(ϵ⁻, ϵ⁺), label=false, ylim = (2650, 2950))





# ---------------------------- 
xs = [0.0, 0.3, 0.6, 0.9]
μs1, σs1 = mean.(eachcol(costUBDF)),   std.(eachcol(costUBDF))
μs2, σs2 = mean.(eachcol(costLBDF)),   std.(eachcol(costLBDF))
# plot ribbon
Plots.scatter( xs, μs1, color=:lightblue, label=false)
Plots.scatter!( xs, μs2, color=:pink, label=false)
# plot mean point
Plots.scatter!(xs, μs1, color=:blue, marker=(:circle, 8, 1.), label="UB")
Plots.scatter!(xs, μs2, color=:red, marker=(:circle, 8, 1.), label="LB")
# plot CI
n = 4; 
YLB = [costLBDF[:, i] for i in 1:n]
Ym = median.(YLB)
ϵ⁻ = Ym .- quantile.(YLB, fill(0.25, n))
ϵ⁺ = quantile.(YLB, fill(0.75, n)) .- Ym
scatter!(xs, Ym, ms=6, yerror=(ϵ⁻, ϵ⁺), label=false)

savefig("src/Experiments/WeightTest/WeightTest.png")
## ------------------------------- Data Processing ------------------------------ #
## compute CI
using Distributions
function t_test(x; conf_level=0.95)
    alpha = (1 - conf_level)
    tstar = quantile(TDist(length(x)-1), 1 - alpha/2)
    SE = std(x)/sqrt(length(x))

    lo, hi = mean(x) .+ [-1, 1] * tstar * SE
    "($lo, $hi)"
end

t_test(costUBDF[:,1])
t_test(costUBDF[:,2])
t_test(costUBDF[:,3])
t_test(costUBDF[:,4])
t_test(costUBDF[:,5])
t_test(costUBDF[:,6])
t_test(costUBDF[:,7])
t_test(costUBDF[:,8])





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
