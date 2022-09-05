using CSV, DataFrames, Printf, Gurobi, JuMP
using JLD2, FileIO

totalCost = load("src/Experiments/ScenariosSizeTest/totalCost.jld2")["totalCost"]
gurobiResultList = load("src/Experiments/ScenariosSizeTest/gurobiResultList.jld2")["gurobiResultList"]



LBDF = Dict()
for (i, j) in keys(gurobiResultList)
    LBDF[j, i] = gurobiResultList[i,j].OPT
end

gapdf = Dict()
for (k, ss) in keys(LBDF)
    if ss != 150
        gapdf[k, ss] = totalCost[k, ss] - LBDF[k, ss]
    end
end
## ------------------------------ Transfer into two dataframes ----------------------------- ##

costMatrix = zeros(20,4)
for (i,j) in keys(gapdf)
  index = 0
  for k in [20, 50, 100, 200]
    if k <= j 
      index = index + 1
    end
  end
  costMatrix[i,index] = gapdf[(i,j)]
end
gapDF = DataFrame(costMatrix, :auto)
gapDF = DataFrames.select(gapDF, :x1 => :sz20, :x2 => :sz50, :x3 => :sz100, :x4 => :sz200)



# https://discourse.julialang.org/t/how-to-plot-mean-deviation-in-plots-jl/1975/15
using Statistics
using Plots; pyplot()
xs = [20, 50, 100, 200, 500]
μs1, σs1 = mean.(eachcol(gapDF)), std.(eachcol(gapDF))
# plot ribbon
Plots.plot( xs, μs1, color=:lightblue, ribbon=σs1,label=false)
# plot mean point
Plots.plot!(xs, μs1, color=:blue, marker=(:circle, 8, 1.), label="Optimality Gap",  xlab = "Sample size", ylab = "Value of Gap")

# plot CI
n = 5; ## the number of columns
YUB = [gapDF[:, i] for i in 1:n]
Ym = mean.(YUB)
ϵ⁻ = Ym .- quantile.(YUB, fill(0.25, n))
ϵ⁺ = quantile.(YUB, fill(0.75, n)) .- Ym
scatter!(xs, Ym, ms=6, yerror=(ϵ⁻, ϵ⁺), label= false, title = "Confidence intervals and point estimates of optimality gaps")


# -------------------- 
# plot ribbon
Plots.scatter( xs, μs1, color=:lightblue,label=false)

# plot mean point

# plot CI
n = 4; ## the number of columns
YUB = [gapDF[:, i] for i in 1:n]
Ym = mean.(YUB)          
ϵ⁻ = Ym .- quantile.(YUB, fill(0.25, n))
ϵ⁺ = quantile.(YUB, fill(0.75, n)) .- Ym                  # 106.9, [25.32498985, 31.32498985]
scatter!(xs, Ym, ms=6, yerror=(ϵ⁻, ϵ⁺), label=false)
# --------------------


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

t_test(gapDF[:,1])
t_test(gapDF[:,2])
t_test(gapDF[:,3])
t_test(gapDF[:,4])
t_test(costDF[:,5])




CI_DF = DataFrame(:ss => [20, 50, 100, 200, 500], :GAP => Ym, :CI => [[831.3, 1081.2], [491.6, 671.1], [365.8, 489.5], [267.0, 371.9], [81.6, 138.3]], :Percentage => (ϵ⁺ .+ ϵ⁻)./2642.63 * 100)
# latexify(CI_DF) ## 2642.63

# [[831.3, 1081.2], [491.6, 671.1], [365.8, 489.5], [267.0, 371.9], [81.6, 138.3]]

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
