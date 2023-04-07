using CSV, DataFrames, Printf, Gurobi, JuMP
using JLD2, FileIO

totalCost = load("src/NumericalResults/computationalPerformance/totalCost.jld2")["totalCost"]
gurobiResultList = load("src/NumericalResults/computationalPerformance/gurobiResultList.jld2")["gurobiResultList"]

UBDF = Dict()
for (i, j) in keys(gurobiResultList)
  UBDF[j, i] = gurobiResultList[i,j].OPT
end


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



# https://discourse.julialang.org/t/how-to-plot-mean-deviation-in-plots-jl/1975/15
using Statistics
using Plots; pyplot()
xs = [20, 50, 100, 200, 500]
μs1, σs1, max1, min1 = mean.(eachcol(costUBDF)),   std.(eachcol(costUBDF)), maximum.(eachcol(costUBDF)), minimum.(eachcol(costUBDF))
μs2, σs2, max2, min2 = mean.(eachcol(costLBDF)),   std.(eachcol(costLBDF)), maximum.(eachcol(costLBDF)), minimum.(eachcol(costLBDF))
# plot ribbon
p = Plots.plot( xs, μs1, color=:lightblue, ribbon= (μs1 .- min1, max1 .- μs1),label=false, guidefontsize=16, xtickfontsize=17, ytickfontsize=17, legendfontsize=15)
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
savefig(p, "/Users/aaron/WildFire_Power/src/NumericalResults/computationalPerformance/ConfidenceInterval.pdf")

