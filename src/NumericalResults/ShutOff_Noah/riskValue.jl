# load packages 
using Pkg
Pkg.activate(".")
using JuMP
using JLD2, FileIO
using PowerModels
using PowerPlots, ColorSchemes, VegaLite

include("src/alg/def.jl")
network_data = parse_file("data/RTS_GMLC/case_RTS_GMLC.m")

indexSets = load("data/testData_RTS/demand/indexSets.jld2")["indexSets"]
paramDemand = load("data/testData_RTS/demand/paramDemand.jld2")["paramDemand"]
Ω_rv = load("data/testData_RTS/Ω_rvList.jld2")["Ω_rvList"][500, 13]


periods = 1:24
b_index = keys(Ω_rv[1].ub)
g_index = keys(Ω_rv[1].ug)
l_index = keys(Ω_rv[1].ul)

risk = Dict{Int,Any}(i=>Dict{Symbol,Any}() for i in periods)
for i in periods
    risk[i][:bus] = Dict{Int,Int}(i=>0 for i in b_index)
    risk[i][:gen] = Dict{Int,Int}(i=>0 for i in g_index)
    risk[i][:line] = Dict{Tuple{Int,Int},Int}(i=>0 for i in l_index)
end


for (i,rv) in Ω_rv
    t = rv.τ
    for (k,v) in rv.ub
        risk[t][:bus][k]+=v
        risk[t][:bus][k]+=length(unique(rv.Ibl[k]))+length(unique(rv.Ibg[k]))+length(unique(rv.Ibb[k]))
    end
    for (k,v) in rv.ug
        risk[t][:gen][k]+=v
        risk[t][:gen][k]+=length(unique(rv.Igl[k]))+length(unique(rv.Igg[k]))+length(unique(rv.Igb[k]))
    end
    for (k,v) in rv.ul
        risk[t][:line][k]+=v
        risk[t][:line][k]+=length(unique(rv.Ill[k]))+length(unique(rv.Ilg[k]))+length(unique(rv.Ilb[k]))
    end
end



save("src/NumericalResults/ShutOff_Noah/risk.jld2", "risk", risk)



