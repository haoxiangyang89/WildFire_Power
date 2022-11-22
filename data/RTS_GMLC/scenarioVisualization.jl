using Pkg
Pkg.activate(".")
using JuMP
using JLD2, FileIO
using PowerModels
using PowerPlots


include("src/alg/def.jl")
@time Ω_rv = load("testData_RTS/Ω_rv5000.jld2")["Ω_rv"]


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

        # if k==(312,323)
        #     println(risk[t][:line][k])
        # end
        risk[t][:line][k]+=v
        risk[t][:line][k]+=length(unique(rv.Ill[k]))+length(unique(rv.Ilg[k]))+length(unique(rv.Ilb[k]))
        # if k==(312,323)
        #     println(risk[t][:line][k])
        # end
    end
end




# b_ignitions |> println
# g_ignitions |> println
# l_ignitions |> println

# using VegaLite

# @vlplot(
#     :circle,
#     x={collect(keys(b_ignitions)), type="n"},
#     y=collect(values(b_ignitions))
# )

# @vlplot(
#     :circle,
#     x={collect(keys(g_ignitions)), type="n"},
#     y=collect(values(g_ignitions))
# )

# @vlplot(
#     :circle,
#     x={collect(keys(l_ignitions)), type="n"},
#     y=collect(values(l_ignitions))
# )


# data = parse_file("RTS_GMLC.m")
network_data = parse_file("data/RTS_GMLC/case_RTS_GMLC.m")
network_data["dcline"]=Dict{String,Any}()

mn_data = replicate(network_data, 24)

for (nwid,nw) in mn_data["nw"]
    t = parse(Int,nwid)
    for (id,branch) in nw["branch"]
        f_idx=branch["f_bus"]
        t_idx=branch["t_bus"]
        branch["power_risk"] = risk[t][:line][(f_idx,t_idx)]
    end
    for (id,gen) in nw["gen"]
        gen["power_risk"] = risk[t][:gen][parse(Int,id)]
    end
    for (id,bus) in nw["bus"]
        bus["power_risk"] =  risk[t][:bus][parse(Int,id)]
    end
end

risk_max = maximum([maximum([maximum(values(risk[t][type])) for t in periods]) for type in [:line,:bus,:gen]])

using PowerPlots, ColorSchemes
function power_risk_plot(data)
    cs=reverse(colorscheme2array(ColorSchemes.colorschemes[:RdYlGn_10]))
    data = layout_network(data, edge_types=["branch"], node_types=["bus","gen","load"])
    p=powerplot(data,
    width=250,height=250,show_flow=false, fixed=true, components=["bus","branch"],
    branch_data="power_risk", branch_color=cs, branch_data_type="quantitative",
    gen_color=["#d5d5d5"], connector_color="#d5d5d5",
    bus_color=["#d5d5d5"],
    node_size=12.5, branch_size=2, connector_size=0.5)
    p.layer[1]["layer"][1]["encoding"]["color"]["scale"]["domainMax"]=maximum(collect(values(l_ignitions)))
    p.layer[1]["layer"][1]["encoding"]["color"]["legend"]=Dict("gradientLength"=>100,"title"=>"Line Risk")
    return p
end

function power_risk_mn_plot(mn_data)
    cs=reverse(colorscheme2array(ColorSchemes.colorschemes[:RdYlGn_10]))

    net_layout = layout_network(mn_data["nw"]["1"], edge_types=["branch"], node_types=["bus","gen","load"])
    for (nwid,nw) in mn_data["nw"]
        for comptype in ["bus","gen","load"]
            for (compid,comp) in nw[comptype]
                comp["xcoord_1"] = net_layout[comptype][compid]["xcoord_1"]
                comp["ycoord_1"] = net_layout[comptype][compid]["ycoord_1"]
            end
        end
        for comptype in ["branch"]
            for (compid,comp) in nw[comptype]
                if haskey(net_layout[comptype][compid], "xcoord_1") # don't deal with parallel lines
                    comp["xcoord_1"] = net_layout[comptype][compid]["xcoord_1"]
                    comp["ycoord_1"] = net_layout[comptype][compid]["ycoord_1"]
                    comp["xcoord_2"] = net_layout[comptype][compid]["xcoord_2"]
                    comp["ycoord_2"] = net_layout[comptype][compid]["ycoord_2"]
                end
            end
        end
    end


    p=powerplot(mn_data,
                width=250,height=250,show_flow=false,
                fixed=true, components=["bus","branch"],
                branch_data="power_risk", branch_color=cs, branch_data_type="quantitative",
                gen_data="power_risk", gen_color=cs, gen_data_type="quantitative",
                bus_data="power_risk", bus_color=cs, bus_data_type="quantitative",
                # gen_color=["#d5d5d5"],
                # bus_color=["#d5d5d5"],
                connector_color="#d5d5d5",
                load_color=["#d5d5d5"],
                node_size=12.5, branch_size=2, connector_size=0.5
                )
    
    p.layer[1]["encoding"]["color"]["scale"]["domainMax"]=risk_max
    p.layer[3]["encoding"]["color"]["scale"]["domainMax"]=risk_max
    p.layer[4]["encoding"]["color"]["scale"]["domainMax"]=risk_max
    p.layer[1]["encoding"]["color"]["scale"]["domainMin"]=0
    p.layer[3]["encoding"]["color"]["scale"]["domainMin"]=0
    p.layer[4]["encoding"]["color"]["scale"]["domainMin"]=0
    p.layer[1]["encoding"]["color"]["legend"]=Dict("gradientLength"=>100,"title"=>"Line Risk")
    p.layer[3]["encoding"]["color"]["legend"]=Dict("gradientLength"=>100,"title"=>"Bus Risk")
    p.layer[4]["encoding"]["color"]["legend"]=Dict("gradientLength"=>100,"title"=>"Gen Risk")
    return p
end

power_risk_mn_plot(mn_data)

using JSON
open("mn_risk_rtsgmlc.json","w") do f
    JSON.print(f, mn_data)
end

ignitions = 0
size=Int[]

for (i,rv) in Ω_rv
    for (k,v) in rv.ul
        if k==(312,323)
            ignitions += v

            a=v
            a+=length(rv.Ill[k])+length(rv.Ilg[k])+length(rv.Ilb[k])
            push!(size, a)
        end
    end
end
ignitions
println(size)