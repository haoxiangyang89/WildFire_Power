# load packages 
using Pkg
Pkg.activate(".")
using JuMP
using JLD2, FileIO
using PowerModels
using PowerPlots, ColorSchemes, VegaLite

include("src/alg/def.jl")
network_data = parse_file("data/RTS_GMLC/case_RTS_GMLC.m")

indexSets = load("testData_RTS_Sparse/indexSets.jld2")["indexSets"]
paramDemand = load("testData_RTS_Sparse/paramDemand.jld2")["paramDemand"]
Ω_rv = load("testData_RTS_Sparse/Ω_rv.jld2")["Ω_rv"]

# indexSets = load("testData_RTS_New/indexSets.jld2")["indexSets"]
# paramDemand = load("testData_RTS_New/paramDemand.jld2")["paramDemand"]
# ignitionList = load("testData_RTS_New/ignitionList.jld2")["ignitionList"]
# Ω_rv = ignitionList[50, 1]

Solution = load("testData_RTS_Sparse/Solution.jld2")["Solution"];
solution = Solution.first_state_variable;
demandSatisfication = Solution.x;

solution = modelInformationAfterShutOff.first_state_variable;
demandSatisfication = modelInformationAfterShutOff.x;
# gurobiResultList = load("testData_RTS_New/Experiments/ProbTest/gurobiResultList.jld2")["gurobiResultList"];
# Solution1 = gurobiResultList[0.1, 1]; Solution2 = gurobiResultList[0.5, 1]; 
# solution1 = Solution1.first_state_variable; demandSatisfication1 = Solution1.x;
# solution2 = Solution2.first_state_variable; demandSatisfication2 = Solution2.x;




# scenarioProcessing function
function scenarioProcessing(; Ω_rv::Dict{Int64, RandomVariables} = Ω_rv, 
                                paramDemand::ParamDemand = paramDemand, 
                                    network_data::Dict{String, Any} = network_data, 
                                    indexSets::IndexSets = indexSets,
                                    demandSatisfication = demandSatisfication,
                                    solution::Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 2, Ax, L} where {Ax, L<:Tuple{JuMP.Containers._AxisLookup, JuMP.Containers._AxisLookup}}} = solution
                                    )
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

    for key in keys(network_data["gen"]) 
        g = network_data["gen"][key]["index"]
        if paramDemand.cg[g] ≤ 50 
            network_data["gen"][key]["gen_type"] = "Wind"
        elseif paramDemand.cg[g] ≤ 1000 
            network_data["gen"][key]["gen_type"] = "Coal"
        elseif paramDemand.cg[g] ≤ 2500 
            network_data["gen"][key]["gen_type"] = "Nuclear"
        end
    end


    for key in keys(network_data["bus"]) 
        b = network_data["bus"][key]["bus_i"]
        for d in indexSets.Dᵢ[b] 
            if paramDemand.w[d] ≥ 500 
                network_data["bus"][key]["load_type"] = "High"
            elseif 200 ≤ paramDemand.w[d] ≤ 450 
                network_data["bus"][key]["load_type"] = "Middle"
            elseif 50 ≤ paramDemand.w[d] ≤ 150 
                network_data["bus"][key]["load_type"] = "Low"
            end
        end
    end

    network_data["dcline"]=Dict{String,Any}()

    mn_data = replicate(network_data, 24)

    for (nwid,nw) in mn_data["nw"]
        t = parse(Int,nwid)
        for (id,branch) in nw["branch"]
            f_idx=branch["f_bus"]
            t_idx=branch["t_bus"]
            branch["power_risk"] = risk[t][:line][(f_idx,t_idx)]
            if t > 1
                branch["ShutOff"] = solution[:zl][(f_idx, t_idx), t - 1] - solution[:zl][(f_idx, t_idx), t]
            else
                branch["ShutOff"] = 1 - solution[:zl][(f_idx, t_idx), t]
            end
            if solution[:zl][(f_idx, t_idx), t] ≥ 0.5 
                branch["Status"] = "Energized"
            else 
                branch["Status"] = "De-energized"
            end
            # branch["ShutOff"] = solution[:zl][(f_idx, t_idx), t] ≈ 1 ? "Energized" : "De-energized"
        end
        for (id,gen) in nw["gen"]
            gen["power_risk"] = risk[t][:gen][parse(Int,id)]
            g = parse(Int,id)
            if t > 1
                gen["ShutOff"] = solution[:zg][g, t - 1] - solution[:zg][g, t]
            else
                gen["ShutOff"] = 1 - solution[:zg][g, t]
            end
            if solution[:zg][g, t] ≥ 0.5 
                gen["Status"] = "Energized"
            else 
                gen["Status"] = "De-energized"
            end
            # gen["ShutOff"] = solution[:zg][g, t] ≈ 1 ? "Energized" : "De-energized"
        end
        for (id,bus) in nw["bus"]
            bus["power_risk"] =  risk[t][:bus][parse(Int,id)]
            b = parse(Int,id)
            if t > 1
                bus["ShutOff"] = solution[:zb][b, t - 1] - solution[:zb][b, t]
            else 
                bus["ShutOff"] = 1 - solution[:zb][b, t]
            end
            if solution[:zb][b, t] ≥ 0.5 
                bus["Status"] = "Energized"
            else 
                bus["Status"] = "De-energized"
            end
            # bus["ShutOff"] = solution[:zb][b, t] ≈ 1 ? "Energized" : "De-energized"

        end

        for (id,load) in nw["load"]
            d = parse(Int,id)
            b = load["load_bus"]
            bus = mn_data["nw"][nwid]["bus"]["$b"]
            bus["Load_Shed"] = round(1 - demandSatisfication[d, t], digits = 4)
            load["Load_Shed"] = round(1 - demandSatisfication[d, t], digits = 4)
        end
    end

    # risk_max = max([maximum([maximum(values(risk[t][type])) for t in periods]) for type in [:line,:bus,:gen]])
    risk_max = Dict()
    for type in [:line,:bus,:gen]
        risk_max[type] = maximum([maximum(values(risk[t][type])) for t in periods])
    end
    return (risk_max = risk_max, mn_data = mn_data)
end 


# plot risk function 
function power_risk_mn_plot(; mn_data::Dict{String, Any} = mn_data, risk_max::Dict{Any, Any} = risk_max)
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
                width=250,height=250,
                show_flow=false,
                fixed=true, components=["bus","branch"],
                branch_data="power_risk", branch_color=cs, branch_data_type="quantitative",
                gen_data="power_risk", gen_color=cs, gen_data_type="quantitative",
                bus_data="power_risk", bus_color=cs, bus_data_type="quantitative",
                # load_data="power_risk", load_color=cs, load_data_type="quantitative",
                # gen_color=["#d5d5d5"],
                # bus_color=["#d5d5d5"],
                connector_color="#d5d5d5",
                load_color=["#d5d5d5"],
                # bus_size = 25, gen_size = 20, branch_size = 3, connector_size = 1
                node_size=12.5, branch_size=2, connector_size=0.5
                )
    
    p.layer[1]["encoding"]["color"]["scale"]["domainMax"]=risk_max[:line]
    p.layer[3]["encoding"]["color"]["scale"]["domainMax"]=risk_max[:bus]
    p.layer[4]["encoding"]["color"]["scale"]["domainMax"]=risk_max[:gen]
    p.layer[1]["encoding"]["color"]["scale"]["domainMin"]=0
    p.layer[3]["encoding"]["color"]["scale"]["domainMin"]=0
    p.layer[4]["encoding"]["color"]["scale"]["domainMin"]=0
    p.layer[1]["encoding"]["color"]["legend"]=Dict("gradientLength"=>100,"title"=>"Line Risk")
    p.layer[3]["encoding"]["color"]["legend"]=Dict("gradientLength"=>100,"title"=>"Bus Risk")
    p.layer[4]["encoding"]["color"]["legend"]=Dict("gradientLength"=>100,"title"=>"Gen Risk")
    return p
end

## plot shut off function
function power_shutoff_mn_plot(; mn_data::Dict{String, Any} = mn_data, AllComponents::Bool = false)
    cs=reverse(colorscheme2array(ColorSchemes.colorschemes[:RdYlGn_10]))

    net_layout = layout_network(mn_data["nw"]["1"], edge_types=["branch"], node_types=["bus","gen","load"])
    # put the coordinate of net_layout into data_mn
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

    if AllComponents
        p=powerplot(mn_data,
                width=250,height=250,
                fixed=true, components=["bus","branch"],
                branch_data= "Status", branch_color=[:lightgrey, :green], branch_data_type="nominal", 
                gen_data= "Status", gen_color=[:lightgrey, "#0047AB"], gen_data_type="nominal", 
                # bus_data= "ShutOff", bus_color=[:lightgreen, :lightgrey], bus_data_type="nominal",
                bus_data= "Load_Shed", bus_color=[:lightgreen, "#DB492A"], bus_data_type="quantitative",
                connector_color="#d5d5d5",
                # bus_size = 25, gen_size = 20, branch_size = 3, connector_size = 1
                node_size=12.5, branch_size=2, connector_size=0.5
                )

        p.layer[3]["encoding"]["color"]["legend"]=Dict("gradientLength"=>100,"title"=>"Load Shed")
    else
        p=powerplot(mn_data,
                width=250,height=250,
                show_flow=false,
                fixed=true, components=["bus","branch", "load"],
                branch_data="ShutOff", branch_color=[:green, :lightgrey], branch_data_type="nominal",
                gen_data= :gen_type, gen_data_type = "nominal", gen_color = colorscheme2array(ColorSchemes.colorschemes[:seaborn_deep]), # tableau_green_blue_white
                bus_data= "load_type", bus_data_type = "nominal", bus_color = [:green, :Red, :grey, :orange],
                connector_color="#d5d5d5",
                # bus_size = 25, gen_size = 20, branch_size = 3, connector_size = 1
                node_size=12.5, branch_size=2, connector_size=0.5
                )
        p.layer[1]["encoding"]["color"]["legend"]=Dict("gradientLength"=>100,"title"=>"Branch ShutOff")
        p.layer[4]["encoding"]["color"]["legend"]=Dict("gradientLength"=>100,"title"=>"Gen Type")
        p.layer[3]["encoding"]["color"]["legend"]=Dict("gradientLength"=>100,"title"=>"Load Type")
        
    end
    # p.layer[1]["encoding"]["color"]["scale"]["domainMax"]=risk_max
    # p.layer[5]["encoding"]["color"]["scale"]["domainMax"]=risk_max
    # p.layer[4]["encoding"]["color"]["scale"]["domainMax"]=risk_max
    # p.layer[1]["encoding"]["color"]["scale"]["domainMin"]=0
    # p.layer[5]["encoding"]["color"]["scale"]["domainMin"]=0
    # p.layer[4]["encoding"]["color"]["scale"]["domainMin"]=0
    # p.layer[1]["encoding"]["color"]["legend"]=Dict("gradientLength"=>100,"title"=>"Branch")
    # p.layer[5]["encoding"]["color"]["legend"]=Dict("gradientLength"=>100,"title"=>"Gen")
    # p.layer[4]["encoding"]["color"]["legend"]=Dict("gradientLength"=>100,"title"=>"Bus")
    return p
end


(risk_max, mn_data) = scenarioProcessing(; Ω_rv = Ω_rv, solution = solution, network_data = network_data);
# (risk_max, mn_data1) = scenarioProcessing(; Ω_rv = Ω_rv, solution = solution1, demandSatisfication = demandSatisfication1, network_data = network_data);
# (risk_max, mn_data2) = scenarioProcessing(; Ω_rv = Ω_rv, solution = solution2, demandSatisfication = demandSatisfication2, network_data = network_data);


p1 = power_risk_mn_plot(; mn_data = mn_data)
p2 = power_shutoff_mn_plot(; mn_data = mn_data, AllComponents = true)
p3 = power_shutoff_mn_plot(; mn_data = mn_data, AllComponents = false)


# p1 = power_risk_mn_plot(; mn_data = mn_data1)
# p2 = power_shutoff_mn_plot(; mn_data = mn_data1, AllComponents = true)
# p3 = power_shutoff_mn_plot(; mn_data = mn_data1, AllComponents = false)

# p12 = power_risk_mn_plot(; mn_data = mn_data2)
# p22 = power_shutoff_mn_plot(; mn_data = mn_data2, AllComponents = true)
# p32 = power_shutoff_mn_plot(; mn_data = mn_data2, AllComponents = false)






# indexSets = load("testData_RTS_Sparse/indexSets.jld2")["indexSets"]
# paramOPF = load("testData_RTS_Sparse/paramOPF.jld2")["paramOPF"]
# paramDemand = load("testData_RTS_Sparse/paramDemand.jld2")["paramDemand"]
# Ω_rv = load("testData_RTS_Sparse/Ω_rv.jld2")["Ω_rv"]
# prob = load("testData_RTS_Sparse/prob.jld2")["prob"]


# modelInformation = gurobiOptimizeTest!(indexSets, 
#                                           paramDemand, 
#                                           paramOPF, 
#                                           Ω_rv,
#                                           prob; 
#                                           mipGap = 1e-2, timelimit = 6000)


# save("testData_RTS_Sparse/modelInformationAfterShutOff.jld2", "modelInformationAfterShutOff", modelInformationAfterShutOff)

