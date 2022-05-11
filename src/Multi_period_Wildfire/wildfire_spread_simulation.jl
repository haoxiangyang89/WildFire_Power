################################################################################################################################################
#####################################################  Cellular Automton  ######################################################################
################################################################################################################################################
using Agents, Random, Distributions
using CSV, Geodesy
using InteractiveDynamics
using CairoMakie
using PowerModels




## data struct that given current wind info
struct WindInfo
    V   ::Float64  ## wind speed
    θw  ::Int64    ## wind direction (0 - 360)
end



mutable struct CellInfo <: AbstractAgent
    id::Int
    pos::NTuple{2, Int}
    # P₀          ::Float64   ## (0.58)the probability that a neighboring cell is burning and garment in the next step of the simulation under conditions of absence of fire and elevation difference between the central cell and neighboring
    Pveg        ::Float64   ## Vegetation density (empty = -1.0, limited/cultivated = -0.4, forests = 0.4, shrub = 0.4)
    Pden        ::Float64   ## Kind of vegetation (empty = -1.0, sparse = -.3, normal = 0.0, dense = 0.3)
    E           ::Float64   ## elevation of the cell
    state       ::Int64     ## 0. The cell does not contain vegetable fuel and, therefore, it can not burn.
                             # 1. It contains fuel (vegetation) that has not ignited.
                             # 2. It contains burning vegetation.
                             # 3. It contains vegetation that has burned completely
end



struct CellEnvironmentInfo
    Pveg        ::Float64   ## Vegetation density (empty = -1.0, limited/cultivated = -0.4, forests = 0.4, shrub = 0.4)
    Pden        ::Float64   ## Kind of vegetation (empty = -1.0, sparse = -.3, normal = 0.0, dense = 0.3)
    E           ::Float64   ## elevation of the cell
    state       ::Int64     ## 0. The cell does not contain vegetable fuel and, therefore, it can not burn.
                             # 1. It contains fuel (vegetation) that has not ignited.
                             # 2. It contains burning vegetation.
                             # 3. It contains vegetation that has burned completely

end





## function to compute the probability of fire for agent when given neighbor info
function Probability_burn(neighbor::CellInfo, 
                            agent::CellInfo;
                            windInfo::WindInfo = windInfo,
                            P₀::Float64 = .58,
                            c1::Float64 = 0.045, c2::Float64 = 0.131, 
                            a::Float64 = 0.078,
                            l::Float64 = 1000. ## the distance between those two cells
                            )

    diff_pos = agent.pos .- neighbor.pos
    if diff_pos[1] == 1
        if diff_pos[2] == 1
            θc = 315
        elseif diff_pos[2] == 0
            θc = 270
        elseif diff_pos[2] == -1
            θc = 225
        end
    elseif diff_pos[1] == 0
        if diff_pos[2] == 1
            θc = 0
        elseif diff_pos[2] == -1
            θc = 180
        end
    elseif diff_pos[1] == -1
        if diff_pos[2] == 1
            θc = 45
        elseif diff_pos[2] == 0
            θc = 90
        elseif diff_pos[2] == -1
            θc = 135
        end
    end
    
    θ = abs(θc - windInfo.θw) * π/180
    ft = exp(windInfo.V * c2 * (cos(θ) - 1))
    Pw = exp(c1 * windInfo.V) * ft
    if θc%90 == 0
        θs = atan((agent.E - neighbor.E)/l)
    else
        θs = atan((agent.E - neighbor.E)/(1.414 * l))
    end
    Pele = exp(a * θs)
    Prob = P₀ * (1 + agent.Pveg) * (1 + agent.Pden) * Pw * Pele
    return Prob
end



function initialize(relative_location::Dict{Int64, Tuple{Int64, Int64}}, line_location_id::Dict{Any, NamedTuple{(:id, :WFPI), Tuple{Int64, Float64}}} ;  griddims::Tuple{Int64, Int64} = (100, 100), 
                        environmentInfo::Dict{Tuple{Int64, Int64}, CellEnvironmentInfo} = environmentInfo, time_span::Int64 = 1)

    space = GridSpace(griddims, periodic = false)
    # The `trees` field is coded such that
    # No ignition = 0, ignition = 1, burnt = 2
    forest = ABM(
        CellInfo, space;
        properties = (
                        ignition = zeros(Int, griddims),
                        fault_WFPI = zeros(Float64, griddims),
                        bus_exist = zeros(Int, griddims),
                        line_exist = zeros(Int, griddims),
                        bus_fired = zeros(Int, griddims), 
                        line_fired = zeros(Int, griddims), 
                        bus_fault = zeros(Int, griddims), 
                        line_fault = zeros(Int, griddims), 
                        windSpeed = [0.1, ],
                        windDirection = [176, ],
                        time_span = time_span
                        )
    )

    # populate the model with agents, adding equal amount of the two types of agents
    # at random positions in the model


    id = 1
    for I in CartesianIndices(forest.ignition)
        ## need to find a logic for the following fields
        Pveg = environmentInfo[I.I].Pveg         ## Vegetation density (limited = -0.4, medium = 0., high = 0.3)
        Pden = environmentInfo[I.I].Pden         ## Kind of vegetation
        E = environmentInfo[I.I].E               ## elevation of the cell
        state = environmentInfo[I.I].state 
        agent = CellInfo(id, (1, 1), Pveg, Pden, E, state)
        add_agent!(agent, I.I, forest)
        id = id + 1

        # forest.ignition[I] = min(state, 0) ## No ignition = 0, ignition = 1, burnt = 2
    end

    for bus_id in keys(relative_location)
        pos = relative_location[bus_id]
        forest.bus_exist[pos...] = 1
    end

    for pos in keys(line_location_id)
        forest.line_exist[pos...] = 1
        forest.fault_WFPI[pos...] = line_location_id[pos].WFPI/100
    end
    return forest
end



function agent_step!(agent, forest)
    windInfo = WindInfo(forest.windSpeed[1], forest.windDirection[1])
    time_span = forest.time_span
    if agent.state == 2
        ## If it is burning previously, then it is burnt now
        agent.state = rand(forest.rng) <= time_span/24 ? 3 : agent.state
    elseif agent.state == 1
        ## If this cell has an ignition, then its state becomes burning
        if forest.ignition[agent.pos...] == 1
            # agent.state = 2
            agent.state = rand(forest.rng) <= time_span/24 ? 2 : 1
        else
            ## If nearby agent is burnt, then its state becomes burning
            for neighbor in nearby_agents(agent, forest)
                if neighbor.state == 2
                    p = Probability_burn(neighbor, agent; windInfo = windInfo, P₀ = .58)
                    agent.state = rand(1)[1] <= p ? 2 : agent.state
                end
            end
        end
    end

    return 
end




function wildfire_ignition_step!(forest::AgentBasedModel)
        # ignition will only occure at the cells that contain fuel (vegetation) but has not ignited w.p. p
        for id in 1:nagents(forest)
            agent = forest[id]
            if (agent.state == 0 || agent.state == 1)
                forest.ignition[agent.pos...] = 0

                # p = Probability_ignition(agent, windInfo)
                p = 1e-4
                forest.ignition[agent.pos...] = (rand(forest.rng) <= p && agent.state == 1) ? 1 : 0
            else 
                forest.ignition[agent.pos...] = agent.state - 1
            end
        end

        for I in findall(isequal(1), forest.bus_exist)
            if forest.ignition[I] == 1  ## if the cell with one bus is burning, we say this bus is fired
                forest.bus_fault[I] = rand(forest.rng) <= 0.1 ? 1 : 0
            end

            if forest.ignition[I] == 2  ## if the cell with one bus is burnt, we say this bus is fired
                forest.bus_fired[I] = 1
            end
        end


        for I in findall(isequal(1), forest.line_exist)
            if forest.ignition[I] == 1  ## if the cell with one bus is burning, we say this bus is fired
                forest.line_fault[I] = rand(forest.rng) <= forest.fault_WFPI[I] ? 1 : 0
            end

            if forest.ignition[I] == 2  ## if the cell with one branch burnt, we say this line is fired
                forest.line_fired[I] = 1
            end
        end
        forest.windSpeed[1] = wsample([.1, .2, .4, .7], [0.2, .4, .3, 0.1], 1)[1]
        forest.windDirection[1] = wsample([10, 80, 120, 190, 250, 300, 340], [0.1, .2, .3, 0.1, 0.1, 0.1, .1], 1)[1]

end



################################################################################################################################################
######################################################  Data Preparation  ######################################################################
################################################################################################################################################




function prepareIndexSets(  network_data::Dict{String, Any} ,
                                                    T::Int64,
                                                    Ω::Int64
                                                    )
    

    D = Vector{Int64}()
    G = Vector{Int64}()
    B = Vector{Int64}()

    L = Vector{Tuple{Int64, Int64}}()

    Dᵢ =    Dict{Int64,Vector{Int64}}()
    Gᵢ =    Dict{Int64,Vector{Int64}}()
    out_L = Dict{Int64,Vector{Int64}}()
    in_L =  Dict{Int64,Vector{Int64}}()


    _b = Dict{Tuple{Int64, Int64}, Float64}()  ## total line charging susceptance
    θmax = network_data["branch"]["1"]["angmax"]
    θmin = network_data["branch"]["1"]["angmin"]
    W = Dict{Tuple{Int64, Int64}, Float64}()
    smax = Dict{Int64, Float64}()
    smin = Dict{Int64, Float64}()

    Demand = Dict{Int64, Dict{Int64, Float64}}()
    for t in 1:T 
        Demand[t] = Dict{Int64, Float64}()
    end
    w = Dict{Int64, Float64}()              ## priority level of load D
    cb = Dict{Int64, Float64}()             ## set of fire damage cost cᵢ at :b ∈ B
    cg = Dict{Int64, Float64}()
    cl = Dict{Tuple{Int64, Int64}, Float64}()

    for i in keys(network_data["bus"])
        b = network_data["bus"][i]["bus_i"]
        push!(B, b)
        Dᵢ[b]    = Vector{Int64}()
        Gᵢ[b]    = Vector{Int64}()
        out_L[b] = Vector{Int64}()
        in_L[b]  = Vector{Int64}()
        cb[b] = 1500. ## 0                ############# need to revise
    end

    for i in keys(network_data["load"])
        d = network_data["load"][i]["index"]
        b = network_data["load"][i]["load_bus"]
        w[d] = network_data["load"][i]["pd"] * 1e5                     ## priority level of load d

        push!(Dᵢ[b], d)
        push!(D, d)
        for t in 1:T 
            demand = network_data["load"][i]["pd"] * (1 + .2 * t)
            Demand[t][d] = demand
        end
    end

    for i in keys(network_data["gen"])
        g = network_data["gen"][i]["index"]
        b = network_data["gen"][i]["gen_bus"]

        push!(G, g)
        push!(Gᵢ[b], g)

        smax[g] = network_data["gen"][i]["pmax"]
        smin[g] = network_data["gen"][i]["pmin"]
        cg[g] = 1200.                               ############# need to revis
    end


    for i in keys(network_data["branch"])
        l = (network_data["branch"][i]["f_bus"], network_data["branch"][i]["t_bus"])

        if l ∉ L 
            push!(L, l) 
            push!(out_L[l[1]], l[2])
            push!(in_L[l[2]], l[1])

            _b[l] = network_data["branch"][i]["b_fr"]   ## total line charging susceptance
            W[l] = network_data["branch"][i]["rate_a"]              
            cl[l] = 500.                                ############# need to revise
        end
    end
    

    paramOPF = ParamOPF(_b, θmax, θmin, W, smax, smin)
    indexSets = IndexSets(D, G, unique(L), B ,T, [1:Ω...], Dᵢ, Gᵢ, out_L, in_L)
    paramDemand = ParamDemand(Demand, w, cb, cg, cl, 1e5)
 
     return (indexSets = indexSets, 
             paramOPF = paramOPF, 
             paramDemand = paramDemand)
end

function prepareSimulation(businfo::DataFrame, branchInfo::DataFrame, WFPI_Info::DataFrame;     
                                                            n::Int64 = 100  ## increase boundary
                                                            )
    
    busLocation = businfo[!, [1, 14, 15]]

    utm_ca = UTMfromLLA(11, true, wgs84)
    bus1_utm = Dict{Int64, Any}()
    row_num = 1
    for id in busLocation[:, 1]
        bus_lla = LLA(busLocation[row_num, 2], busLocation[row_num, 3])
        bus1_utm[id] = utm_ca(bus_lla)
        row_num += 1
    end



    xmin = floor(minimum(bus1_utm[i].x for i in keys(bus1_utm)))
    xmax = ceil(maximum(bus1_utm[i].x for i in keys(bus1_utm)))

    ymin = floor(minimum(bus1_utm[i].y for i in keys(bus1_utm)))
    ymax = ceil(maximum(bus1_utm[i].y for i in keys(bus1_utm)))

    
    bus_id_location = Dict{Int64, Tuple{Int64, Int64}}()
    bus_location_id = Dict{Tuple{Int64, Int64}, Int64}()
    row_num = 1
    for id in keys(bus1_utm)
        bus_lla = bus1_utm[id]
        bus_id_location[id] = Int.(ceil.((bus_lla.x - xmin, bus_lla.y - ymin) ./ 1000)) .+ n
        bus_location_id[bus_id_location[id]] = id
        row_num += 1
    end



    line_id_location = Dict{Int64, Any}()
    line_location_id = Dict{Any, NamedTuple{(:id, :WFPI), Tuple{Int64, Float64}}}()
    line_id_bus = Dict{Int64, Any}()
    row_num = 1
    for id in 1: nrow(branchInfo)
        WFPI_value = WFPI_Info[id, 3]/5
        from_bus = branchInfo[id, 2]
        to_bus = branchInfo[id, 3]
        line_id_bus[id] = (from_bus, to_bus)
        l1 = bus_id_location[from_bus]
        l2 = bus_id_location[to_bus]
        l = l2 .- l1
        l_25 = Int.(ceil.(l1 .+ .25 .* l))
        l_5 = Int.(ceil.(l1 .+ .5 .* l))
        l_75 = Int.(ceil.(l1 .+ .75 .*l))
        line_id_location[id] = [l1, l_25, l_5, l_75, l2]
        line_location_id[l1] = (id = id, WFPI = WFPI_value)
        line_location_id[l_25] = (id = id, WFPI = WFPI_value)
        line_location_id[l_5] = (id = id, WFPI = WFPI_value)
        line_location_id[l_75] = (id = id, WFPI = WFPI_value)
        line_location_id[l2] = (id = id, WFPI = WFPI_value)
        row_num += 1
    end



    xdiff = xmax - xmin
    x_grid_num = maximum(bus_id_location[i][1] for i in keys(bus_id_location)) + n

    ydiff = ymax - ymin
    y_grid_num = maximum(bus_id_location[i][2] for i in keys(bus_id_location)) + n

    return (x_grid_num = x_grid_num, 
            y_grid_num = y_grid_num, 
            line_location_id = line_location_id, 
            line_id_location = line_id_location, 
            line_id_bus = line_id_bus,
            bus_id_location = bus_id_location, 
            bus_location_id = bus_location_id)
end


function prepareScenarios( ;period_span::Int64 = 1, 
                            T::Int64 = 48, 
                            Ω::Int64 = 5, 
                            indexSets::IndexSets = indexSets, prob_fault::Float64 = 0.1, line_id_bus::Dict = line_id_bus)

    function check_vector_exist(key::Union{Tuple{Int64, Int64}, Int64}, dict::Dict)
        if key ∈ keys(dict)
            push!(dict[key], key)
        else
            dict[key] = []
            push!(dict[key], key)
        end
    end


    function check_vector_exist(key::Union{Tuple{Int64, Int64}, Int64}, dict::Dict, G::Dict)
        if length(key) == 1
            if key ∈ keys(dict)
                push!(dict[key], G[key]...)
            else
                dict[key] = []
                push!(dict[key], G[key]...)
            end
        else
            if key ∈ keys(dict)
                push!(dict[key], G[key[1]]...)
                push!(dict[key], G[key[2]]...)
            else
                dict[key] = []
                push!(dict[key], G[key[1]]...)
                push!(dict[key], G[key[2]]...)
            end
        end
    end 


    environmentInfo = Dict{Tuple{Int64, Int64}, CellEnvironmentInfo}()

    for i in 1:x_grid_num
        for j in 1:y_grid_num
            Pveg = wsample([-1.0, -.4, .4, .4], [0.1, .3, .5, 0.1], 1)[1]
            Pden = Pveg == -1.0 ? -1.0 : wsample([-1.0, -.3, 0.0, .3], [0.15, .25, .45, 0.15], 1)[1]
            E = wsample([10., 20., 30., 40.], [0.1, .3, .5, 0.1], 1)[1]
            # state = Pveg == -1.0 ? 0 : wsample([1, 2], [0.99995, 5e-5], 1)[1]
            state = Pveg == -1.0 ? 0 : 1
            environmentInfo[(i, j)] = CellEnvironmentInfo(Pveg, Pden, E, state)
        end
    end

    Ω_rv = Dict{Int64, RandomVariables}()
    Gᵢ = indexSets.Gᵢ
    B = indexSets.B
    G = indexSets.G
    L = indexSets.L
    for ω in 1:Ω 
        τ = rand(2:T)

        ub = Dict{Int64, Int64}()
        ug = Dict{Int64, Int64}()
        ul = Dict{Tuple{Int64, Int64}, Int64}() 

        vb = Dict{Int64, Int64}()
        vg = Dict{Int64, Int64}()
        vl = Dict{Tuple{Int64, Int64}, Int64}()

        Ibb = Dict{Int64, Vector{Int64}}()
        Ibg = Dict{Int64, Vector{Int64}}()
        Ibl = Dict{Int64, Vector{Tuple{Int64, Int64}}}()

        Igb = Dict{Int64, Vector{Int64}}()
        Igg = Dict{Int64, Vector{Int64}}()
        Igl = Dict{Int64, Vector{Tuple{Int64, Int64}}}()

        Ilb = Dict{Tuple{Int64, Int64}, Vector{Int64}}()
        Ilg = Dict{Tuple{Int64, Int64}, Vector{Int64}}()
        Ill = Dict{Tuple{Int64, Int64}, Vector{Tuple{Int64, Int64}}}()

        for b in B 
            ub[b] = 0
            vb[b] = 0
            Ibb[b] = []
            Ibg[b] = []
            Ibl[b] = []
        end

        for g in G 
            ug[g] = 0
            vg[g] = 0
            Igb[g] = []
            Igg[g] = []
            Igl[g] = []
        end

        for l in L 
            ul[l] = 0
            vl[l] = 0
            Ilb[l] = []
            Ilg[l] = []
            Ill[l] = []
        end


        forest = initialize(bus_id_location, line_location_id;  griddims = (x_grid_num, y_grid_num), 
                        environmentInfo = environmentInfo, time_span = 1)

        ## generate wildfire random variables
        for i in 1:floor((T/period_span))
            Agents.step!(forest, agent_step!, wildfire_ignition_step!, period_span)

            if sum(forest.line_fired) > 0 && sum(forest.bus_fired) > 0
                τ = i * period_span

                if sum(forest.line_fired) > 0
                    for I in findall(isequal(1), forest.line_fired)
                        id = line_location_id[I.I].id
                        vl[line_id_bus[id]] = 1
                        check_vector_exist(line_id_bus[id], Ill)
                        # check_vector_exist(line_id_bus[id], Ilb)
                        check_vector_exist(line_id_bus[id], Ilg, Gᵢ)

                        # push!(Ill[line_id_bus[id]], line_id_bus[id])
                        if line_id_bus[id] ∈ keys(Ilb)
                            push!(Ilb[line_id_bus[id]], line_id_bus[id]...)
                        else 
                            Ilb[line_id_bus[id]] = []
                            push!(Ilb[line_id_bus[id]], line_id_bus[id]...)
                        end
                        # push!(Ilg[line_id_bus[id]], Gᵢ[line_id_bus[id][1]]...)
                        # push!(Ilg[line_id_bus[id]], Gᵢ[line_id_bus[id][2]]...)
                    end
                end
                
                
                if sum(forest.bus_fired) > 0
                    for I in findall(isequal(1), forest.bus_fired)
                        bus_id = bus_location_id[I.I]
                        vb[bus_id] = 1
                        for id in Gᵢ[bus_id] 
                            vg[id] = 1 
                            check_vector_exist(id, Igg)
                            check_vector_exist(bus_id, Igb)
                        end
                        check_vector_exist(bus_id, Ibb)
                        check_vector_exist(bus_id, Ibg, Gᵢ)
                    end
                end 

                if sum(forest.line_fault) > 0
                    for I in findall(isequal(1), forest.line_fault)
                        id = line_location_id[I.I].id
                        ul[line_id_bus[id]] = 1
                    end
                end

                if sum(forest.bus_fault) > 0
                    for I in findall(isequal(1), forest.bus_fault)
                        bus_id = bus_location_id[I.I]
                        ub[bus_id] = 1
                    end
                end



                

                break
            end
        end

        ## generate fault random variables
        for b in B 
            fault = rand(Binomial(1,prob_fault), 1)[1]
            if fault == 1
                ub[b]  = fault
                for g in Gᵢ[b]
                    ug[g] = rand(Binomial(1,prob_fault), 1)[1]
                end
            end
        end

        # for l in L 
        #     fault = rand(Binomial(1,prob_fault), 1)[1]
        #     if fault == 1
        #         ub[b]  = 
        #         for g in Gᵢ[b]
        #             ug[g] = rand(Binomial(1,prob_fault), 1)[1]
        #         end
        #     end
        # end


        Ω_rv[ω] = RandomVariables(τ, ub, ug, ul, vb, vg, vl, Ibb, Ibg, Ibl, Igb, Igg, Igl, Ilb, Ilg, Ill)
    end

    return Ω_rv
end

