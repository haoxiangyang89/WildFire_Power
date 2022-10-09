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

## function to compute the probability of outage
function Probability_fault(lightningDensity::Float64, maxWind::Float64, L::Union{Float64, Int64})
    lightningDensity = 0.9
    maxWind = 8
    W = [1, lightningDensity, lightningDensity * maxWind, lightningDensity^2, maxWind^2]
    B = [-13.2719 2.8091 -.2515 -.1438 .0963; 
    -13.9225 3.1331 -.0791 -.2645 .0727;
    -11.9084 1.6476 -.0387 -.1069 .0595;
    -11.8423 1.0980 .2043 -.2692 .0515;
    -13.5802 3.2945 -.2511 -.2008 .1078; ]

    L = 73
    prob = 1 .- exp.(- L * exp.(B * W))

    return prob
end



function initialize(relative_location::Dict{Int64, Tuple{Int64, Int64}}, line_location_id::Dict{Any, NamedTuple{(:id, :WFPI, :Length), Tuple{Int64, Float64, Union{Float64, Int64}}}} ;  griddims::Tuple{Int64, Int64} = (100, 100), 
                        environmentInfo::Dict{Tuple{Int64, Int64}, CellEnvironmentInfo} = environmentInfo, time_span::Int64 = 1)

    space = GridSpace(griddims, periodic = false)
    # The `trees` field is coded such that
    # No ignition = 0, ignition = 1, burnt = 2
    forest = ABM(
        CellInfo, space;
        properties = (
                        ignition = zeros(Int, griddims),  ## 0, 1, 2
                        fault_WFPI = zeros(Float64, griddims),  
                        busExist = zeros(Int, griddims),  ## binary
                        lineExist = zeros(Int, griddims),  ## binary
                        lineLength = zeros(Union{Float64, Int64}, griddims),  ## binary * length
                        busFired = zeros(Int, griddims),     ## ## binary
                        lineFired = zeros(Int, griddims),   ## binary, by natural wildfire
                        busFault = zeros(Int, griddims),    ## binary
                        lineFault = zeros(Int, griddims),   ## binary, fault caused by weather, wildlife, human, lightning
                        windInfo = [1.0, 0.0],               ## [windSpeed, windDirection, lightningDensity]
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
        forest.busExist[pos...] = 1
    end

    for pos in keys(line_location_id)
        forest.lineExist[pos...] = 1
        forest.fault_WFPI[pos...] = line_location_id[pos].WFPI/100    
        forest.lineLength[pos...] = line_location_id[pos].Length/3   
    end
    return forest
end



function agent_step!(agent, forest)
    windInfo = WindInfo(forest.windInfo[1], forest.windInfo[2])
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
                if neighbor.state ≥ 2
                    p = Probability_burn(neighbor, agent; windInfo = windInfo, P₀ = .58)
                    agent.state = rand(1)[1] <= p ? 2 : 1
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


        for I in findall(isequal(1), forest.lineExist)
            if forest.ignition[I] == 1  ## if the cell with one bus is burning, we say this bus is fired
                forest.lineFault[I] = rand(forest.rng) <= forest.fault_WFPI[I] ? 1 : 0

                dataClassification = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
                windDays = [5, 25, 200, 250, 150, 75, 30, 10, 3, 1]
                lightningDensity = wsample(dataClassification .* 0.09, [0.75, 0.2, 0.15, 0.15, 0.1, 0.08, 0.05, 0.03, 0.02, 0.01], 1)[1]
                maxWind = wsample(dataClassification .* 1.4, windDays./sum(windDays), 1)[1]

                forest.lineFault[I] = rand(forest.rng) <= minimum(Probability_fault(lightningDensity, maxWind, forest.lineLength[I])) ? 1 : 0
            end

            if forest.ignition[I] == 2  ## if the cell with one branch burnt, we say this line is fired
                forest.lineFired[I] = 1
            end
        end


        for I in findall(isequal(1), forest.busExist)
            if forest.lineFault[I] == 1  ## if the cell with one bus is burning, we say this bus is fired
                forest.busFault[I] = 1
            end

            if forest.ignition[I] == 2  ## if the cell with one bus is burnt, we say this bus is fired
                forest.busFired[I] = 1
            end
        end

        dataClassification = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        windDays = [5, 25, 200, 250, 150, 75, 30, 10, 3, 1]
        forest.windInfo[1] = wsample(dataClassification * 1.1, windDays./sum(windDays), 1)[1]
        forest.windInfo[2] = wsample([0, 80, 150, 200, 240, 280, 310, 350], [0.1, .12, .11, 0.13, 0.14, 0.09, 0.16, 0.15], 1)[1]
end