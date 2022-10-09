## based on Agent.jl v4.5
using Agents, Random
using InteractiveDynamics
using CairoMakie
using Distributions


include("CelluarAutomaton/data_struct.jl")
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



function initialize(relative_location::Dict{Int64, Tuple{Int64, Int64}}, line_location_id::Dict{Any, Int64} ;  griddims::Tuple{Int64, Int64} = (100, 100), 
                        environmentInfo::Dict{Tuple{Int64, Int64}, CellEnvironmentInfo} = environmentInfo, time_span::Int64 = 1)

    space = GridSpace(griddims, periodic = false)
    # The `trees` field is coded such that
    # No ignition = 0, ignition = 1, burnt = 2
    windInfo = WindInfo(.1, 176)
    forest = ABM(
        CellInfo, space;
        properties = (
                        ignition = zeros(Int, griddims),
                        bus_exist = zeros(Int, griddims),
                        line_exist = zeros(Int, griddims),
                        bus_fired = zeros(Int, griddims), 
                        line_fired = zeros(Int, griddims), 
                        windInfo = windInfo,
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
    end
    return forest
end






function agent_step!(agent, forest)
    windInfo = forest.windInfo
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
        windInfo = forest.windInfo
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
            if forest.ignition[I] >= 1  ## if the cell with one bus has burning or burnt fire, we say this bus is fired
                forest.bus_fired[I] = 1
            end
        end


        for I in findall(isequal(1), forest.line_exist)
            if forest.ignition[I] == 2  ## if the cell with one branch has burnt fire, we say this line is fired
                forest.line_fired[I] = 1
            end
        end
        wind_speed = wsample([.1, .2, .4, .7], [0.2, .4, .3, 0.1], 1)[1]
        wind_direction = wsample([10, 80, 120, 190, 250, 300, 340], [0.1, .2, .3, 0.1, 0.1, 0.1, .1], 1)[1]
        forest.windInfo = WindInfo(wind_speed, wind_direction)

end







