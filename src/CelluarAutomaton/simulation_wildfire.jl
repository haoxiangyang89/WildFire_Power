using Agents, Random
using InteractiveDynamics
using CairoMakie

function initialize( ;  griddims::Tuple{Int64, Int64} = (100, 100), 
                        environmentInfo::Dict{Tuple{Int64, Int64}, CellEnvironmentInfo} = environmentInfo)

    space = GridSpace(griddims, periodic = false)
    # The `trees` field is coded such that
    # Empty = 0, No wildfire = 1, ingition = 2
    forest = ABM(
        CellInfo, space;
        properties = (wildfire = zeros(Int, griddims),
        scheduler = Schedulers.randomly)
    )

    # populate the model with agents, adding equal amount of the two types of agents
    # at random positions in the model


    id = 1
    for I in CartesianIndices(forest.wildfire)
        ## need to find a logic for the following fields
        Pveg = environmentInfo[I.I].Pveg         ## Vegetation density (limited = -0.4, medium = 0., high = 0.3)
        Pden = environmentInfo[I.I].Pden         ## Kind of vegetation
        E = environmentInfo[I.I].E               ## elevation of the cell
        state = environmentInfo[I.I].initial_state 
        agent = CellInfo(id, (1, 1), Pveg, Pden, E, state)
        add_agent!(agent, I.I, forest)
        id = id + 1

        forest.wildfire[I] = min(state, 2) ## Empty = 0, No wildfire = 1, ingition = 2
    end
    return forest
end






function agent_step!(agent, forest; windInfo::WindInfo = windInfo)
    if agent.state >= 2
        ## If it is burning previously, then it is burnt now
        agent.state = 3
    elseif agent.state == 1

        ## If this cell has a wildfire, then its state becomes burning
        if forest.wildfire[agent.pos...] == 2
            agent.state = 2
        else
            ## If nearby agent is burnt, then its state becomes burning
            for neighbor in nearby_agents(agent, forest)
                if neighbor.state == 3
                    p = Probability_burn(neighbor, agent; windInfo = windInfo, P₀ = .9)
                    agent.state = rand(forest.rng) <= p ? 2 : agent.state
                end
            end
        end
    end

    forest.wildfire[agent.pos...] = min(agent.state, 2)

    return 
end


function wildfire_step!(forest::AgentBasedModel; windInfo::WindInfo = windInfo)
    # Find trees that are burning (coded as 2)
    for id in 1:nagents(forest)
        agent = forest[id]
        if agent.state ≥ 2
            forest.wildfire[agent.pos...] == 2
        elseif agent.state == 0
            forest.wildfire[agent.pos...] == 0
        elseif agent.state == 1
            ## at time t, this cell has vegetation, and no ingition
            prob = []
            for neighbor in nearby_agents(agent, forest)
                if neighbor.state == 0  ## 
                    p = 0
                elseif neighbor.state == 1
                    p = Probability_burn(neighbor, agent; windInfo = windInfo, P₀ = .3)
                elseif neighbor.state == 2
                    p = Probability_burn(neighbor, agent; windInfo = windInfo, P₀ = .58)
                elseif neighbor.state == 3
                    p = Probability_burn(neighbor, agent; windInfo = windInfo, P₀ = .8)
                end
                push!(prob, p)
            end
            forest.wildfire[agent.pos...] = rand(forest.rng) <= sum(prob)/length(prob) ? 2 : 1
        end  
    end
end












