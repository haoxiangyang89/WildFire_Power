windInfo = WindInfo(.2, 176)


environmentInfo = Dict{Tuple{Int64, Int64}, CellEnvironmentInfo}()

for i in 1:n
    for j in 1:n
        environmentInfo[(i, j)] = CellEnvironmentInfo(.1, 1.1, 100.0, 1)
    end
end





forest = initialize( ;  griddims = (n, n), 
                        environmentInfo = environmentInfo)
wildfire_step!(forest; windInfo = windInfo)

Agents.step!(forest, agent_step!, wildfire_step!, 1)
step!(forest, agent_step!, wildfire_step!, 1)






Random.seed!(2)
forest = initialize(griddims = (n, n))
burnt_percentage(f) = count(t == 2 for t in f.wildfire) / prod(size(f.wildfire))
mdata = [burnt_percentage]

_, data = run!(forest, agent_step!, wildfire_step!, 10; mdata)
data


forest = initialize(griddims = (n, n))
Agents.step!(forest, agent_step!, wildfire_step!, 1)

plotkwargs = (
    add_colorbar = false,
    heatarray = :wildfire,
    heatkwargs = (
        colorrange = (0, 3),
        colormap = cgrad([:white, :green, :red, :darkred]; categorical = true),
    ),
)
fig, _ = abmplot(forest; plotkwargs...)
fig


Random.seed!(10)
forest = forest_fire(griddims = (n, n))
add_agent!(forest) # Add one dummy agent so that abm_video will allow us to plot.
abmvideo(
    "forest.mp4",
    forest,
    agent_step!,
    wildfire_step!;
    as = 0,
    framerate = 5,
    frames = 20,
    spf = 5,
    title = "Forest Fire",
    plotkwargs...,
)

