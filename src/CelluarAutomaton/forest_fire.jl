using Agents, Random
using InteractiveDynamics
using CairoMakie

@agent Automata GridAgent{2} begin end

function forest_fire(; density = 0.7, griddims = (100, 100))
    space = GridSpace(griddims; periodic = false, metric = :euclidean)
    # The `trees` field is coded such that
    # Empty = 0, Green = 1, Burning = 2, Burnt = 3
    forest = ABM(Automata, space; properties = (trees = zeros(Int, griddims),))
    for I in CartesianIndices(forest.trees)
        if rand(forest.rng) < density
            # Set the trees at the left edge on fire
            forest.trees[I] = I[1] == 1 ? 2 : 1
        end
    end
    return forest
end

forest = forest_fire()

function tree_step!(forest)
    # Find trees that are burning (coded as 2)
    for I in findall(isequal(2), forest.trees)
        for idx in nearby_positions(I.I, forest)
            # If a neighbor is Green (1), set it on fire (2)
            if forest.trees[idx...] == 1
                forest.trees[idx...] = 2
            end
        end
        # Finally, any burning tree is burnt out (2)
        forest.trees[I] = 3
    end
end

Agents.step!(forest, dummystep, tree_step!, 1)
count(t == 3 for t in forest.trees) # Number of burnt trees on step 1

Agents.step!(forest, dummystep, tree_step!, 10)
count(t == 3 for t in forest.trees) # Number of burnt trees on step 11

Random.seed!(2)
forest = forest_fire(griddims = (20, 20))
burnt_percentage(f) = count(t == 3 for t in f.trees) / prod(size(f.trees))
mdata = [burnt_percentage]

_, data = run!(forest, dummystep, tree_step!, 10; mdata)
data


forest = forest_fire()
Agents.step!(forest, dummystep, tree_step!, 1)

plotkwargs = (
    add_colorbar = false,
    heatarray = :trees,
    heatkwargs = (
        colorrange = (0, 3),
        colormap = cgrad([:white, :green, :red, :darkred]; categorical = true),
    ),
)
fig, _ = abmplot(forest; plotkwargs...)
fig


Random.seed!(10)
forest = forest_fire(density = 0.6)
add_agent!(forest) # Add one dummy agent so that abm_video will allow us to plot.
abmvideo(
    "forest.mp4",
    forest,
    dummystep,
    tree_step!;
    as = 0,
    framerate = 5,
    frames = 20,
    spf = 5,
    title = "Forest Fire",
    plotkwargs...,
)