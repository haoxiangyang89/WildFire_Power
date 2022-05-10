businfo = CSV.read("/Users/aaron/RTS-GMLC/RTS_Data/SourceData/bus.csv", DataFrame)
branchInfo = CSV.read("/Users/aaron/RTS-GMLC/RTS_Data/SourceData/branch.csv", DataFrame)

function prepareSimulation(businfo::DataFrame, branchInfo::DataFrame;     
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
    line_location_id = Dict{Any, Int64}()
    line_id_bus = Dict{Int64, Any}()
    row_num = 1
    for id in 1: nrow(branchInfo)
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
        line_location_id[l1] = id
        line_location_id[l_25] = id
        line_location_id[l_5] = id
        line_location_id[l_75] = id
        line_location_id[l2] = id
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
            bus_id_location = bus_id_location, 
            bus_location_id = bus_location_id)
end

(x_grid_num, y_grid_num, line_location_id, 
                            line_id_location, 
                            bus_id_location, 
                            bus_location_id) = prepareSimulation(businfo, branchInfo)

 




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


forest = initialize(bus_id_location, line_location_id;  griddims = (x_grid_num, y_grid_num), 
                        environmentInfo = environmentInfo)
wildfire_ignition_step!(forest)



step!(forest, agent_step!, wildfire_ignition_step!, 1)
sum(forest.ignition)
sum(forest.line_fired)
sum(forest.bus_fired)


Random.seed!(2)
forest = initialize(bus_id_location; griddims = (x_grid_num, y_grid_num), environmentInfo = environmentInfo)
burnt_percentage(f) = count(t == 2 for t in f.ignition) / prod(size(f.ignition))
mdata = [burnt_percentage]

_, data = run!(forest, agent_step!, wildfire_ignition_step!, 48; mdata)
data


for I in findall(isequal(1), forest.wildfire_exist)
    bus_location_id[I.I]
end


Random.seed!(2)
forest = initialize(griddims = (n, n))
burnt_percentage(f) = count(t == 2 for t in f.ignition) / prod(size(f.ignition))
mdata = [burnt_percentage]

_, data = run!(forest, agent_step!, wildfire_ignition_step!, 10; mdata)
data


forest = initialize(griddims = (n, n))
Agents.step!(forest, agent_step!, wildfire_ignition_step!, 1)

plotkwargs = (
    add_colorbar = false,
    heatarray = :ignition,
    heatkwargs = (
        colorrange = (0, 3),
        colormap = cgrad([:white, :green, :red, :darkred]; categorical = true),
    ),
)
fig, _ = abmplot(forest; plotkwargs...)
fig


Random.seed!(10)
forest = initialize( ;  griddims = (n, n), 
environmentInfo = environmentInfo)
add_agent!(forest) # Add one dummy agent so that abm_video will allow us to plot.
abmvideo(
    "forest.mp4",
    forest,
    agent_step!,
    wildfire_ignition_step!;
    as = 0,
    framerate = 5,
    frames = 20,
    spf = 5,
    title = "Forest Fire",
    plotkwargs...,
)





