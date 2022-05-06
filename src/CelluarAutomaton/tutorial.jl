# https://grimmel.github.io/posts/2020/10/blog-post-1/
using Distributions, DelimitedFiles, Random, StatsBase, Plots


struct CellularAutomaton
    pa                   ::Array          # 2D array that can be 0 (un-occupied) or 1 (occupied)
    caIndex              ::Array          # Index to reference cells in the landscape
    suitability          ::Array          # 2D array containing scaled suitability values
    dispersalProbability ::Float64        # Probability that a cell will disperse at each time step
    meanDispersal        ::Float64        # Mean dispersal distance
end


function selectProportion(pa, caIndex, dispersalProbability)
    nPresences = Int(sum(pa))
    total = Array{CartesianIndex}(undef, nPresences)
    counter = 1
    for idx in caIndex
        if pa[idx] === 1.0
            total[counter] = idx
            counter+=1
        end
    end
    numberDispersing = sample(total,rand(Binomial(nPresences,dispersalProbability)))
    return numberDispersing
end


function newPos(meanDispersal, cartIdx)
    distance = rand(Exponential(meanDispersal),1)
    # + 0.75 ensures dispersal outside of the initial cell
    distance = distance[1] + 0.75
    angle = 360.0*rand()
    # Remember 1 = y, 2 = x
    x = Int(round(cos(deg2rad(angle))*distance,digits=0)) + cartIdx[2]
    y = Int(round(sin(deg2rad(angle))*distance,digits=0)) + cartIdx[1]
    return (y,x)
end


function colonise(cellularAutomaton::CellularAutomaton)
    shape = size(cellularAutomaton.pa)
    dCells = selectProportion(cellularAutomaton.pa, cellularAutomaton.caIndex, cellularAutomaton.dispersalProbability)
    for cartIdx in dCells
        newXY = newPos(cellularAutomaton.meanDispersal,cartIdx)
        if newXY[2]>=1 && newXY[2] <= shape[2] && newXY[1] >=1 && newXY[1]<=shape[1]
            cellularAutomaton.pa[newXY[1],newXY[2]] = 1
        end
    end
    return cellularAutomaton
end


function test_exponential(dispersal)
    graph = zeros(50,50)
    initx = 25
    inity = 20
    for i in 1:1000000
        pos = newPos(dispersal, CartesianIndex(initx, inity))
        if pos[1]>=1 && pos[1]<=50 && pos[2]>=1 && pos[2]<=50
            graph[pos[2],pos[1]] +=1
        end
    end
    return(graph)
end
test1 = test_exponential(1)
test3 = test_exponential(3)
test6 = test_exponential(6)
test12 = test_exponential(12)
plot(heatmap(test1,c = :viridis), heatmap(test3,c = :viridis), heatmap(test6,c = :viridis), heatmap(test12,c = :viridis))






function extinction(cellularAutomaton::CellularAutomaton)
    for idx in cellularAutomaton.caIndex
        if cellularAutomaton.pa[idx] === 1.0
            survived = rand(Bernoulli(cellularAutomaton.suitability[idx]),1)
            if survived[1] === false
                cellularAutomaton.pa[idx] = 0.0
            end
        end
    end
    return cellularAutomaton
end

function simulate(cellularAutomaton::CellularAutomaton, iterations)
    for i in 1:iterations
        cellularAutomaton = colonise(cellularAutomaton)
        cellularAutomaton = extinction(cellularAutomaton)
    end
    return cellularAutomaton
end


landscape = [1. 0 0 1; 1 1 0 0; 0 0 0 1]
ls_dimension = size(landscape)
suit = landscape./100
pa = ones(ls_dimension)
caIndex = CartesianIndices(pa)
dispersalProbability = 0.8
distance = 1.0
simulationModel = CellularAutomaton(pa, caIndex, suit, dispersalProbability, distance)
n_iterations = 100
@time simulate(simulationModel,n_iterations)
# open("D:/test/outJULIA.asc","w") do io
#     write(io,"NCOLS 400\nNROWS 400\nXLLCORNER 0\nYLLCORNER 0\nCELLSIZE 100\nNODATA_value -9999\n")
#     writedlm(io,simulationModel.pa)
# end



test = simulate(simulationModel,n_iterations)
plot(heatmap(test.pa,c = :viridis))
