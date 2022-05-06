using Distributions, DelimitedFiles, Random, StatsBase, Plots


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
test = simulate(simulationModel,n_iterations)
plot(heatmap(test.pa,c = :viridis))

