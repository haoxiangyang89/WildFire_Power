## data struct for cellular automaton

struct CellularAutomaton
    pa                   ::Array          # 2D array that can be 0 (un-occupied) or 1 (occupied)
    caIndex              ::Array          # Index to reference cells in the landscape
    suitability          ::Array          # 2D array containing scaled suitability values
    dispersalProbability ::Float64        # Probability that a cell will disperse at each time step
    meanDispersal        ::Float64        # Mean dispersal distance
end


## data struct for wildfire

struct Wildfire
    


end

