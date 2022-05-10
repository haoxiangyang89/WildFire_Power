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







