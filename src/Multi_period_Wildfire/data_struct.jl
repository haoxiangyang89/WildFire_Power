#############################################################################################
################################### Indices & Index Sets ####################################
#############################################################################################
struct IndexSets
    D       :: Vector{Int64}                  ## set of load demand
    G       :: Vector{Int64}                  ## set of generators
    L       :: Vector{Tuple{Int64, Int64}}    ## set of transmission lines
    B       :: Vector{Int64}                  ## set of buses
    T       :: Int64  ## set of time periods  1:T
    Ω       :: Vector{Int64}
    _D
    _G
    out_L
    in_L
end
## _D, _G, out_L, in_L

struct ParamDemand
    demand        :: Dict{Int64, Dict{Int64, Float64}}       ## set of load demand  Dict{t, Dict{d, demand}}
    w             :: Dict{int64, Float64}                    ## priority level of load D
    cb            :: Dict{int64, Float64}                    ## set of fire damage cost cᵢ at i ∈ B
    cg            :: Dict{int64, Float64}                    ## set of fire damage cost cᵢ at i ∈ G
    cl            :: Vector{Tuple{Int64, Float64}}           ## set of fire damage cost cᵢ at i ∈ L
end


struct ParamOPF  ## included in a period dict
    b       :: Dict{Tuple{Int64, Int64}, Float64}            
    θmax    :: Float64
    θmin    :: Float64
    W       :: Dict{Tuple{Int64, Int64}, Float64}
    smax    :: Dict{Int64, Float64}
    smin    :: Dict{Int64, Float64}
end







#############################################################################################
####################################   Data Structure   #####################################
#############################################################################################
struct CutCoefficient
    v               ::Dict{Int64,Dict{Int64, Float64}} # [i][k] where i is the iteration index, and k is the scenario index
    πb               ::Dict{Int64,Dict{Int64, Vector{Float64}}}  # [[1.,2.,3.],[1.,2.,3.]]  -- push!(Π, π) to add a new element
    πg               ::Dict{Int64,Dict{Int64, Vector{Float64}}}
    πl               ::Dict{Int64,Dict{Int64, Vector{Float64}}}
end



## for period t with realization ω
struct RandomVariables  
    τ   ::Int64

    ub   ::Dict{Int64, Bool}                                                ## whether there exists a fault
    ug   ::Dict{Int64, Bool}   
    ul   ::Dict{Tuple{Int64, Float64}, Bool} 

    vb   ::Dict{Int64, Bool}                                                ## whether there exists a fire caused by natural condition
    vg   ::Dict{Int64, Bool}   
    vl   ::Dict{Tuple{Int64, Float64}, Bool} 

    Ibb   ::Dict{Int64, Int64}                                              ## the set of buses which is affected by a bus
    Ibg   ::Dict{Int64, Int64}                                              ## the set of generators which is affected by a bus
    Ibl   ::Dict{Int64, Tuple{Int64, Float64}}                              # the set of lines which is affected by a bus

    Igb   ::Dict{Int64, Int64}                                              ## the set of buses which is affected by a generators
    Igg   ::Dict{Int64, Int64}                                              ## the set of generators which is affected by a generators
    Igl   ::Dict{Int64, Tuple{Int64, Float64}}                              ## the set of lines which is affected by a generators

    Ilb   ::Dict{Tuple{Int64, Float64}, Int64}                              ## the set of buses which is affected by a line
    Ilg   ::Dict{Tuple{Int64, Float64}, Int64}                              ## the set of generators which is affected by a line
    Ill   ::Dict{Tuple{Int64, Float64}, Tuple{Int64, Float64}}              ## the set of lines which is affected by a line
end

# Ω_rv = Dict{Int64,RandomVariables}()
# Ω_rv[ω]= RandomVariables(...) # ω is a number
# Ω_rv[ω].τ




struct GurobiModelInfo
    model           :: Model
    x               :: Array{VariableRef} ## for current state, x is the number of generators
    y               :: Array{VariableRef} ## amount of electricity
    slack           :: Matrix{VariableRef}
    num_Ω           :: Int64
end


## data structure for levelset method
mutable struct FunctionInfo
    x_his        :: Dict{Int64, Dict{Symbol, Vector}}  ## record every x_j point
    G_max_his    :: Dict{Int64, Float64}          ## record max(g[k] for k in 1:m)(x_j)
    f_his        :: Dict{Int64, Float64}          ## record f(x_j)
    df           :: Dict{Symbol, Vector{Float64}}
    dG           :: Dict{Int64, Vector{Float64}}  ## actually is a matrix.  But we use dict to store it
    G            :: Dict{Int64, Float64}          
end




struct ModelInfo
    model :: Model
    x     :: Vector{VariableRef}
    y     :: VariableRef
    z     :: VariableRef
end







struct LevelSetMethodParam
    μ             ::Float64   ## param for adjust α
    λ             ::Float64   ## param for adjust level
    threshold     ::Float64   ## threshold for Δ
    nxt_bound     ::Float64   ## lower bound for solving next iteration point π
    max_iter      ::Int64     
    Output        ::Int64     ## Gurobi Output parameter
    Output_Gap    ::Bool      ## if True will print Δ info
    Adj           ::Bool      ## whether adjust oracle lower bound
end



################### nonanticipativity for multistage problem #######################

function recursion_scenario_tree(pathList::Vector{Int64}, P::Float64, scenario_sequence::Dict{Int64, Dict{Int64, Any}}, t::Int64;   
    Ω::Dict{Int64,Dict{Int64,RandomVariables}} = Ω, prob::Dict{Int64,Vector{Float64}} = prob, T::Int64 = 2)

    if t <= T
        for ω_key in keys(Ω[t])

            pathList_copy = copy(pathList)
            P_copy = copy(P)

            push!(pathList_copy, ω_key)
            P_copy = P_copy * prob[t][ω_key]

            recursion_scenario_tree(pathList_copy, P_copy, scenario_sequence, t+1, Ω = Ω, prob = prob, T = T)
        end
    else
        if haskey(scenario_sequence, 1)
            scenario_sequence[maximum(keys(scenario_sequence))+1] = Dict(1 => pathList, 2 => P)
        else
            scenario_sequence[1] = Dict(1 => pathList, 2 => P)
        end
        return scenario_sequence
    end

end

# scenario_sequence = Dict{Int64, Dict{Int64, Any}}()  ## the first index is for scenario index, the second one is for stage
# pathList = Vector{Int64}()
# push!(pathList, 1)
# P = 1.0

# recursion_scenario_tree(pathList, P, scenario_sequence, 2, T = T)
# scenario_tree = scenario_sequence



## sampling function 
function DrawSamples(scenario_sequence::Dict{Int64, Dict{Int64, Any}})
    # draw f, A, B, C, b from Ωₜ according to distribution P
    P = Vector{Float64}()
    for key in keys(scenario_sequence)
        push!(P, scenario_sequence[key][2])
    end
    items = [i for i in keys(scenario_sequence)]
    weights = Weights(P)
    j = sample(items, weights)
    return j
end


## form scenarios
function SampleScenarios(scenario_sequence::Dict{Int64, Dict{Int64, Any}}; T::Int64 = 5, M::Int64 = 30)
    ## a dict to store realization for each stage t in scenario k
    scenarios = Dict{Int64, Int64}()
    for k in 1:M
          scenarios[k] = DrawSamples(scenario_sequence)
    end
    return scenarios
end





function round!(a::Float64)               ## a = 1.3333e10
    b = floor(log10(a))                   ## b = 10
    c = round(a/10^b,digits = 2)          ## c = 1.33
    d = c * 10^b                          ## d = 1.33e10
    return [b, c, d]
end

