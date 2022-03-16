using JuMP, Gurobi, Random

const GRB_ENV = Gurobi.Env()


include("data_struct.jl")

## input data

# include("generationTest.jl")
# include("runtests_small2.jl")
# include("runtests_small3.jl")


########################################################################################################################################################
############################################  auxiliary function: nonanticipativity for multistage problem #############################################
########################################################################################################################################################
function recursion_scenario_constraint(pathList::Vector{Int64}, P::Float64, scenario_sequence::Dict{Int64, Dict{Int64, Any}}, t::Int64;   
                    Ω::Dict{Int64,Dict{Int64,RandomVariables}} = Ω, prob::Dict{Int64,Vector{Float64}} = prob, T::Int64 = 2, gurobiModelInfo::GurobiModelInfo = gurobiModelInfo)

    if t <= T
        for ω_key in keys(Ω[t])

            pathList_copy = copy(pathList)
            P_copy = copy(P)

            push!(pathList_copy, ω_key)
            P_copy = P_copy * prob[t][ω_key]

            ## nonanticipativity for multi-stage problem
            if t < T
                if haskey(scenario_sequence, 1)
                    first =  maximum(keys(scenario_sequence)) + 1
                    last  =  maximum(keys(scenario_sequence)) + gurobiModelInfo.num_Ω^(T-t)
                    @constraint(gurobiModelInfo.model, [i = first, j in (first + 1): last], gurobiModelInfo.x[:, t, i] .== gurobiModelInfo.x[:, t, j]) 
                    @constraint(gurobiModelInfo.model, [i = first, j in (first + 1): last], gurobiModelInfo.y[:, t, i] .== gurobiModelInfo.y[:, t, j]) 
                    @constraint(gurobiModelInfo.model, [i = first, j in (first + 1): last], gurobiModelInfo.slack[t, i] .== gurobiModelInfo.slack[t, j]) 

                else
                    @constraint(gurobiModelInfo.model, [i = 1, j in 2:gurobiModelInfo.num_Ω^(T-t)], gurobiModelInfo.x[:, t, i] .== gurobiModelInfo.x[:, t, j]) 
                    @constraint(gurobiModelInfo.model, [i = 1, j in 2:gurobiModelInfo.num_Ω^(T-t)], gurobiModelInfo.y[:, t, i] .== gurobiModelInfo.y[:, t, j]) 
                    @constraint(gurobiModelInfo.model, [i = 1, j in 2:gurobiModelInfo.num_Ω^(T-t)], gurobiModelInfo.slack[t, i] .== gurobiModelInfo.slack[t, j]) 
                end
            end
            recursion_scenario_constraint(pathList_copy, P_copy, scenario_sequence, t+1, Ω = Ω, prob = prob, T = T, gurobiModelInfo = gurobiModelInfo)
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







################################################################################################################################################
############################################################     Gurobi function   #############################################################
################################################################################################################################################

function gurobiOptimize!(Ω::Dict{Int64,Dict{Int64,RandomVariables}}, prob::Dict{Int64,Vector{Float64}}, StageCoefficient::Dict{Int64, StageData}; 
                        binaryInfo::BinaryInfo = binaryInfo)

    (A, n, d) = (binaryInfo.A, binaryInfo.n, binaryInfo.d)
    T = length(Ω); num_Ω = length(Ω[1]);
    W = num_Ω^(T-1) # number of scenarios


    model = Model( optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV), 
                                                "OutputFlag" => 1, 
                                                "Threads" => 1) 
                                                )
    @variable(model, x[i = 1:d, t = 1:T, ω in 1:W] >= 0, Int)   ## for current state, x is the number of generators will be built in this stage
    @variable(model, y[i = 1:d, t = 1:T, ω in 1:W] >= 0)        ## amount of electricity
    @variable(model, slack[t = 1:T, ω in 1:W] >=0 )
    @variable(model, S[i = 1:d, t = 1:T, ω in 1:W] >=0 )

    @constraint(model, [t in 1:T, ω in 1:W], S[:, t, ω] .== sum(x[:, j, ω] for j in 1:t ) )
    @constraint(model, [t in 1:T], S[:, t, :] .<= StageCoefficient[t].ū)
    @constraint(model, [t in 1:T, ω in 1:W], y[:,t, ω] .<= StageCoefficient[t].h * StageCoefficient[t].N * (S[:, t, ω] + StageCoefficient[t].s₀))

    ## nonanticipativity for multistage problem 
    @constraint(model, [ω in 2:W], x[:, 1, 1] .== x[:, 1, ω])  ## nonanticipativity for 2-stage problem
    gurobiModelInfo = GurobiModelInfo(model, x, y, slack, num_Ω)
    scenario_sequence = Dict{Int64, Dict{Int64, Any}}()  ## the first index is for scenario index, the second one is for stage
    pathList = Vector{Int64}()
    push!(pathList, 1)

    recursion_scenario_constraint(pathList, 1.0, scenario_sequence, 2, Ω = Ω, prob = prob, T = T, gurobiModelInfo = gurobiModelInfo)
    
    scenario_tree = scenario_sequence
    #############################################################################################################################
    @constraint(model, [t in 1:T, ω in 1:W], sum(y[i, t, ω] for i in 1:d) + slack[t, ω] >= Ω[t][scenario_tree[ω][1][t]].d[1] )

    @objective(model, Min, sum( sum( scenario_tree[ω][2] * (StageCoefficient[t].c1' * x[:, t, ω] + StageCoefficient[t].c2' * y[:, t, ω] + StageCoefficient[t].penalty * slack[t, ω]) for t in 1:T ) for ω in 1:W) )
    optimize!(model)
    ####################################################### solve the model and display the result ###########################################################

    # JuMP.objective_value(model)
    # JuMP.value.(x[:, 1, :])
    # JuMP.value.(x[:, 2, :])
    # JuMP.value.(x[:, 3, :])

    # JuMP.value.(y[:, 2, :])

    # JuMP.value.(slack)
    gurobiResult = Dict(1=> JuMP.objective_value(model), 2 =>  JuMP.value.(x[:, 1, 1]), 3=>  JuMP.value.(y[:, 1, 1]))

    return gurobiResult
end



# gurobiOptimize!(Ω, prob, StageCoefficient,
#                         binaryInfo = binaryInfo)