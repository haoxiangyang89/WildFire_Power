"""
    Enhanded_LC converges to optimal solution with sub-optimal value within 8 iterations
    However, LC cannot arrive to optimal solution within 200 iterations.

"""

##########################################################################################
############################  To generate Stage Data  ####################################
##########################################################################################
ū = [3., 5., 3.] # ū = [4.,10.,10.,1.,45.,4.]

binaryInfo = binarize_gen(ū)
(A, n, d) = (binaryInfo.A, binaryInfo.n, binaryInfo.d)


c1 = [[15000., 13000, 10000], [10000., 6000, 4000]]
c2 = [[1.9, 4.3, 5.0], [2.5, 3.1, 4.0]]

StageCoefficient = Dict{Int64,StageData}()

s₀ = [1,1,1]
penalty = 1e3
N = Array{Float64,2}(undef, d, d)
N =     [150. 0  0;
        0  90  0;
        0  0  70;]
h = 300.
T = 2

for t in 1:T 
    StageCoefficient[t] = StageData(c1[t], c2[t], ū, h, N, s₀, penalty)
end




##########################################################################################
############################  To generate random variable  ###############################
##########################################################################################
T = 2
N_rv = Vector{Int64}()  # the number of realization of each stage
num_Ω = 4
N_rv = [num_Ω for t in 1:T]  ## xxxx 需要update


Random.seed!(12345)

Ω = Dict{Int64,Dict{Int64,RandomVariables}}()   # each stage t, its node are included in Ω[t]
initial_demand = 1e4  #  5.685e8

for t in 1:T 
    Ω[t] = Dict{Int64,RandomVariables}()
    for i in 1:N_rv[t]
        if t == 1
            Ω[t][i]= RandomVariables([initial_demand])
        else
            # Ω[t][i]= RandomVariables( (rand(1)[1]/5+1)*Ω[t-1][i].d )
            Ω[t][i]= RandomVariables( (.3*i+1)*Ω[t-1][i].d )
        end
    end
end




prob = Dict{Int64,Vector{Float64}}()  # P(node in t-1 --> node in t ) = prob[t]
for t in 1:T 
    prob[t] = [0.25 for i in 1:N_rv[t]]
end

# prob[2][1] = .3
# prob[2][2] = .7













scenario_sequence = Dict{Int64, Dict{Int64, Any}}()  ## the first index is for scenario index, the second one is for stage
pathList = Vector{Int64}()
push!(pathList, 1)
P = 1.0

recursion_scenario_tree(pathList, P, scenario_sequence, 2, T = T)
scenario_tree = scenario_sequence


