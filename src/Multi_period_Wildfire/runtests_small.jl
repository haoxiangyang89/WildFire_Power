function prepareIndexSets(network_data::Dict{String, Any} ,T::Int64, Ω::Int64; 
    prob_fault::NTuple{6, Float64} = (.1, .22, .09, .05, 0.3, 0.1))
    (pub, pug, pul, pvb, pvg, pvl) = prob_fault

    D = Vector{Int64}()
    G = Vector{Int64}()
    B = Vector{Int64}()

    L = Vector{Tuple{Int64, Int64}}()

    Dᵢ =    Dict{Int64,Vector{Int64}}()
    Gᵢ =    Dict{Int64,Vector{Int64}}()
    out_L = Dict{Int64,Vector{Int64}}()
    in_L =  Dict{Int64,Vector{Int64}}()


    _b = Dict{Tuple{Int64, Int64}, Float64}()  ## total line charging susceptance
    θmax = network_data["branch"]["1"]["angmax"]
    θmin = network_data["branch"]["1"]["angmin"]
    W = Dict{Tuple{Int64, Int64}, Float64}()
    smax = Dict{Int64, Float64}()
    smin = Dict{Int64, Float64}()

    Demand = Dict{Int64, Dict{Int64, Float64}}()
    for t in 1:T 
        Demand[t] = Dict{Int64, Float64}()
    end
    w = Dict{Int64, Float64}()              ## priority level of load D
    cb = Dict{Int64, Float64}()             ## set of fire damage cost cᵢ at :b ∈ B
    cg = Dict{Int64, Float64}()
    cl = Dict{Tuple{Int64, Int64}, Float64}()

    for i in keys(network_data["bus"])
        b = network_data["bus"][i]["bus_i"]
        push!(B, b)
        Dᵢ[b]    = Vector{Int64}()
        Gᵢ[b]    = Vector{Int64}()
        out_L[b] = Vector{Int64}()
        in_L[b]  = Vector{Int64}()
        cb[b] = 300. ## 0                ############# need to revise
    end

    for i in keys(network_data["load"])
        d = network_data["load"][i]["index"]
        b = network_data["load"][i]["load_bus"]
        w[d] = network_data["load"][i]["pd"] * 1e5                     ## priority level of load d

        push!(Dᵢ[b], d)
        push!(D, d)
        for t in 1:T 
            demand = network_data["load"][i]["pd"] * (1 + .05 * t)
            Demand[t][d] = demand
        end
    end

    for i in keys(network_data["gen"])
        g = network_data["gen"][i]["index"]
        b = network_data["gen"][i]["gen_bus"]

        push!(G, g)
        push!(Gᵢ[b], g)

        smax[g] = network_data["gen"][i]["pmax"]
        smin[g] = network_data["gen"][i]["pmin"]
        cg[g] = 200.                               ############# need to revis
    end


    for i in keys(network_data["branch"])
        l = (network_data["branch"][i]["f_bus"], network_data["branch"][i]["t_bus"])

        push!(L, l)
        push!(out_L[l[1]], l[2])
        push!(in_L[l[2]], l[1])

        _b[l] = network_data["branch"][i]["b_fr"]   ## total line charging susceptance
        W[l] = network_data["branch"][i]["rate_a"]              
        cl[l] = 50.                                ############# need to revise
    end

    paramOPF = ParamOPF(_b, θmax, θmin, W, smax, smin)
    indexSets = IndexSets(D, G, L, B ,T, [1:Ω...], Dᵢ, Gᵢ, out_L, in_L)
    paramDemand = ParamDemand(Demand, w, cb, cg, cl, 1e5)


    ## construct random variables
    Ω_rv = Dict{Int64, RandomVariables}()
    for ω in 1:Ω 
        τ = rand(2:T)

        ub = Dict{Int64, Int64}()
        ug = Dict{Int64, Int64}()
        ul = Dict{Tuple{Int64, Int64}, Int64}() 

        vb = Dict{Int64, Int64}()
        vg = Dict{Int64, Int64}()
        vl = Dict{Tuple{Int64, Int64}, Int64}()

        Ibb = Dict{Int64, Vector{Int64}}()
        Ibg = Dict{Int64, Vector{Int64}}()
        Ibl = Dict{Int64, Vector{Tuple{Int64, Int64}}}()
        
        Igb = Dict{Int64, Vector{Int64}}()
        Igg = Dict{Int64, Vector{Int64}}()
        Igl = Dict{Int64, Vector{Tuple{Int64, Int64}}}()

        Ilb = Dict{Tuple{Int64, Int64}, Vector{Int64}}()
        Ilg = Dict{Tuple{Int64, Int64}, Vector{Int64}}()
        Ill = Dict{Tuple{Int64, Int64}, Vector{Tuple{Int64, Int64}}}()


        for b in B
            Ibb[b] = [b]
            Ibg[b] = []
            Ibl[b] = []

            Igb[b] = []
            Igg[b] = []
            Igl[b] = []

            ub[b]  = rand(Binomial(1,pub), 1)[1]
            vb[b]  = rand(Binomial(1,pvb), 1)[1]
        
            for g in Gᵢ[b]
                ug[g]  = rand(Binomial(1,pug), 1)[1]
                vg[g]  = rand(Binomial(1,pvg), 1)[1]
                push!(Ibg[b], g)
            end
        
        end

        
        for l in L
            Ilb[l] = [l[1], l[2]]
            Ilg[l] = []
            Ill[l] = []
            
            ul[l]  = rand(Binomial(1,pul), 1)[1]
            vl[l]  = rand(Binomial(1,pvl), 1)[1]
            (l1, l2) = l 
    
            ## itself
            push!(Ill[l], l)

            ## bus to bus
            push!(Ibb[l1], l2)
            push!(Ibb[l2], l1)

            
            for g in Gᵢ[l1]
                push!(Ibg[l2], g)  ## fire in l2 can impact generators in l1
                push!(Igb[g], l2)  ## fire on generators n l1 can impact generators in l2

                push!(Ilg[l], g)
                push!(Igl[g], l)

                push!(Igg[g], g)
                push!(Igb[g], l1)
            end

            for g in Gᵢ[l2]
                push!(Ibg[l1], g)  ## fire in l1 can impact generators in l2
                push!(Igb[g], l1)  ## fire on generators in l2 can impact generators in l1

                push!(Ilg[l], g)
                push!(Igl[g], l)

                push!(Igg[g], g)
                push!(Igb[g], l2)

                Ibg[g] = unique(Ibg[g])
                Igg[g] = unique(Igg[g])
                Igl[g] = unique(Igl[g])
                Igb[g] = unique(Igb[g])
            end

            ## bus to line
            push!(Ibl[l1], l)
            push!(Ibl[l2], l)
            Ibl[l1] = unique(Ibl[l1])
            Ibl[l2] = unique(Ibl[l2])
    
        end

        Ω_rv[ω] = RandomVariables(τ, ub, ug, ul, vb, vg, vl, Ibb, Ibg, Ibl, Igb, Igg, Igl, Ilb, Ilg, Ill)
    end

    return (indexSets = indexSets, 
            paramOPF = paramOPF, 
            paramDemand = paramDemand, 
            Ω_rv = Ω_rv)
end



network_data = PowerModels.parse_file("/Users/aaron/matpower7.1/data/case30.m")
# display(network_data) # raw dictionary
# PowerModels.print_summary(network_data) # quick table-like summary
# PowerModels.component_table(network_data, "bus", ["vmin", "vmax"]) # component data in matrix form


## construct _prepareIndexSets = prepareIndexSets(D, G, L, B ,3, [1,2,3,4])
T = 5
Ω = 4
pub = .1
pug = .1
pul = .1
pvb = .1
pvg = .1
pvl = .1
prob_fault = (pub, pug, pul, pvb, pvg, pvl) ## prob_faultility of u (v) on a bus (generator, line)
# _prepareIndexSets = prepareIndexSets(network_data, T, Ω; prob_fault = prob_fault)

(indexSets, paramOPF, paramDemand, Ω_rv) = prepareIndexSets(network_data, T, Ω; prob_fault = prob_fault)


prob = Dict{Int64, Float64}()
for ω in indexSets.Ω 
    prob[ω] = .25
end






#############################################################################################################

max_iter = 200; ϵ = 1e-2; Enhanced_Cut = true；

λ_value = .1; Output = 0; Output_Gap = false; Adj = false; Enhanced_Cut = true; threshold = 1e2; 
levelSetMethodParam = LevelSetMethodParam(0.95, λ_value, threshold, 1e14, 3e3, Output, Output_Gap, Adj)



resultDict = SDDiP_algorithm(Ω_rv, prob, 
                    indexSets, 
                    paramDemand, 
                    paramOPF; 
                    levelSetMethodParam = levelSetMethodParam,
                    ϵ = 0.001, M = 1, max_iter = 200, 
                    Enhanced_Cut = true)






using JLD2, FileIO, DataFrames
result_enhanced = copy(Dict(:solHistory => sddipResult, :solution => _LB, :gapHistory => gapList) )
cut_enhanced = copy(cut_collection)

@save "runtests_small2_enhanced.jld2" result_enhanced cut_enhanced
# # @load "runtests_small2_enhanced.jld2" result_enhanced cut_enhanced

# result_LC = copy(Dict(:solHistory => sddipResult, :solution => _LB, :gapHistory => gapList) )
# cut_LC = copy(cut_collection)

# @save "runtests_small2_LC.jld2" result_LC cut_LC
# # @load "runtests_small2_LC.jld2" result_LC cut_LC



# using DataFrames
# using Latexify
# df = DataFrame(A = 'x':'z', B = ["M", "F", "F"])
# latexify(df; env=:table, latex=false)




















