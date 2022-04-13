using PowerModels
using Ipopt
using Distributions 

function prepareIndexSets(network_data::Dict{String, Any} ,T::Int64, Ω::Int64; 
    prob_fault::NTuple{6, Float64} = (.1, .22, .09, .05, 0.3, 0.1))
    (pub, pug, pul, pvb, pvg, pvl) = prob_fault

    D = Vector{Int64}()
    G = Vector{Int64}()
    B = Vector{Int64}()

    L = Vector{Tuple{Int64, Int64}}()

    _D =    Dict{Int64,Vector{Int64}}()
    _G =    Dict{Int64,Vector{Int64}}()
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
    w = Dict{Int64, Float64}()
    cb = Dict{Int64, Float64}()  
    cg = Dict{Int64, Float64}()
    cl = Dict{Tuple{Int64, Int64}, Float64}()

    for i in keys(network_data["bus"])
        b = network_data["bus"][i]["bus_i"]
        push!(B, b)
        _D[b]    = Vector{Int64}()
        _G[b]    = Vector{Int64}()
        out_L[b] = Vector{Int64}()
        in_L[b]  = Vector{Int64}()
        push!(_D[b],0)
        push!(_G[b],0)
        push!(out_L[b],0)
        push!(in_L[b],0)
        cb[b] = 456.            ############# need to revise
    end

    for i in keys(network_data["load"])
        d = network_data["load"][i]["load_bus"]
        w[d] = 1
        cb = Dict{Int64, Float64}()  
        cd = Dict{Int64, Float64}()
        cl = Dict{Tuple{Int64, Int64}, Float64}()
        push!(D, d)
        for t in 1:T 
            demand = network_data["load"][i]["pd"] * (1 + .05 * t)
            Demand[t][d] = demand
        end
    end

    for i in keys(network_data["gen"])
        g = network_data["gen"][i]["gen_bus"]
        push!(G, g)
        smax[g] = network_data["gen"][i]["pmax"]
        smin[g] = network_data["gen"][i]["pmin"]
        cg[g] = 456.              ############# need to revis
    end


    for i in keys(network_data["branch"])
        l = (network_data["branch"][i]["f_bus"], network_data["branch"][i]["t_bus"])
        # f_bus = network_data["branch"][i]["f_bus"]
        # t_bus = network_data["branch"][i]["t_bus"]
        _b[l] = network_data["branch"][i]["b_fr"]  ## total line charging susceptance
        push!(L, l)
        W[l] = network_data["branch"][i]["rate_a"]              
        cl[l] = 456.              ############# need to revise
    end

    for l in L 
        (l1, l2) = l
        
        if out_L[l1][1] == 0 
            deleteat!(out_L[l1], 1)
        end

        if in_L[l2][1] == 0 
            deleteat!(in_L[l2], 1)
        end

        push!(out_L[l1], l2)
        push!(in_L[l2], l1)


        if _D[l1][1] == 0 
            deleteat!(_D[l1], 1)
        end 
        push!(_D[l1], l2)

        if _G[l1][1] == 0 
            deleteat!(_G[l1], 1)
        end 
        push!(_G[l1], l2)
    end

    paramOPF = ParamOPF(_b, θmax, θmin, W, smax, smin)
    indexSets = IndexSets(D, G, L, B ,T, [1:Ω...], _D, _G, out_L, in_L)
    paramDemand = ParamDemand(Demand, w, cb, cg, cl)



    Ω_rv = Dict{Int64, RandomVariables}()
    for ω in 1:Ω 
        τ = 2 

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
            Ibb[b] = []
            Ibg[b] = []
            Ibl[b] = []

            Igb[b] = []
            Igg[b] = []
            Igl[b] = []

            ub[b]  = rand(Binomial(1,pub), 1)[1]
            vb[b]  = rand(Binomial(1,pvb), 1)[1]
            if b ∈ keys(Ibb)
                push!(Ibb[b], b)
            else 
                Ibb[b] = []
                push!(Ibb[b], b)
            end
        end
    
        for g in G

            ug[g]  = rand(Binomial(1,pug), 1)[1]
            vg[g]  = rand(Binomial(1,pvg), 1)[1]
            if g ∈ keys(Igg)
                push!(Igg[g], g)
            else 
                Igg[g] = []
                push!(Igg[g], g)
            end
        end
        
        for l in L
            Ilb[l] = []
            Ilg[l] = []
            Ill[l] = []
            
            ul[l]  = rand(Binomial(1,pul), 1)[1]
            vl[l]  = rand(Binomial(1,pvl), 1)[1]
            (l1, l2) = l 
    
            if l ∈ keys(Ill)
                push!(Ill[l], l)
            else 
                Ill[l] = []
                push!(Ill[l], l)
            end
    
            if l1 ∈ B 
                push!(Ilb[l], l1)
                if l2 ∈ B 
                    push!(Ibb[l1], l2)
                end
                if l2 ∈ G 
                    push!(Ibg[l1], l2)
                end
                push!(Ibl[l1], l)
            elseif l1 ∈ G 
                push!(Ilg[l], l1)
                if l2 ∈ B 
                    push!(Igb[l1], l2)
                end
                if l2 ∈ G 
                    push!(Igg[l1], l2)
                end
                push!(Igl[l1], l)
            end
    
            if l2 ∈ B 
                push!(Ilb[l], l2)
                if l1 ∈ B 
                    push!(Ibb[l2], l1)
                end
                if l1 ∈ G 
                    push!(Ibg[l2], l1)
                end
                push!(Ibl[l2], l)
            elseif l2 ∈ G 
                push!(Ilg[l], l2)
                if l1 ∈ B 
                    push!(Igb[l2], l1)
                end
                if l1 ∈ G 
                    push!(Igg[l2], l1)
                end
                push!(Igl[l2], l)
            end
        end
        Ω_rv[ω] = RandomVariables(τ, ub, ug, ul, vb, vg, vl, Ibb, Ibg, Ibl, Igb, Igg, Igl, Ilb, Ilg, Ill)
    end

    return (indexSets = indexSets, 
            paramOPF = paramOPF, 
            paramDemand = paramDemand, 
            Ω_rv = Ω_rv)
end








network_data = PowerModels.parse_file("/Users/aaron/matpower7.1/data/case30.m")
display(network_data) # raw dictionary
PowerModels.print_summary(network_data) # quick table-like summary
PowerModels.component_table(network_data, "bus", ["vmin", "vmax"]) # component data in matrix form


## construct _prepareIndexSets = prepareIndexSets(D, G, L, B ,3, [1,2,3,4])
T = 3
Ω = 4
prob_fault = (pub, pug, pul, pvb, pvg, pvl) ## prob_faultility of u (v) on a bus (generator, line)
_prepareIndexSets = prepareIndexSets(network_data, T, Ω; prob_fault = prob_fault)






























