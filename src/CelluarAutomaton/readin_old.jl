################################################################################################################################################
######################################################  Data Preparation  ######################################################################
################################################################################################################################################
function prepareIndexSets(  network_data::Dict{String, Any} ,
                                                    T::Int64,
                                                    Ω::Int64
                                                    )
    

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
        cb[b] = rand(1)[1] * 1e6 ## 0                ############# need to revise
    end

    for i in keys(network_data["load"])
        d = network_data["load"][i]["index"]
        b = network_data["load"][i]["load_bus"]
        w[d] = network_data["load"][i]["pd"] * 1e4                     ## priority level of load d

        push!(Dᵢ[b], d)
        push!(D, d)
        for t in 1:T 
            demand = network_data["load"][i]["pd"] * (1 + .2 * t)
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
        cg[g] = rand(1)[1] * 8e5                             ############# need to revis
    end


    for i in keys(network_data["branch"])
        l = (network_data["branch"][i]["f_bus"], network_data["branch"][i]["t_bus"])

        if l ∉ L 
            push!(L, l) 
            push!(out_L[l[1]], l[2])
            push!(in_L[l[2]], l[1])

            _b[l] = 1/network_data["branch"][i]["br_x"]   ## total line charging susceptance
            W[l] = network_data["branch"][i]["rate_a"]              
            cl[l] = rand(1)[1] * 1e4                               ############# need to revise
        end
    end
    

    paramOPF = ParamOPF(_b, θmax, θmin, W, smax, smin)
    indexSets = IndexSets(D, G, unique(L), B ,T, [1:Ω...], Dᵢ, Gᵢ, out_L, in_L)
    paramDemand = ParamDemand(Demand, w, cb, cg, cl, 1e8)
 
     return (indexSets = indexSets, 
             paramOPF = paramOPF, 
             paramDemand = paramDemand)
end

function prepareSimulation(businfo::DataFrame, branchInfo::DataFrame, WFPI_Info::DataFrame;     
                                                            n::Int64 = 100, ## increase boundary
                                                            grid_length::Int64 = 3000
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
        bus_id_location[id] = Int.(ceil.((bus_lla.x - xmin, bus_lla.y - ymin) ./ grid_length)) .+ n
        bus_location_id[bus_id_location[id]] = id
        row_num += 1
    end



    line_id_location = Dict{Int64, Any}()
    line_location_id = Dict{Any, NamedTuple{(:id, :WFPI, :Length), Tuple{Int64, Float64, Union{Float64, Int64}}}}()
    line_id_bus = Dict{Int64, Any}()
    row_num = 1
    for id in 1: nrow(branchInfo)
        length = branchInfo[id, 14]
        WFPI_value = WFPI_Info[id, 3]/5
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
        line_location_id[l1] = (id = id, WFPI = WFPI_value, Length = length)
        line_location_id[l_25] = (id = id, WFPI = WFPI_value, Length = length)
        line_location_id[l_5] = (id = id, WFPI = WFPI_value, Length = length)
        line_location_id[l_75] = (id = id, WFPI = WFPI_value, Length = length)
        line_location_id[l2] = (id = id, WFPI = WFPI_value, Length = length)
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
            line_id_bus = line_id_bus,
            bus_id_location = bus_id_location, 
            bus_location_id = bus_location_id)
end


function prepareScenarios( ;period_span::Int64 = 1, 
                            T::Int64 = 48, 
                            Ω::Int64 = 5, 
                            indexSets::IndexSets = indexSets, prob_fault::Float64 = 0.1, line_id_bus::Dict = line_id_bus)

    function check_vector_exist(key::Union{Tuple{Int64, Int64}, Int64}, dict::Dict)
        if key ∈ keys(dict)
            push!(dict[key], key)
        else
            dict[key] = []
            push!(dict[key], key)
        end
    end


    function check_vector_exist(key::Union{Tuple{Int64, Int64}, Int64}, dict::Dict, G::Dict)
        if length(key) == 1
            if key ∈ keys(dict)
                push!(dict[key], G[key]...)
            else
                dict[key] = []
                push!(dict[key], G[key]...)
            end
        else
            if key ∈ keys(dict)
                push!(dict[key], G[key[1]]...)
                push!(dict[key], G[key[2]]...)
            else
                dict[key] = []
                push!(dict[key], G[key[1]]...)
                push!(dict[key], G[key[2]]...)
            end
        end
    end 


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

    Ω_rv = Dict{Int64, RandomVariables}()
    Gᵢ = indexSets.Gᵢ
    B = indexSets.B
    G = indexSets.G
    L = indexSets.L
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
            ub[b] = 0
            vb[b] = 0
            Ibb[b] = []
            Ibg[b] = []
            Ibl[b] = []
        end

        for g in G 
            ug[g] = 0
            vg[g] = 0
            Igb[g] = []
            Igg[g] = []
            Igl[g] = []
        end

        for l in L 
            ul[l] = 0
            vl[l] = 0
            Ilb[l] = []
            Ilg[l] = []
            Ill[l] = []
        end


        forest = initialize(bus_id_location, line_location_id;  griddims = (x_grid_num, y_grid_num), 
                        environmentInfo = environmentInfo, time_span = 1)

        ## generate wildfire random variables
        for i in 1:floor((T/period_span))
            Agents.step!(forest, agent_step!, wildfire_ignition_step!, period_span)

            if sum(forest.lineFired) > 0 && sum(forest.busFired) > 0
                τ = i * period_span

                if sum(forest.lineFired) > 0
                    for I in findall(isequal(1), forest.lineFired)
                        id = line_location_id[I.I].id
                        vl[line_id_bus[id]] = 1
                        check_vector_exist(line_id_bus[id], Ill)
                        # check_vector_exist(line_id_bus[id], Ilb)
                        check_vector_exist(line_id_bus[id], Ilg, Gᵢ)

                        # push!(Ill[line_id_bus[id]], line_id_bus[id])
                        if line_id_bus[id] ∈ keys(Ilb)
                            push!(Ilb[line_id_bus[id]], line_id_bus[id]...)
                        else 
                            Ilb[line_id_bus[id]] = []
                            push!(Ilb[line_id_bus[id]], line_id_bus[id]...)
                        end
                        # push!(Ilg[line_id_bus[id]], Gᵢ[line_id_bus[id][1]]...)
                        # push!(Ilg[line_id_bus[id]], Gᵢ[line_id_bus[id][2]]...)
                    end
                end
                
                
                if sum(forest.busFired) > 0
                    for I in findall(isequal(1), forest.busFired)
                        bus_id = bus_location_id[I.I]
                        vb[bus_id] = 1
                        for id in Gᵢ[bus_id] 
                            vg[id] = 1 
                            check_vector_exist(id, Igg)
                            check_vector_exist(bus_id, Igb)
                        end
                        check_vector_exist(bus_id, Ibb)
                        check_vector_exist(bus_id, Ibg, Gᵢ)
                    end
                end 

                if sum(forest.lineFault) > 0
                    for I in findall(isequal(1), forest.lineFault)
                        id = line_location_id[I.I].id
                        ul[line_id_bus[id]] = 1
                    end
                end

                if sum(forest.busFault) > 0
                    for I in findall(isequal(1), forest.busFault)
                        bus_id = bus_location_id[I.I]
                        ub[bus_id] = 1
                    end
                end

                break
            end
        end

        ## generate fault random variables
        for l in L 
            if ul[l] == 1
                ub[l[1]] = 1
                ub[l[2]] = 1
            end
        end


        for b in B 
            if ub[b] == 1
                for g in Gᵢ[b]
                    ug[g] = rand(Binomial(1, .8), 1)[1]
                end
            end
        end

        


        Ω_rv[ω] = RandomVariables(τ, ub, ug, ul, vb, vg, vl, Ibb, Ibg, Ibl, Igb, Igg, Igl, Ilb, Ilg, Ill)
    end

    return Ω_rv
end