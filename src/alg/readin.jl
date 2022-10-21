################################################################################################################################################
######################################################  Data Preparation  ######################################################################
################################################################################################################################################


function prepareIndexSets(  network_data::Dict{String, Any},
                                                    T::Int64,
                                                    Ω::Int64; 
                                                    branchInfo::DataFrame = branchInfo
                                                    )
    

    D = Vector{Int64}()
    G = Vector{Int64}()
    B = Vector{Int64}()

    L = Vector{Tuple{Int64, Int64}}() 
    multiLines = Dict() ## record the number of lines between two nodes (especailly for double lines)

    Dᵢ =    Dict{Int64,Vector{Int64}}()
    Gᵢ =    Dict{Int64,Vector{Int64}}()
    out_L = Dict{Int64,Vector{Int64}}()
    in_L =  Dict{Int64,Vector{Int64}}()


    _b = Dict{Tuple{Int64, Int64}, Float64}()                   ## total line charging susceptance
    θmax = 100
    θmin = - 100
    W = Dict{Tuple{Int64, Int64}, Float64}()
    smax = Dict{Int64, Float64}()
    smin = Dict{Int64, Float64}()

    Demand = Dict{Int64, Dict{Int64, Float64}}()
    for t in 1:T 
        Demand[t] = Dict{Int64, Float64}()
    end
    w = Dict{Int64, Float64}()                                  ## priority level of load D
    cb = Dict{Int64, Float64}()                                 ## set of fire damage cost cᵢ at :b ∈ B
    cg = Dict{Int64, Float64}()
    cl = Dict{Tuple{Int64, Int64}, Float64}()

    for i in keys(network_data["bus"])
        b = network_data["bus"][i]["bus_i"]
        push!(B, b)
        Dᵢ[b]    = Vector{Int64}()
        Gᵢ[b]    = Vector{Int64}()
        out_L[b] = Vector{Int64}()
        in_L[b]  = Vector{Int64}()
        cb[b] = 5
    end

    for i in keys(network_data["load"])
        d = network_data["load"][i]["index"]
        b = network_data["load"][i]["load_bus"]
        w[d] = round(100^rand(Uniform(0.1, 1)))                         ## priority level of load d

        push!(Dᵢ[b], d)
        push!(D, d)
        for t in 1:T 
            demand = network_data["load"][i]["pd"]
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
        cg[g] = wsample([50, 1000, 2500], [0.2, .75, 0.05], 1)[1]                  
    end


    for i in keys(network_data["branch"])
        l = (network_data["branch"][i]["f_bus"], network_data["branch"][i]["t_bus"])

        if l ∉ L 
            multiLines[l] = 1
            push!(L, l) 
            push!(out_L[l[1]], l[2])
            push!(in_L[l[2]], l[1])

            _b[l] = - 1/network_data["branch"][i]["br_x"]                          ## total line charging susceptance
            W[l] = network_data["branch"][i]["rate_a"]         
            cl[l] = 0.285 *  branchInfo[parse(Int64,i), :Length] 
        else
            multiLines[l] = multiLines[l] + 1
            _b[l] = _b[l] - 1/network_data["branch"][i]["br_x"] 
            W[l] = W[l] + network_data["branch"][i]["rate_a"] 
            cl[l] = cl[l] + 0.285 *  branchInfo[parse(Int64,i), :Length] 
        end

    end
    

    paramOPF = ParamOPF(_b, θmax, θmin, W, smax, smin)
    indexSets = IndexSets(D, G, L, B ,T, [1:Ω...], Dᵢ, Gᵢ, out_L, in_L)
    paramDemand = ParamDemand(Demand, w, cb, cg, cl, 1e4)
 
     return (indexSets = indexSets, 
             paramOPF = paramOPF, 
             paramDemand = paramDemand, 
             multiLines = multiLines)
end

function prepareSimulation(businfo::DataFrame, branchInfo::DataFrame, WFPI_Info::DataFrame;     
                                                            n::Int64 = 100,  ## increase boundary
                                                            grid_length::Int64 = 1000
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
        WFPI_value = WFPI_Info[id, 3]/sum(WFPI_Info[:,3])   ## normalized WFPI
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
                            indexSets::IndexSets = indexSets, line_id_bus::Dict = line_id_bus,
                            bus_id_location::Dict = bus_id_location, line_location_id::Dict = line_location_id
                            )
                
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

    Ω_rv = Dict{Int64, RandomVariables}();
    Gᵢ = indexSets.Gᵢ;
    B = indexSets.B;
    G = indexSets.G;
    L = indexSets.L;
    in_L = indexSets.in_L;
    out_L = indexSets.out_L;
    for ω in 1:Ω 
        τ = 24

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
                        environmentInfo = environmentInfo, time_span = 1);

        ## generate wildfire random variables
        # disruption_not_occur = true
        for i in 1:floor((T/period_span))
            Agents.step!(forest, agent_step!, wildfire_ignition_step!, period_span)
            # sum(forest.ignition)
            # sum(forest.busFired)
            # sum(forest.lineFired)
            if sum(forest.lineFired) > 0 && sum(forest.busFired) > 0 
                τ = i * period_span
                if sum(forest.lineFired) > 0
                    for I in findall(isequal(1), forest.lineFired)
                        id = line_location_id[I.I].id
                        vl[line_id_bus[id]] = 1
                    end
                end
                
                
                if sum(forest.busFired) > 0
                    for I in findall(isequal(1), forest.busFired)
                        bus_id = bus_location_id[I.I]
                        vb[bus_id] = 1
                    end
                end 

                if sum(forest.lineFault) > 0
                    for I in findall(isequal(1), forest.lineFault)
                        id = line_location_id[I.I].id
                        ul[line_id_bus[id]] = rand(Binomial(1, .4), 1)[1]
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
                ub[l[1]] = rand(Binomial(1, .1), 1)[1]
                ub[l[2]] = rand(Binomial(1, .1), 1)[1]
            end
        end


        for b in B 
            if ub[b] == 1
                for g in Gᵢ[b]
                    ug[g] = rand(Binomial(1, .2), 1)[1]
                end
            end
        end


        firedBus = []
        for b in keys(vb)
            if vb[b] == 1
                push!(firedBus, b)
            end 
        end


        # to generate the random sets
        (x_grid_num2, y_grid_num2, line_location_id2, 
                                    line_id_location2, 
                                    line_id_bus2,
                                    bus_id_location2, 
                                    bus_location_id2) = prepareSimulation(businfo, branchInfo, WFPI_Info; n = 10, grid_length = 5000);

        for firedID in firedBus 
            firedPos = bus_id_location2[firedID]
            forest_single_component = initialize_single_component!(firedPos, bus_id_location2, line_location_id2;  
                                                                    griddims = (x_grid_num2, y_grid_num2), 
                                                                    environmentInfo = environmentInfo, 
                                                                    time_span = 1, line_id_location = line_id_location2)
            Agents.step!(forest_single_component, agent_step!, wildfire_ignition_step!, 3)
            if sum(forest_single_component.lineFired) > 0
                for I in findall(isequal(1), forest_single_component.lineFired)
                    push!(Ibl[firedID], line_id_bus2[line_location_id2[I.I].id])
                end
            end

            if sum(forest_single_component.busFired) > 0
                for I in findall(isequal(1), forest_single_component.busFired)
                    push!(Ibb[firedID], bus_location_id2[I.I])
                end
            end 

        end

        for l in L 
            if vl[l] == 1
                push!(Ilb[l], l[rand(1:end)]) 
                
                if !isempty(out_L[l[1]])
                    b2 = out_L[l[1]][rand(1:end)]
                    push!(Ill[l], (l[1], b2))
                end

                if !isempty(out_L[l[2]])
                    b2 = in_L[l[2]][rand(1:end)]
                    push!(Ill[l], (b2, l[2]))
                end
            end
        end


        for b in firedBus 
            for g in Gᵢ[b]
                if rand(Binomial(1, .3), 1)[1] == 1
                    push!(Ibg[b], g)
                    push!(Igb[g], b)
                end 

                for g2 in Gᵢ[b]
                    if rand(Binomial(1, .3), 1)[1] == 1
                        push!(Igg[g], g2)
                    end 
                end

                if rand(Binomial(1, .3), 1)[1] == 1
                    if !isempty(out_L[b])
                        b2 = out_L[b][rand(1:end)]
                        push!(Igl[g], (b, b2))
                    end
                end
            end

        end


        Ω_rv[ω] = RandomVariables(τ, ub, ug, ul, vb, vg, vl, Ibb, Ibg, Ibl, Igb, Igg, Igl, Ilb, Ilg, Ill)
    end

    return Ω_rv
end

