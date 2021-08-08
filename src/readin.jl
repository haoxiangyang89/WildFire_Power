# read in wildfire/grid data

function readMP(fileAdd::String)
    # read in the power model
    network_data = parse_file(fileAdd);
    return network_data;
end

function readScen(scenInfo_add::String, scenSpread_add::String)
    # function to read in scenarios information
    scenInfo_raw = readdlm(scenInfo_add, ',');
    mt,nt = size(scenInfo_raw);
    scenList = [i for i in 1:(nt - 1)];
    first_column_i = scenInfo_raw[:,1];
    tInfo = scenInfo_raw[findfirst(isequal("t"), first_column_i),2:nt];
    uInfo = [eval(Meta.parse(item)) for item in scenInfo_raw[findfirst(isequal("u"), first_column_i),2:nt]];
    wInfo = [eval(Meta.parse(item)) for item in scenInfo_raw[findfirst(isequal("w"), first_column_i),2:nt]];

    scenSpread_raw = readdlm(scenSpread_add, ',');
    ms,ns = size(scenSpread_raw);
    first_column_s = scenSpread_raw[:,1];
    componentList = scenSpread_raw[2:ms,1];
    spreadInfo = Dict();
    for ic in componentList
        spreadInfo[ic] = [eval(Meta.parse(item)) for item in scenSpread_raw[findfirst(isequal(ic), first_column_s),2:ns]];
    end

    scenDict = Dict();
    for ω in scenList
        spread_scen = Dict();
        for ic in componentList
            spread_scen[ic] = spreadInfo[ic][ω];
        end
        # for each scenario, construct the data structure
        scenDict[ω] = scenInfo(ω, tInfo[ω], uInfo[ω], wInfo[ω], spread_scen);
    end

    return scenDict;
end
