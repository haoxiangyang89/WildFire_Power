# test incidence
include("loadMod.jl");

scenInfo_add = "../data/case13/scens_info.csv"
scenSpread_add = "../data/case13/scens_spread.csv"

scenDict = readScen(scenInfo_add,scenSpread_add)
