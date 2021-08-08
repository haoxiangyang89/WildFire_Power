# define data structure

struct scenInfo
    # scenario information
    scenID :: Int64
    τ :: Int64
    u :: Array{Array{Any,1},1}
    w :: Array{Array{Any,1},1}

    spread :: Dict{Any,Any}
end
