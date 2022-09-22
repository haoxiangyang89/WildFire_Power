function sparsity(a::Vector)
    total = size(a)[1]
    B = []
    for i=1:total
        if a[i] != 0
            push!(B,i)
        end
    end
    return size(B)[1]/total * 100
end

using Distributions
function t_test(x; conf_level=0.95)
    alpha = (1 - conf_level)
    tstar = quantile(TDist(length(x)-1), 1 - alpha/2)
    SE = std(x)/sqrt(length(x))

    lo, hi = mean(x) .+ [-1, 1] * tstar * SE
    "($lo, $hi)"
end


sparsityList = [] 
angleHorizontal = []
diversityDict = Dict()
for j in 1:Ω    
    for i in keys(cut_collection[1].v) 
        a = cut_collection[j].πb[i][1]
        b = cut_collection[j].πg[i][1]
        c = cut_collection[j].πl[i][1]
        d = vcat(a,b,c)
        for k in keys(cut_collection[1].v)
            a = cut_collection[j].πb[k][1]
            b = cut_collection[j].πg[k][1]
            c = cut_collection[j].πl[k][1]
            d2 = vcat(a,b,c) 
            diversityDict[j, i, k] = round(abs(sum(d .* d2) + 1)/((sqrt(sum(d.^2)+1)) * sqrt(sum(d2.^2)+1)), digits = 5)
        end
        push!(sparsityList, sparsity(d))
        push!(angleHorizontal, 1/sqrt(sum(d.^2)+1) )
    end
end

t_test(sparsityList)
t_test(acos.(angleHorizontal) .* 180 / π)
t_test(acos.(collect(values(diversityDict))) .* 180 / π)

