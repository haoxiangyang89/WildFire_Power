using JLD2, FileIO

max_iter = 100; ϵ = 1e-2; 


sddipResult1 = SDDiP_algorithm(; ϵ = 1e-4, max_iter = 20) 
@save "RTS_test1.jld2" sddipResult1
# @load "RTS_test1.jld2" sddipResult1


sddipResult2 = SDDiP_algorithm(; ϵ = 1e-3, max_iter = 100) 
@save "RTS_test2.jld2" sddipResult2