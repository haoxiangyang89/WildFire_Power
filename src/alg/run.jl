using JLD2, FileIO

max_iter = 100; ϵ = 1e-4; 


sddipResult1 = SDDiP_algorithm() 
@save "RTS_test1.jld2" sddipResult1
# @load "RTS_test1.jld2" sddipResult1


sddipResult2 = SDDiP_algorithm() 
@save "RTS_test2.jld2" sddipResult2