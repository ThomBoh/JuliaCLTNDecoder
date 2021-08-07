using Printf
using PastaQ
using ITensors
using Random
import PastaQ: gate
import PastaQ.ITensors: array
using Test

n = 7
circuit = randomcircuit(n,4)

U = runcircuit(circuit; process = true)
sites = firstsiteinds(U)
@show U
# qubit to be "traced" over
x = 4


tracetensor = ITensor([1,1], sites[x])


M = ITensor[]
for ix in 1:length(U)
    if ix == x
        T = U[ix] * tracetensor
        T = T * tracetensor'
        M[ix-1] = M[ix-1] * T
    else
        push!(M, U[ix])
    end
end
newU = MPO(M)
@show newU
for ix in x:length(U)-1
    newU[ix] = replacetags(newU[ix],"n=$(ix+1)", "n=$(ix)")
end
@show newU

@test length(newU) == length(U)-1
for ix in 1:length(newU)
    @test hastags(newU[ix],"n=$(ix)" )
end






#U[]
