using Printf
using PastaQ
using ITensors
using Random
import PastaQ: gate
import PastaQ.ITensors: array
using Test

bd=100
n = 14
circuit = randomcircuit(n,4)

U = runcircuit(circuit; process = true)
sites = firstsiteinds(U)


M = ITensor[]
j=0
for ix in 1:length(U)
    if ix%2==0
        tracetensor=ITensor([1,1],sites[ix])
        T = U[ix] * tracetensor
        T = T * tracetensor'
        M[ix-1-j] = M[ix-1-j] * T
        global j=j+1
    else
        push!(M, U[ix])
    end
end
newU = MPO(M)

println("old times")
psi=productstate(U)
@time A=noprime(*(U,psi;cutoff=1e-15,maxdim=bd))
@time A=noprime(*(U,A;cutoff=1e-15,maxdim=bd))
@time A=noprime(*(U,A;cutoff=1e-15,maxdim=bd))
@time A=noprime(*(U,A;cutoff=1e-15,maxdim=bd))


println("new times")
psi=productstate(newU)
@time A=noprime(*(newU,psi;cutoff=1e-15,maxdim=bd))
@time A=noprime(*(newU,A;cutoff=1e-15,maxdim=bd))
@time A=noprime(*(newU,A;cutoff=1e-15,maxdim=bd))
@time A=noprime(*(newU,A;cutoff=1e-15,maxdim=bd))
