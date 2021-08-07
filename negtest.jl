using Printf
using PastaQ
using ITensors
using Random
import PastaQ: gate
import PastaQ.ITensors: array
using Test

sites = siteinds("S=1/2",1)

ket=productstate(sites,["Z+",])

bra0=MPS(sites,[1,0])
bra0=MPS(sites,[0,1])

plz=inner(ket,bra0)
println(plz)
plz=inner(ket,bra1)
println(plz)


ket=productstate(sites,["Z-",])

bra0=MPS(sites,[1,0])
bra0=MPS(sites,[0,1])


plz=inner(ket,bra0)
println(plz)
plz=inner(ket,bra1)
println(plz)


ket=productstate(sites,["X+",])

bra0=MPS(sites,[1,0])
bra0=MPS(sites,[0,1])


plz=inner(ket,bra0)
println(plz)
plz=inner(ket,bra1)
println(plz)

ket=productstate(sites,["X-",])

bra0=MPS(sites,[1,0])
bra0=MPS(sites,[0,1])


plz=inner(ket,bra0)
println(plz)
plz=inner(ket,bra1)
println(plz)





#U[]
