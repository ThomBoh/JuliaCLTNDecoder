using Printf
using PastaQ
using ITensors
using Random
import PastaQ: gate
import PastaQ.ITensors: array
using Test


ket=productstate(1,["Z+",])

println(ket[1])

ket=productstate(1,["Z-",])

println(ket[1])

ket=productstate(1,["X+",])

println(ket[1])

ket=productstate(1,["X-",])

println(ket[1])



#U[]
