using Printf
using PastaQ
using ITensors
using Random
import PastaQ: gate
import PastaQ.ITensors: array
using DelimitedFiles

function gate(::GateName"1q")
    return [1 1
            1 1]
end

function gate(::GateName"2q")
    return [1 1 1 1
            1 1 1 1
            1 1 1 1
            1 1 1 1]
end

function gate(::GateName"3q")
    return [1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1]
end

function gate(::GateName"4q")
    return [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
end

function check()

    circ1=Tuple[]
    push!(circ1,("1q",1))
    TOut=runcircuit(circ1; process =  true)

    circ2=Tuple[]
    push!(circ2,("2q",(1,2)))
    TOut=runcircuit(circ2; process =  true)

    circ3=Tuple[]
    push!(circ3,("3q",(1,2,3)))
    TOut=runcircuit(circ3; process =  true)

    circ4=Tuple[]
    push!(circ4,("4q",(1,2,3,4)))
    TOut=runcircuit(circ4; process =  true)

end

check()
