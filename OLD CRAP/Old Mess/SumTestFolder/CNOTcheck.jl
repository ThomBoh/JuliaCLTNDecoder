using Printf
using PastaQ
using ITensors
using Random
import PastaQ: gate
import PastaQ.ITensors: array
using DelimitedFiles
using Test
ITensors.disable_warn_order()
#@set_warn_order 1000
#Use DoMCCLTN to do monte carlo simulation for ML decoding of repetition code

#To compare results, use the function Pfail in PyCLTNQuimb.py which does an
#exact calculation.

#single qubit gate failure location
function gate(::GateName"pz"; p::Number)
  return [1-p  p
          p    1-p]
end

# CNOT gate failure -- first qubit is control, second qubit is target
function gate(::GateName"CNOT"; pz1::Number, pz2::Number, pz12::Number)
  return [1-pz1-pz2-pz12  pz12            pz1             pz2
          pz2             pz1             pz12            1-pz1-pz2-pz12
          pz1             pz2             1-pz1-pz2-pz12  pz12
          pz12            1-pz1-pz2-pz12  pz2             pz1]
end

#function gate(::GateName"CNOT"; pz1::Number, pz2::Number, pz12::Number)
#  return [1-pz1-pz2-pz12  pz2             pz1             pz12
#          pz12            pz1             pz2             1-pz1-pz2-pz12
#          pz1             pz12            1-pz1-pz2-pz12  pz2
#          pz2             1-pz1-pz2-pz12  pz12            pz1]
#end

#X-type stabilizer ancilla measurement and re-initialization gate -- use pi0 or pi1 depending on
#measurement outcome obtained for that ancilla. The probability of a Z-failure
#on the re-initialized ancilla after the measurement is p.
function gate(::GateName"pi0"; p::Number)
  return [ 1-p 0
           p   0]
end

function gate(::GateName"pi1"; p::Number)
  return [ 0   1-p
           0   p]
end

#This gate "resets" Z-type stabilizer anicllae after measurement. In the biased
#error models that we're simulating, we're only actually interested in phase-flip
#noise and therefore results of X-type stabilizer measurements. So, this gate
#essentially traces out the the degrees of freedom for the Z-type stabilizer measurement
#outcomes that we don't care about.

function gate(::GateName"pizres"; p::Number)
  return [ 1-p 1-p
           p   p]
end

function gate(::GateName"pizresfinal"; p::Number)
  return [ 1 1
           1 1]
end


function CNcheck()
    gates=Tuple[]
    push!(gates,("CNOT",(2,1),(pz1=10,pz2=100,pz12=1000,)))
    str=["Z-","Z+"]

    instate=productstate(2,str)
    MPS=runcircuit(instate,gates)

    oot=PastaQ.array(MPS)
    println(oot)

end

CNcheck()
