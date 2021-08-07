using Printf
using PastaQ
using ITensors
import PastaQ: gate
import PastaQ.ITensors: array

function gate(::GateName"pz"; p::Number)
  return [1-p  p
          p    1-p]
end

# CNOT gate
function gate(::GateName"CNOT"; pz1::Number, pz2::Number, pz12::Number)
  return [1-pz1-pz2-pz12  pz2             pz1             pz12
          pz12            pz1             pz2             1-pz1-pz2-pz12
          pz1             pz12            1-pz1-pz2-pz12  pz2
          pz2             1-pz1-pz2-pz12  pz12            pz1]
end


function CLTNRepHWFull(nq,p)

    #First initialize all hardware noise parameters. Note that p=k1/k2
    #these are just different probabilities for different failures that I'll plug into the gates

    pxm=0.004 #X readout failure probability
    kappa2=10^7 #2-photon dissipation rate
    alpha2=8 #photon number
    p1=10*p #fail prob for data qubit wait location during very first ancilla initialization
    p2=0.299*(p^0.5) #fail prob for data qubit wait location during execution of a CNOT on other qubits
    p3=p*kappa2*alpha2*(350*(10^(-9)) +(10/(kappa2*alpha2))) #fail prob for data qubit wait location during readout of ancilla
    pz1=0.845*(p^0.5) #prob of z1 otimes I failure after CNOT (qubit 1 is control, 2 is target)
    pz2=0.133*(p^0.5) #prob of I otimes z2 failure after CNOT (qubit 1 is control, 2 is target)
    pz12=0.133*(p^0.5) #prob of z1 otimes z2 failure after CNOT
    pai=15*p/2 #ancilla qubit initialization failure

    #array of data qubits
    qub = [1:nq;]

    #empty circuit
    gates=Tuple[]

    #build initial layer
    for q in qub
        push!(gates,("pz",q,(p=p1,))) #data qubit
        if q < nq
            push!(gates,("pz",nq+q,(p=pai,))) #ancilla qubit (there should only be n-1 of these)
        end
    end

    #build first layer of CNOTs
    for q in qub
        if q == nq
            push!(gates,("pz",nq,(p=p2,))) #wait location for final data qubit, which doesn't participate in CNOT
        else
            push!(gates,("CNOT",(nq+q,q),(pz1=pz1,pz2=pz2,pz12=pz12,)))
        end
    end

    #build second layer of CNOTs
    for q in qub
        if q == 1
            push!(gates,("pz",1,(p=p2,))) #first data qubit wait (no CNOT)
        else
            push!(gates,("CNOT",(nq+q-1,q),(pz1=pz1,pz2=pz2,pz12=pz12,)))
        end
    end

    #final wait locations
    for q in qub
        push!(gates,("pz",q,(p=p3,))) #data
        if q < nq
            push!(gates,("pz",nq+q,(p=pxm,))) #ancilla
        end
    end


    N=2*nq-1
    T=qubits(N)
    TOut=runcircuit(T, gates,cutoff=1e-10,maxdim=5000)

    return TOut
end

#produce MPS and contract it
T=CLTNRepHWFull(7,0.00001)
T=prod(T)

#print out a negative element :'[
println(T[1,2,1,1,2,2,1,1,2,2,2,1,1])
