using Printf
using PastaQ
using ITensors
using Random
import PastaQ: gate
import PastaQ.ITensors: array

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
  return [1-pz1-pz2-pz12  pz2             pz1             pz12
          pz12            pz1             pz2             1-pz1-pz2-pz12
          pz1             pz12            1-pz1-pz2-pz12  pz2
          pz2             1-pz1-pz2-pz12  pz12            pz1]
end

# ancilla measurement and re-initialization gate -- b is the measurement outcome
#of an X-basis measurement performed on the ancilla before it is re-initialized
#in the |+> state. Ancilla |+> state initialization failure probability is built-in as p
gate(::GateName"pi01"; p::Number, b::Number) =
  [ (1-b)*(1-p) (1-b)*p
    b*(1-p)     b*(p)]



#This function builds the measurement circuit for a nq-qubit repetition code
#with nr rounds of syndrome measurement and uses PastaQ's runcircuit to turn it into an MPS.
#Note that in order to be able to re-use ancillas and keep qubit connectivity of
#the circuit under control, measurement results for the first nr-1 rounds of ancilla
#measurements (provided in 2d array mbits) are pre-inserted into the circuit.
#The output MPS has 2*nq-1 legs -- the nq odd-numbered legs are for specifying
#the final accumulated error on the data qubits, and the nq-1 even-numbered legs
#are for specifying the measurement outcomes for the final (nr'th) round of syndrome measurements.

#acc specifies cutoff, and bd specifies max bond dimension for the MPS

# nq --> number of data qubits (distance) for repetition code (this means the circuit acts 2*nq-1 qubits because there are nq-1 ancillas)
# nr --> number of repeated rounds of syndrome measurements
# p --> failure parameter, kappa_1/kappa_2
# mbits[nr][nq-1] --> array containing meausrement outcomes (1 or 0) for all rounds of syndrome measurement
# acc --> accuracy cutoff for runcircuit
# bd --> bond dimension cutoff for runcircuit

function CLTNRepHWFull(nq::Int,nr::Int,p::Float64,mbits::Array{Int},acc::Number,bd::Int)

    #First initialize all hardware noise parameters and failure rates. Note that p=k1/k2

    pxm=0.004 #X readout failure probability
    kappa2=10^7 #2-photon dissipation rate
    alpha2=8 #photon number
    p1=10*p #fail prob for data qubit wait location during very first ancilla initialization
    p2=0.299*(p^0.5) #fail prob for data qubit wait location during execution of a CNOT on other qubits
    p3=p*kappa2*alpha2*(350*(10^(-9)) +(10/(kappa2*alpha2))) #fail prob for data qubit wait location during readout of ancilla
    pz1=0.845*(p^0.5) #prob of z1 otimes I failure after CNOT (qubit 1 is control, 2 is target)
    pz2=0.133*(p^0.5) #prob of I otimes z2 failure after CNOT (qubit 1 is control, 2 is target)
    pz12=0.133*(p^0.5) #prob of z1 otimes z2 failure after CNOT
    pai=15*p/2 #failure rate of |+> state prep of ancilla qubit

    qub = [1:nq;] #enumeration of data qubits
    rep = [1:nr;] #enumeration of number of rounds

    gates=Tuple[]

    #build the circuit out of gates
    #each measurement round consists of four time steps

    for r in rep

        #first time step -- initialization/idling

        for q in qub

            #data qubits
            if r==1
                #for first round, data qubits idle while ancillas initialize
                push!(gates,("pz",2*q-1,(p=p1,)))
            elseif r>1
                #for subsequent rounds use "perfect" location since idle error
                #is pushed into failure rate for previous timestep.
                push!(gates,("pz",2*q-1,(p=0,)))
            end

            #ancilla qubits (remember there are nq-1 of them)
            if q < nq
                if r==1
                    #for first round, standard |+> state prep failure
                    push!(gates,("pz",2*q,(p=pai,)))
                elseif r>1
                    #for subsequent rounds, use "projection" gate that takes into
                    #account meas outcome from previous timestep. Then re-initialize.
                    push!(gates,("pi01",2*q,(p=pai,b=mbits[r-1,q],)))
                end
            end
        end

        #second time step -- CNOTs between data [1,..,nq-1] and all ancilla

        for q in qub
            if q == nq
                #final data qubit idles, no CNOT
                push!(gates,("pz",2*q-1,(p=p2,)))
            else
                #CNOT, ancilla qubit i controls on data qubit i
                push!(gates,("CNOT",(2*q,2*q-1),(pz1=pz1,pz2=pz2,pz12=pz12,)))
            end
        end

        #third time step -- CNOTs between data [2,..,nq] and all ancilla

        for q in qub
            if q == 1
                #first data qubit idles, no CNOT
                push!(gates,("pz",q,(p=p2,)))
            else
                #CNOT, ancilla qubit i controls on data qubit i+1
                push!(gates,("CNOT",(2*(q-1),2*q-1),(pz1=pz1,pz2=pz2,pz12=pz12,)))
            end
        end

        #fourth time step -- data qubits idle, ancillas are measured (imperfectly)

        #ancillas have measurement failure with probability pxm
        for q in qub
            if q < nq
                push!(gates,("pz",2*q,(p=pxm,)))
            end
        end

        #data qubits idle with failure probability p3
        for q in qub
            push!(gates,("pz",2*q-1,(p=p3,)))
        end
    end

    #println(gates)
    N=2*nq-1 #total number of qubits

    T=productstate(N) #all zeros initial state sets initial error configuration to be trivial

    TOut=runcircuit(T,gates,cutoff=acc,maxdim=bd)

    return TOut #return the MPS, whose physical indices are as described at top of function
end



#The following function takes an MPS (TNin), number of data qubits (nq), number of measurement
#rounds (nr), error on data qubits (dat), logical error on data qubits (l),
#and measurement outcomes for the final round of ancilla measurements (anc)
#and determines if Maximum Likelihood Decoding either succeeds or fails for that
#error configuration

#TNin-->MPS from runcircuit
#nq-->number of data qubits
#nr-->number of measurement round repetitions
#dat-->array describing actual error on data qubits (1 means Z error on that qubit, 0 means nothing)
#l-->1 if actual error on data qubits has logical Z error in its LTS decomposition, 0 if no
#(in the case of the repetition code, this is the same as whether there is a Z error on the middle qubit)
#anc-->array describing outcomes for final round of syndrome measurements (of which there are nq-1)

function MLDec(TNin::MPS,nq::Int,nr::Int,dat::Array{Int},l::Int,anc::Array{Int})


    pe=basepe(nq,dat,l) #extract the pure error on the data qubits. In this case
                        #it's dat + l*[1,1,...,1]

    #build two different product states to overlap the MPS with
    #mpsi specifies the pure error on the data qubits and no logical Z op
    #mpsz specifies the pure error on the data qubits and also a logical Z op
    #both also specify the outcomes of the final syndrome measurements

    stringi=buildstring(nq,nr,pe,0,anc)
    stringz=buildstring(nq,nr,pe,1,anc)
    mpsi=productstate(siteinds(TNin),stringi)
    mpsz=productstate(siteinds(TNin),stringz)

    #compute "probability" that this pure error and meas outcomes occurred, and no logical z
    pli=inner(TNin,mpsi)

    #compute "probability" that this pure error and meas outcomes occurred, as well as logical z
    plz=inner(TNin,mpsz)

    if pli>plz
        #if it was more likely that no logical z occurred, decoder does nothing
        ml=0
    else
        #if it was more likely that logical z occurred, decoder applies logical z
        ml=1
    end

    if pli<0

        println("IDENTITY COSET PROBABILITY IS NEGATIVE")
        println(pli)
        println("DATA PURE ERROR")
        println(pe)
        println("SYNDROME")
        println(anc)


    elseif plz < 0
        println("LOGICAL Z COSET PROBABILITY IS NEGATIVE")
        println(plz)
        println("DATA PURE ERROR")
        println(pe)
        println("SYNDROME")
        println(anc)

    end


    return ml


end


#the following function takes number of data qubits nq, number of syndrome measurement
#repetitions nr, and p=k1/k2, and randomly generates an accumulated error pattern on
#the data and ancilla qubits from the distribution determined by the measurement
#circuit and the circuit error model.

#nq-->number of data qubits
#nr-->number of meas repetitions
#p-->noise parameter kappa_1/kappa_2

function MCCLErr(nq::Int,nr::Int,p::Number)

    pxm=0.004 #X readout failure probability
    kappa2=10^7 #2-photon dissipation rate
    alpha2=8 #photon number
    p1=10*p #fail prob for data qubit wait location during very first ancilla initialization
    p2=0.299*(p^0.5) #fail prob for data qubit wait location during execution of a CNOT on other qubits
    p3=p*kappa2*alpha2*(350*(10^(-9)) +(10/(kappa2*alpha2))) #fail prob for data qubit wait location during readout of ancilla
    pz1=0.845*(p^0.5) #prob of z1 otimes I failure after CNOT (qubit 1 is control, 2 is target)
    pz2=0.133*(p^0.5) #prob of I otimes z2 failure after CNOT (qubit 1 is control, 2 is target)
    pz12=0.133*(p^0.5) #prob of z1 otimes z2 failure after CNOT
    pai=15*p/2

    data=zeros(Int,nq)
    anc=zeros(Int,nr,nq-1)

    #t=0

    for q in [1:nq;]
        r=rand(Float64)

        if r<p1
            data[q]=(data[q]+1)%2
        end
        if q<nq
            r=rand(Float64)
            if r<pai
                anc[1,q]=(anc[1,q]+1)%2
            end
        end
    end

    #t=1

    for q in [1:nq;]
        if q<nq
            r=rand(Float64)
            if r<pz1
                anc[1,q]=(anc[1,q]+data[q]+1)%2
            elseif (r>pz1)&&(r<(pz1+pz2))
                anc[1,q]=(anc[1,q]+data[q])%2
                data[q]=(data[q]+1)%2

            elseif (r>(pz1+pz2))&&(r<(pz1+pz2+pz12))

                anc[1,q]=(anc[1,q]+data[q]+1)%2
                data[q]=(data[q]+1)%2

            else
                anc[1,q]=(anc[1,q]+data[q])%2
            end
        else
            r=rand(Float64)
            if r<p2

                data[q]=(data[q]+1)%2
            end
        end
    end

    #t=2

    for q in [1:nq;]
        if q<nq
            r=rand(Float64)
            if r<pz1

                anc[1,q]=(anc[1,q]+data[q+1]+1)%2
            elseif (r>pz1)&&(r<(pz1+pz2))

                anc[1,q]=(anc[1,q]+data[q+1])%2
                data[q+1]=(data[q+1]+1)%2

            elseif (r>(pz1+pz2))&&(r<(pz1+pz2+pz12))

                anc[1,q]=(anc[1,q]+data[q+1]+1)%2
                data[q+1]=(data[q+1]+1)%2

            else
                anc[1,q]=(anc[1,q]+data[q+1])%2
            end
        else
            r=rand(Float64)
            if r<p2

                data[1]=(data[1]+1)%2
            end
        end
    end

    #t=3

    for q in [1:nq;]
        r=rand(Float64)
        if r<p3

            data[q]=(data[q]+1)%2
        end
        if q<nq
            r=rand(Float64)
            if r<pxm

                anc[1,q]=(anc[1,q]+1)%2
            end
        end
    end

    if nr>1
        for rep in [2:nr;]

            #t=0

            for q in [1:nq;]
                r=rand(Float64)
                if r<0
                    data[q]=(data[q]+1)%2
                end
                if q<nq
                    r=rand(Float64)
                    if r<pai
                        anc[rep,q]=(anc[rep,q]+1)%2
                    end
                end
            end

            #t=1
            for q in [1:nq;]
                if q<nq
                    r=rand(Float64)
                    if r<pz1
                        anc[rep,q]=(anc[rep,q]+data[q]+1)%2
                    elseif (r>pz1)&&(r<pz1+pz2)
                        anc[rep,q]=(anc[rep,q]+data[q])%2
                        data[q]=(data[q]+1)%2
                    elseif (r>pz1+pz2)&&(r<pz1+pz2+pz12)
                        anc[rep,q]=(anc[rep,q]+data[q]+1)%2
                        data[q]=(data[q]+1)%2
                    else
                        anc[rep,q]=(anc[rep,q]+data[q])%2
                    end
                else
                    r=rand(Float64)
                    if r<p2
                        data[q]=(data[q]+1)%2
                    end
                end
            end

            #t=2
            for q in [1:nq;]
                if q<nq
                    r=rand(Float64)
                if r<pz1
                        anc[rep,q]=(anc[rep,q]+data[q+1]+1)%2
                    elseif (r>pz1)&&(r<pz1+pz2)
                        anc[rep,q]=(anc[rep,q]+data[q+1])%2
                        data[q+1]=(data[q+1]+1)%2
                    elseif (r>pz1+pz2)&&(r<pz1+pz2+pz12)
                        anc[rep,q]=(anc[rep,q]+data[q+1]+1)%2
                        data[q+1]=(data[q+1]+1)%2
                    else
                        anc[rep,q]=(anc[rep,q]+data[q+1])%2
                    end
                else
                    r=rand(Float64)
                    if r<p2
                        data[1]=(data[1]+1)%2
                    end
                end
            end


            #t=3

            for q in [1:nq;]
                r=rand(Float64)
                if r<p3
                    data[q]=(data[q]+1)%2
                end
                if q<nq
                    r=rand(Float64)
                    if r<pxm
                        anc[rep,q]=(anc[rep,q]+1)%2
                    end
                end
            end

        end

    end
    return data, anc
end

#the next three helper functions convert binary arrays to decimal numbers used for indexing
#in various ways. I don't use them in this iteration of the code.

function b2dec(nq,dat)
    out=0
    for i in [1:nq;]
        out=out+(dat[i])*2^(nq-1 - (i-1))
    end
    return out
end

function b2dec1(nq,dat,l)
    out=0
    for i in [1:nq;]
        out=out+((dat[i]+l)%2)*2^(nq-1 - (i-1))
    end
    return out
end

function b2dec2(nr,nq,anc)
    out=0
    for i in [1:nr;]
        for j in [1:(nq-1);]
            out=out+(anc[i,j])*2^(nr*(nq-1)-1 - ((i-1)*(nq-1)+j-1))
        end
    end
    return out
end

#the following function takes an error configuration on the data qubits and returns
#an error configuration containing ONLY the pure error componenet of the input error
#any logical error is removed

#nq-->number of data qubits
#dat-->error on data qubits
#l-->logical error component on data qubits (same as middle bit of dat)

function basepe(nq::Int,dat::Array{Int},l::Int)
    out=zeros(nq)
    #println(dat)
    for i in [1:nq;]
        out[i]=(dat[i]+l)%2
    end
    return out
end

#the following function takes error configurations on data and ancilla qubits
#and returns an array of strings that can be used to create the appropriate
#MPS to overlap with the circuit tensor network in order to retrieve the desired
#element of the MPS

function buildstring(nq::Int,nr::Int,pe::Array{Float64},z::Int,anc::Array{Int})

    if (pe[1]+z)%2==0
        out=["Z+"]
    else
        out=["Z-"]
    end

    for i in [2:nq;]

        if anc[i-1]==0
            push!(out,"Z+")
        else
            push!(out,"Z-")
        end

        if (pe[i]+z)%2==0
            push!(out,"Z+")
        else
            push!(out,"Z-")
        end

    end

    return out

end

#this function uses Monte Carlo sampling to estimate the average logical failure
#rate of Maximum Likelihood decoding using the tensor network decoder. nq is number
#of data qubits, nr is number of syndrome measurement repetitions, p=k1/k2 controls
#hardware failure rates, cut is the cutoff parameter for runcircuit, bd is Maximum
#bond dimension for runcircuit, err is the desired standard error percentage for the reported
#failure rate, and buff controls the size of the buffer of decoding results to store
#in an array in order to speed up the common decoding cases.

#nq-->number of data qubits (distance of repetition code)
#nr-->number of repetitions of syndrome measurements
#p-->noise parameter kappa_1/kappa_2
#cut-->MPS accuracy cutoff
#bd-->MPS bond dimension cutoff
#err-->threshold for standard error on logical failure rate, as percentage of reported value

function DoMCCLTN(nq::Int,nr::Int,p::Float64,cut::Float64,bd::Int,err::Float64)

    nfail=0 #number of decoder failures
    ntr=0 #number of MC trials
    breakflag=0 #flag for while loop
    mu=0 #average failure rate
    pct=0 #counter for deciding when to print interim results

    #pm=0
    #pct=0
    #mu=-1
    #edge=1
    #bufferkeys=zeros(Int,buff)-ones(Int,buff)
    #bufferl=zeros(Int,buff)-ones(Int,buff)
    #win=0

    while breakflag==0

        #generate data error pattern and syndrome history from distribution
        data, anc = MCCLErr(nq,nr,p)

        #middle bit indicates logical operator in data's LTS decomposition
        l=data[Int(((nq-1)/2)+1)]

        #build MPS from runcircuit for this syndrome history
        MPS=CLTNRepHWFull(nq,nr,p,anc,cut,bd)

        #do maximum likelihood decoding using this MPS and error
        ml=MLDec(MPS,nq,nr,data,l,anc[nr,:])

        #ignore the bunch of commented stuff below, it's there to enable caching
        #decoder results for speedup

        #measno=b2dec2(nr,nq,anc)

        #intanc=b2dec(nq-1,anc[nr,:])
        #if intanc>0
        #    println("anc is")
        #    println(anc)
        #end
        #println(data)

        #pe1=basepe(nq,data,l)
        #deleteat!(pe1,Int(((nq-1)/2)+1))
        #b2d=b2dec1(nq-1,pe1,0)

        #if b2d>0
        #    println("data is")
        #    println(data)
        #end

        #errno=measno+b2d*2^(nr*(nq-1))
        #println(bufferkeys)
        #si=findall( x -> x == errno, bufferkeys )
        #println(errno)
        #sil=length(si)

        #if sil==0
            #println(anc[nr,:])
            #MPS=CLTNRepHWFull(nq,nr,p,anc,cut,bd)
            #println(typeof(MPS))
            #ml=lemaxfinalmem(MPS,nq,nr,data,l,anc[nr,:])
            #if edge<buff+1
                #bufferkeys[edge]=errno
                #bufferl[edge]=ml
                #edge=edge+1
            #end
        #elseif sil>0
            #println(si)
            #si=si[1]
            #ml=bufferl[si]
            #win=win+1
        #end

        #ml=lemaxlike(FullMPS,TNm,b2d,nq,nr)


        #decide if the decoder guessed the right correction
        if Int(l)==Int(ml)
            #println("succ")
            ntr=ntr+1
            pct=pct+1
        else
            #println("fail")
            ntr=ntr+1
            nfail=nfail+1
            pct=pct+1
        end

        mu=nfail/ntr #average logical failure rate for decoder

        #println(mu)

        #calculate standard error and decide if it's below the set threshold
        if (ntr>1)&&(nfail>4)
            stdev=sqrt((nfail*(1-mu)*(1-mu) + (ntr-nfail)*(mu)*(mu))/(ntr-1))
            stderr=stdev/sqrt(ntr)
            if stderr<(err*mu)
                breakflag=1
            end
        end

        #print out statistics every 1000 trials until threshold is reached
        if (pct>1000) && (nfail>4)
            pct=0
            println("mu is")
            println(mu)
            println(nfail/ntr)
            println("number of trials is")
            println(ntr)
            println("number of failures is")
            println(nfail)
            println("stderr is")
            println(stderr)
            println("target is")
            println(err*mu)
            #println("short circuit rate is")
            #println(win/ntr)
        end

    end
    println(nfail)
    println(ntr)
    println(nfail/ntr)
    return mu

end


DoMCCLTN(3,3,3e-4,1e-15,10,0.001)
