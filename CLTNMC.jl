using Printf
using PastaQ
using ITensors
using Random
import PastaQ: gate
import PastaQ.ITensors: array


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
#in the |+> state. Ancilla |+> state initialization failure probability is built-in
gate(::GateName"pi01"; p::Number, b::Number) =
  [ (1-b)*(1-p) (1-b)*p
    b*(1-p)     b*(p)]



#This function builds the measurement circuit for a nq-qubit repetition code
#with nr rounds of syndrome measurement and uses PastaQ's runcircuit to turn it into an MPS.
#Note that in order to be able to re-use ancillas and keep qubit connectivity of
#the circuit under control, measurement results for the first nr-1 rounds of ancilla
#measurements (provided in 2d array mbits) are pre-inserted into the circuit.
#The output MPS has 2*nq-1 legs -- the first nq specify the accumulated error on
#the data qubits and the remaining nq-1 legs specify the measurement results
#for the final (nr'th) round of measurements on the ancilla qubits.

#acc specifies cutoff, and bd specifies max bond dimension for the MPS

function CLTNRepHWFull(nq,nr,p,mbits,acc,bd)

    #First initialize all hardware noise parameters. Note that p=k1/k2

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

    qub = [1:nq;]
    rep = [1:nr;]

    gates=Tuple[]


    for r in rep
        for q in qub
            if r==1
                push!(gates,("pz",2*q-1,(p=p1,)))
            elseif r>1
                push!(gates,("pz",2*q-1,(p=0,)))
            end
            if q < nq
                if r==1
                    push!(gates,("pz",2*q,(p=pai,)))
                elseif r>1
                    push!(gates,("pi01",2*q,(p=pai,b=mbits[r-1,q],)))
                end
            end
        end
        for q in qub
            if q == nq
                push!(gates,("pz",2*q-1,(p=p2,)))
            else
                push!(gates,("CNOT",(2*q,2*q-1),(pz1=pz1,pz2=pz2,pz12=pz12,)))
            end
        end
        for q in qub
            if q == 1
                push!(gates,("pz",q,(p=p2,)))
            else
                push!(gates,("CNOT",(2*(q-1),2*q-1),(pz1=pz1,pz2=pz2,pz12=pz12,)))
            end
        end
        for q in qub
            if q < nq
                push!(gates,("pz",2*q,(p=pxm,)))
            end
        end
        for q in qub
            push!(gates,("pz",2*q-1,(p=p3,)))
        end
    end

    #println(gates)
    N=2*nq-1
    T=qubits(N)
    TOut=runcircuit(T,gates,cutoff=acc,maxdim=bd)

    return TOut
end

#The following function takes an MPS (TNin), number of qubits (nq), number of measurement
#rounds (nr), pure error on data qubits (dat), logical error on data qubits (l),
#and measurement outcomes for the final round of ancilla measurements (anc)
#and determines if Maximum Likelihood Decoding either succeeds or fails for that
#error configuration

function lemaxfinal(TNin,nq,nr,dat,l,anc)
    b0=[1. 0.]
    b1=[0. 1.]

    #first build iTensors of ancilla qubits with final round measurement results
    #contracted in

    if anc[1]==0
        TNm=TNin[2]*ITensor(b0,firstsiteind(TNin,2))
    elseif anc[1]==1
        TNm=TNin[2]*ITensor(b1,firstsiteind(TNin,2))
    end

    for q in [2:nq-1;]
        abit=anc[q]
        if abit==0
            TNm=TNm*TNin[2*q]*ITensor(b0,firstsiteind(TNin,2*q))
        elseif abit==1
            TNm=TNm*TNin[2*q]*ITensor(b1,firstsiteind(TNin,2*q))
        else
            println("control flow problem!!")
        end
    end

    #next, contract in bit configuration on data qubits corresponding to no
    #logical Z on data to obtain relative probability of that event

    TNi=TNm
    for q in [1:Int((nq-1)/2);]
        dbit=(dat[q]+l)%2

        if dbit==0

            TNi=TNi*TNin[2*q-1]*ITensor(b0,firstsiteind(TNin,2*q-1))
        elseif dbit==1

            TNi=TNi*TNin[2*q-1]*ITensor(b1,firstsiteind(TNin,2*q-1))
        else
            println("control flow problem!!")
        end

    end

    TNi=TNi*TNin[nq]*ITensor(b0,firstsiteind(TNin,nq))

    for q in [Int(((nq-1)/2))+1:nq-1;]
        dbit=(dat[q]+l)%2

        if dbit==0

            TNi=TNi*TNin[2*(q+1)-1]*ITensor(b0,firstsiteind(TNin,2*(q+1)-1))
        elseif dbit==1

            TNi=TNi*TNin[2*(q+1)-1]*ITensor(b1,firstsiteind(TNin,2*(q+1)-1))
        else
            println("control flow problem!!")
        end

    end

    #next, contract in bit configuration on data qubits corresponding to
    #logical Z on data to obtain relative probability of that event

    TNz=TNm
    for q in [1:Int((nq-1)/2);]
        dbit=(dat[q]+l)%2

        if dbit==0

            TNz=TNz*TNin[2*q-1]*ITensor(b1,firstsiteind(TNin,2*q-1))
        elseif dbit==1

            TNz=TNz*TNin[2*q-1]*ITensor(b0,firstsiteind(TNin,2*q-1))
        else
            println("control flow problem!!")
        end
    end
    TNz=TNz*TNin[nq]*ITensor(b1,firstsiteind(TNin,nq))

    for q in [Int(((nq-1)/2))+1:nq-1;]
        dbit=(dat[q]+l)%2

        if dbit==0

            TNz=TNz*TNin[2*(q+1)-1]*ITensor(b1,firstsiteind(TNin,2*(q+1)-1))
        elseif dbit==1

            TNz=TNz*TNin[2*(q+1)-1]*ITensor(b0,firstsiteind(TNin,2*(q+1)-1))
        else
            println("control flow problem!!")
        end
    end

    #determine which was more likely given supplied error pattern: no logical z,
    #or logical z?

    if TNi[] > TNz[]
        out=0
    else
        out=1
    end

    #return result

    return out

end

#the following function takes number of data qubits nq, number of syndrome measurement
#repetitions nr, and p=k1/k2, and randomly generates an accumulated error pattern on
#the data and ancilla qubits from the distribution determined by the measurement
#circuit and the circuit error model.

function MCCLErr(nq,nr,p)

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

#the next three functions convert binary arrays to decimal numbers used for indexing
#in various ways

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

#this function uses Monte Carlo sampling to estimate the average logical failure
#rate of Maximum Likelihood decoding using the tensor network decoder. nq is number
#of data qubits, nr is number of syndrome measurement repetitions, p=k1/k2 controls
#hardware failure rates, cut is the cutoff parameter for runcircuit, bd is Maximum
#bond dimension for runcircuit, err is the desired standard error percentage for the reported
#failure rate, and buff controls the size of the buffer of decoding results to store
#in an array in order to speed up the common decoding cases.

function DoMCCLTN(nq,nr,p,cut,bd,err,buff)
    #FullMPS=CLTNRepHWFull(nq,nr,p,cut,bd)
    nfail=0
    ntr=0
    pm=0
    breakflag=0
    pct=0
    mu=-1
    edge=1
    bufferkeys=zeros(Int,buff)-ones(Int,buff)
    bufferl=zeros(Int,buff)-ones(Int,buff)
    win=0
    while breakflag==0
        #println("HI!")
        data, anc = MCCLErr(nq,nr,p)
        measno=b2dec2(nr,nq,anc)
        #intanc=b2dec(nq-1,anc[nr,:])
        #if intanc>0
        #    println("anc is")
        #    println(anc)
        #end
        #println(data)
        l=data[Int(((nq-1)/2)+1)]
        deleteat!(data,Int(((nq-1)/2)+1))
        b2d=b2dec1(nq-1,data,l)
        #if b2d>0
        #    println("data is")
        #    println(data)
        #end

        errno=measno+b2d*2^(nr*(nq-1))
        #println(bufferkeys)
        si=findall( x -> x == errno, bufferkeys )
        #println(errno)
        sil=length(si)

        if sil==0
            #println(anc[nr,:])
            MPS=CLTNRepHWFull(nq,nr,p,anc,cut,bd)
            println("made a new MPS")
            ml=lemaxfinal(MPS,nq,nr,data,l,anc[nr,:])
            if edge<buff+1
                bufferkeys[edge]=errno
                bufferl[edge]=ml
                edge=edge+1
            end
        elseif sil>0
            #println(si)
            si=si[1]
            ml=bufferl[si]
            win=win+1
        end

        #ml=lemaxlike(FullMPS,TNm,b2d,nq,nr)

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

        mu=nfail/ntr
        #println(mu)
        if (ntr>1)&&(nfail>4)
            stdev=sqrt((nfail*(1-mu)*(1-mu) + (ntr-nfail)*(mu)*(mu))/(ntr-1))
            stderr=stdev/sqrt(ntr)
            if stderr<(err*mu)
                breakflag=1
            end
        end

        if (pct>10000) && (nfail>4)
            pct=0
            println("mu is")
            println(mu)
            #println(nfail/ntr)
            println("number of trials is")
            println(ntr)
            println("number of failures is")
            println(nfail)
            println("stderr is")
            println(stderr)
            println("target is")
            println(err*mu)
            println("short circuit rate is")
            println(win/ntr)
        end
    end
    println(nfail)
    println(breakflag)
    println(nfail/ntr)
    return mu

end
