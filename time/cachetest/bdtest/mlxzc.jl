using Printf
using PastaQ
using ITensors
using Random
import PastaQ: gate
import PastaQ.ITensors: array
using DelimitedFiles

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


#create schedules of CNOTs between ancilla qubits and data

function schedulemaker(dz::Integer,dx::Integer)

    #figure out how many z checks per row

    zcounts=zeros(Int,dz+1)
    for i in [1:(dz+1);]
        if i%2==1
            zcounts[i]=floor((dx-1)/2)
        else
            zcounts[i]=floor((dx)/2)
        end
    end

    #build z ancilla qubit schedule

    zsched=[]

    #i indexes which row of checks
    for i in [1:size(zcounts)[1];]
        #j indexes which check in row i
        for j in [1:zcounts[i];]

            if i==1

                sch=[]

                push!(sch,[-1,-1])
                push!(sch,[-1,-1])
                push!(sch,[1,2*(j-1)+2])
                push!(sch,[1,2*(j-1)+3])

                push!(zsched,sch)

            elseif i%2==0

                sch=[]

                if i==(size(zcounts)[1])

                    push!(sch,[i-1,2*(j-1)+1])
                    push!(sch,[i-1,2*(j-1)+2])
                    push!(sch,[-1,-1])
                    push!(sch,[-1,-1])

                else

                    push!(sch,[i-1,2*(j-1)+1])
                    push!(sch,[i-1,2*(j-1)+2])
                    push!(sch,[i,2*(j-1)+1])
                    push!(sch,[i,2*(j-1)+2])
                end

                push!(zsched,sch)

            elseif i%2==1

                sch=[]

                if i==(size(zcounts)[1])

                    push!(sch,[i-1,2*(j-1)+2])
                    push!(sch,[i-1,2*(j-1)+3])
                    push!(sch,[-1,-1])
                    push!(sch,[-1,-1])

                else

                    push!(sch,[i-1,2*(j-1)+2])
                    push!(sch,[i-1,2*(j-1)+3])
                    push!(sch,[i,2*(j-1)+2])
                    push!(sch,[i,2*(j-1)+3])
                end

                push!(zsched,sch)
            end
        end
    end

    #figure out how many x checks per column
    xcounts=zeros(Int,dx+1)
    for i in [1:(dx+1);]
        if i%2==1
            xcounts[i]=floor(dz/2)
        else
            xcounts[i]=floor((dz-1)/2)
        end
    end

    #build x ancilla qubit schedule

    xsched=[]

    #i indexes which row of checks
    for i in [1:size(xcounts)[1];]
        #j indexes which check in row i
        for j in [1:xcounts[i];]

            if i==1

                sch=[]

                push!(sch,[-1,-1])
                push!(sch,[2*(j-1)+1,1])
                push!(sch,[-1,-1])
                push!(sch,[2*(j-1)+2,1])

                push!(xsched,sch)

            elseif i%2==0

                sch=[]

                if i==(size(xcounts)[1])

                    push!(sch,[2*(j-1)+2,i-1])
                    push!(sch,[-1,-1])
                    push!(sch,[2*(j-1)+3,i-1])
                    push!(sch,[-1,-1])

                else

                    push!(sch,[2*(j-1)+2,i-1])
                    push!(sch,[2*(j-1)+2,i])
                    push!(sch,[2*(j-1)+3,i-1])
                    push!(sch,[2*(j-1)+3,i])
                end

                push!(xsched,sch)

            elseif i%2==1

                sch=[]

                if i==(size(xcounts)[1])

                    push!(sch,[2*(j-1)+1,i-1])
                    push!(sch,[-1,-1])
                    push!(sch,[2*(j-1)+2,i-1])
                    push!(sch,[-1,-1])

                else

                    push!(sch,[2*(j-1)+1,i-1])
                    push!(sch,[2*(j-1)+1,i])
                    push!(sch,[2*(j-1)+2,i-1])
                    push!(sch,[2*(j-1)+2,i])
                end

                push!(xsched,sch)
            end
        end
    end


    return zsched, xsched
end

#create schedule of when data qubits on the boundary have wait locations instead of CNOTs

function bdrysched(dz::Integer,dx::Integer)

    bsch=zeros(Int,4,dz,dx)

    for iz in [1:dz;]

        for ix in [1:dx;]

            if (iz==1) && (ix==1)

                bsch[4,iz,ix]=1
                bsch[3,iz,ix]=1

            elseif iz==(dz) && (ix==1)

                if dz%2==0
                    bsch[1,iz,ix]=1
                    bsch[2,iz,ix]=1
                else
                    bsch[2,iz,ix]=1
                    bsch[4,iz,ix]=1
                end

            elseif (ix==(dx)) && (iz==1)

                if dx%2==0

                    bsch[4,iz,ix]=1
                    bsch[3,iz,ix]=1

                else

                    bsch[1,iz,ix]=1
                    bsch[3,iz,ix]=1
                end

            elseif (ix==dx) && (iz==dz)

                if (dx+dz)%2==0
                    bsch[1,iz,ix]=1
                    bsch[2,iz,ix]=1
                else
                    bsch[1,iz,ix]=1
                    bsch[3,iz,ix]=1
                end

            elseif (iz==1) && ((ix>1) && (ix<dx))

                if ix%2==1
                    bsch[3,iz,ix]=1
                else
                    bsch[4,iz,ix]=1
                end

            elseif (ix==1) && ((iz>1) && (iz<dz))

                if iz%2==1
                    bsch[4,iz,ix]=1
                else
                    bsch[2,iz,ix]=1
                end


            elseif (iz==dz) && ((ix>1) && (ix<dx))

                if dz%2==0

                    if ix%2==1
                        bsch[1,iz,ix]=1
                    else
                        bsch[2,iz,ix]=1
                    end
                else
                    if ix%2==1
                        bsch[2,iz,ix]=1
                    else
                        bsch[1,iz,ix]=1
                    end
                end

            elseif (ix==dx) && ((iz>1) && (iz<dz))

                if dx%2==0

                    if iz%2==1
                        bsch[3,iz,ix]=1
                    else
                        bsch[1,iz,ix]=1
                    end
                else
                    if iz%2==1
                        bsch[1,iz,ix]=1
                    else
                        bsch[3,iz,ix]=1
                    end
                end
            end
        end
    end

    return bsch
end

#reorder the x ancilla qubits for easier 1D layout

function xreorder(xsch,dz,dx)
    #assume dz,dx both odd
    #then there are dx+1 "columns" of (dz-1)/2 x checks each
    #and dz-1 "rows" of (dx+1)/2 x checks each
    nxc=size(xsch)[1]
    cpc=(dz-1)/2
    cpr=(dx+1)/2
    neworder=zeros(Int,nxc)
    for i in [1:nxc;]
        k=2*cpc
        c=floor((i-1)/k)+1 #this tells us which pair of columns it lives in
        j=(i-1)%k
        j=j+1 #this tells us which one it is in that pair of columns

        if j<=cpc
            #it's in the first column

            r=2*j-1 #is the row it lives in
            inew=(r-1)*cpr+c


        else
            #it's in the second column
            r=2*(j-cpc) #is the row it lives in
            inew=(r-1)*cpr+c

        end
        neworder[i]=inew
    end

    xschedout=[]

    for i in [1:nxc;]
        blank=[]
        for j in [1:4;]
            push!(blank,[0,0])
        end
        push!(xschedout,blank)

    end

    for i in [1:nxc;]
        xschedout[neworder[i]]=xsch[i]
    end

    return xschedout

end

#create 1D layout of data and ancilla qubits, this uses a simple "back and forth" layout

function SurfLayout(dz::Integer,dx::Integer,zsch,xsch)


    nzc=size(zsch)[1]
    nxc=size(xsch)[1]
    zcpr=(dx-1)/2
    xcpr=(dz+1)/2
    qtot=2*dz*dx-1


    count=0
    mode=0
    zc=1
    xc=1
    drow=1
    layout=[]
    while count<(2*dz*dx-1)

        if mode==0
            #first/last row of z checks, type 1
            for i in [1:zcpr;]

                push!(layout,[1,[zc,zc]])
                zc=zc+1
                count=count+1
            end
            mode=1

        elseif mode==1
            #dx data qubits, type 0
            for i in [1:dx;]

                push!(layout,[0,[drow,i]])
                count=count+1
            end
            drow=drow+1
            if drow>dz
                mode=0
            else
                mode=2
            end

        elseif mode==2
            #dx checks, starting with an x check
            for i in [1:dx;]
                if i%2==1
                    push!(layout,[2,[xc,xc]])
                    xc=xc+1
                    count=count+1
                else
                    push!(layout,[1,[zc,zc]])
                    zc=zc+1
                    count=count+1
                end
            end
            mode=1
        end

    end
    return layout

end

#optimize the "back and forth" layout slightly

function optlayout(layout,dz,dx)

    nrows=2*dz+1
    rowcount=0
    rowmode=0
    zcpr=convert(Int64,(dx-1)/2)
    xcpr=convert(Int64,(dz+1)/2)
    optlay=[]
    q=0
    while rowcount<(2*dz+1)

        if rowcount==0
            #initial check layer -- z checks only, of which there are zcpr

            for i in [1:zcpr;]

                push!(optlay,layout[i])
                q=q+1
            end

            rowcount=rowcount+1
            rowmode=(rowmode+1)%4

        elseif rowcount==(2*dz)

            for i in [1:zcpr;]

                push!(optlay,layout[q+i])

            end
            q=q+2
            rowcount=rowcount+1
            rowmode=(rowmode+1)%4

        elseif rowmode==0

            for i in [1:dx;]

                push!(optlay,layout[q+i])

            end
            q=q+dx
            rowcount=rowcount+1
            rowmode=(rowmode+1)%4

        elseif rowmode==1

            for i in [1:dx;]

                push!(optlay,layout[q+dx-(i-1)])

            end
            q=q+dx
            rowcount=rowcount+1
            rowmode=(rowmode+1)%4

        elseif rowmode==2

            for i in [1:dx;]

                push!(optlay,layout[q+i])

            end
            q=q+dx
            rowcount=rowcount+1
            rowmode=(rowmode+1)%4

        elseif rowmode==3

            push!(optlay,layout[q+1])

            for i in [1:(dx-1);]

                push!(optlay,layout[q+dx-(i-1)])

            end
            q=q+dx
            rowcount=rowcount+1
            rowmode=(rowmode+1)%4

        end

    end

    return optlay

end

#figure out mappings of where each data and ancilla qubit end up in the 1D linear layout

function linemaps(layout,dz,dx,nzc,nxc)

    qlinemap=zeros(Int,dz,dx)
    zlinemap=zeros(Int,nzc)
    xlinemap=zeros(Int,nxc)
    N=2*dx*dz-1
    for i in [1:N;]

        typ=layout[i][1]

        if typ==0
            #data qubit
            addr=layout[i][2]
            zad=addr[1]
            xad=addr[2]

            qlinemap[zad,xad]=i

        elseif typ==1
            #z check qubit
            addr=layout[i][2][1]

            zlinemap[addr]=i

        elseif typ==2
            #x check qubit
            addr=layout[i][2][1]

            xlinemap[addr]=i

        end

    end

    return qlinemap,zlinemap,xlinemap
end

#simulate a CNOT failure location

function CNOTfail(r,p,al2,tmeas,k2,nth,pmz,pmx)

    pz1=0.91*((p)^(0.5))
    pz2=0.15*((p)^(0.5))
    pz1z2=0.15*((p)^(0.5))
    px1=0.93*((p)^(0.5))*exp(-2*al2)
    px2=0.93*((p)^(0.5))*exp(-2*al2)
    px1x2=0.93*((p)^(0.5))*exp(-2*al2)
    pz1x2=0.93*((p)^(0.5))*exp(-2*al2)
    py1=0.93*((p)^(0.5))*exp(-2*al2)
    py1x2=0.93*((p)^(0.5))*exp(-2*al2)
    py2=0.28*p*exp(-2*al2)
    py1z2=0.28*p*exp(-2*al2)
    px1z2=0.28*p*exp(-2*al2)
    pz1y2=0.28*p*exp(-2*al2)
    py1y2=0.28*p*exp(-2*al2)
    px1y2=0.28*p*exp(-2*al2)

    if r < pz1
        zout=[1,0]
        xout=[0,0]
    elseif (r>pz1) && (r<(pz1+pz2))
        zout=[0,1]
        xout=[0,0]
    elseif (r>(pz1+pz2)) && (r<(pz1+pz2+pz1z2))
        zout=[1,1]
        xout=[0,0]
    elseif (r>(pz1+2*pz2)) && (r<(pz1+2*pz2+px1))
        zout=[0,0]
        xout=[1,0]
    elseif (r>(pz1+2*pz2+px1)) && (r<(pz1+2*pz2+px1+px2))
        zout=[0,0]
        xout=[0,1]
    elseif (r>(pz1+2*pz2+2*px1)) && (r<(pz1+2*pz2+2*px1+px1x2))
        zout=[0,0]
        xout=[1,1]
    elseif (r>(pz1+2*pz2+3*px1)) && (r<(pz1+2*pz2+3*px1+pz1x2))
        zout=[1,0]
        xout=[0,1]
    elseif (r>(pz1+2*pz2+4*px1)) && (r<(pz1+2*pz2+4*px1+py1))
        zout=[1,0]
        xout=[1,0]
    elseif (r>(pz1+2*pz2+5*px1)) && (r<(pz1+2*pz2+5*px1+py1x2))
        zout=[1,0]
        xout=[1,1]
    elseif (r>(pz1+2*pz2+6*px1)) && (r<(pz1+2*pz2+6*px1+py2))
        zout=[0,1]
        xout=[0,1]
    elseif (r>(pz1+2*pz2+6*px1+py2)) && (r<(pz1+2*pz2+6*px1+py2+py1z2))
        zout=[1,1]
        xout=[1,0]
    elseif (r>(pz1+2*pz2+6*px1+2*py2)) && (r<(pz1+2*pz2+6*px1+2*py2+px1z2))
        zout=[0,1]
        xout=[1,0]
    elseif (r>(pz1+2*pz2+6*px1+3*py2)) && (r<(pz1+2*pz2+6*px1+3*py2+pz1y2))
        zout=[1,1]
        xout=[0,1]
    elseif (r>(pz1+2*pz2+6*px1+4*py2)) && (r<(pz1+2*pz2+6*px1+4*py2+py1y2))
        zout=[1,1]
        xout=[1,1]
    elseif (r>(pz1+2*pz2+6*px1+5*py2)) && (r<(pz1+2*pz2+6*px1+5*py2+px1y2))
        zout=[0,1]
        xout=[1,1]
    else
        zout=[0,0]
        xout=[0,0]
    end
    return zout,xout
end

#obtain a single circuit error sample
#obtain a single circuit error sample

function SurfMCError(dz,dx,nr,zsched,xsched,bsch,p,al2,tmeas,k2,nth,pmz,pmx)

    ppax=(15/2)*p
    ppaz=0.39*exp(-4*al2)
    pzip=10*p#*(1+exp(-4*al2))
    pxip=10*p*exp(-4*al2)
    pzic=0.31*((p)^(0.5))#*(1+exp(-4*al2))
    pxic=0.31*((p)^(0.5))*exp(-4*al2)
    pzim=k2*p*al2*(tmeas+(10/(k2*al2)) )*(1+2*nth)#*(1+exp(-4*al2))
    pxim=k2*p*al2*(tmeas+(10/(k2*al2)) )*exp(-4*al2)
    pz1=0.91*((p)^(0.5))#+3*0.93*((p)^(0.5))*exp(-2*al2)
    pz2=0.15*((p)^(0.5))#+3*0.28*p*exp(-2*al2)
    pz1z2=0.15*((p)^(0.5))#+3*0.28*p*exp(-2*al2)
    px1=0.93*((p)^(0.5))*exp(-2*al2)
    px2=0.93*((p)^(0.5))*exp(-2*al2)
    px1x2=0.93*((p)^(0.5))*exp(-2*al2)
    pz1x2=0.93*((p)^(0.5))*exp(-2*al2)
    py1=0.93*((p)^(0.5))*exp(-2*al2)
    py1x2=0.93*((p)^(0.5))*exp(-2*al2)
    py2=0.28*p*exp(-2*al2)
    py1z2=0.28*p*exp(-2*al2)
    px1z2=0.28*p*exp(-2*al2)
    pz1y2=0.28*p*exp(-2*al2)
    py1y2=0.28*p*exp(-2*al2)
    px1y2=0.28*p*exp(-2*al2)

    nzc=size(zsched)[1]
    nxc=size(xsched)[1]

    zancz=zeros(Int,nzc)
    zancx=zeros(Int,nzc)
    xancz=zeros(Int,nxc)
    xancx=zeros(Int,nxc)
    dataz=zeros(Int,dz,dx)
    datax=zeros(Int,dz,dx)

    zcmeas=zeros(Int,nr,nzc)
    xcmeas=zeros(Int,nr,nxc)


    for rep in [1:nr;]

        #t=0
        for zc in [1:nzc;]
            r=rand(Float64)
            if r<pzip
                #z check ancilla gets z error
                zancz[zc]=(zancz[zc]+1)%2
            #elseif (r>pzip) && (r<(pzip+pxip))
                #z check ancilla gets x error
            #    zancx[zc]=(zancx[zc]+1)%2
            elseif (r>(pzip)) && (r<(pzip+pxip))
                #z check ancilla gets y error
                zancz[zc]=(zancz[zc]+1)%2
                zancx[zc]=(zancx[zc]+1)%2
            end
        end

        for xc in [1:nxc;]
            r=rand(Float64)
            if r<pzip
                #x check ancilla gets z error
                xancz[xc]=(xancz[xc]+1)%2
            #elseif (r>pzip) && (r<(pzip+pxip))
                #x check ancilla gets x error
            #    xancx[xc]=(xancx[xc]+1)%2
            elseif (r>(pzip)) && (r<(pzip+pxip))
                #x check ancilla gets y error
                xancz[xc]=(xancz[xc]+1)%2
                xancx[xc]=(xancx[xc]+1)%2
            end
        end

        for zq in [1:dz;]
            for xq in [1:dx;]
                if rep==1
                    r=rand(Float64)
                    if r<pzip
                        #data qubit gets z error
                        dataz[zq,xq]=(dataz[zq,xq]+1)%2
                    #elseif (r>pzip) && (r<(pzip+pxip))
                        #data qubit gets x error
                    #    datax[zq,xq]=(datax[zq,xq]+1)%2
                    elseif (r>(pzip)) && (r<(pzip+pxip))
                        #data qubit gets y error
                        dataz[zq,xq]=(dataz[zq,xq]+1)%2
                        datax[zq,xq]=(datax[zq,xq]+1)%2
                    end
                end
            end
        end

        #t=1-4, CNOTs

        for cn in [1:4;]

            for zc in [1:nzc;]
                #z-type check is to identify x errors on data
                #data is CNOT control, ancilla is target
                qub=zsched[zc][cn]
                if qub[1]>-1

                    r=rand(Float64)

                    zf,xf=CNOTfail(r,p,al2,tmeas,k2,nth,pmz,pmx)

                    #z error of the data qubit (control) inherits prior z error from ancilla (target)
                    #as well as the control z error from the CNOT operation
                    dataz[qub[1],qub[2]]=(zf[1]+dataz[qub[1],qub[2]]+zancz[zc])%2
                    #z error of the ancilla qubit (target) inherits target z error from CNOT operation
                    zancz[zc]=(zancz[zc]+zf[2])%2

                    #x error of the ancilla qubit (target) inherits prior x error from data (control)
                    #as well as the target x error from the CNOT operation
                    zancx[zc]=(zancx[zc]+xf[2]+datax[qub[1],qub[2]])%2
                    #x error of the data qubit (control) inherits control x error from CNOT operation
                    datax[qub[1],qub[2]]=(datax[qub[1],qub[2]]+xf[1])%2

                else

                    #ancilla wait

                    r=rand(Float64)

                    if r<pzic
                        zancz[zc]=(zancz[zc]+1)%2
                    #elseif (r>pzic) && (r<(pzic+pxic))
                    #    zancx[zc]=(zancx[zc]+1)%2
                    elseif (r>(pzic)) && (r<(pzic+pxic))
                        zancz[zc]=(zancz[zc]+1)%2
                        zancx[zc]=(zancx[zc]+1)%2
                    end
                end
            end

            for xc in [1:nxc;]
                #x-type check is to identify z errors on data
                #ancilla is CNOT control, data is target
                qub=xsched[xc][cn]
                if qub[1]>-1

                    r=rand(Float64)

                    zf,xf=CNOTfail(r,p,al2,tmeas,k2,nth,pmz,pmx)

                    #z error of the ancilla qubit (control) inherits prior z error from data (target)
                    #as well as the control z error from the CNOT operation
                    xancz[xc]=(xancz[xc]+dataz[qub[1],qub[2]]+zf[1])%2
                    #z error of the data qubit (target) inherits target z error from CNOT operation
                    dataz[qub[1],qub[2]]=(dataz[qub[1],qub[2]]+zf[2])%2


                    #x error of the data qubit (target) inherits prior x error from ancilla (control)
                    #as well as the target x error from the CNOT operation
                    datax[qub[1],qub[2]]=(datax[qub[1],qub[2]]+xancx[xc]+xf[2])%2
                    #x error of the ancilla qubit (control) inherits control x error from CNOT operation
                    xancx[xc]=(xancx[xc]+xf[1])%2

                else

                    #ancilla wait

                    r=rand(Float64)

                    if r<pzic
                        xancz[xc]=(xancz[xc]+1)%2
                    #elseif (r>pzic) && (r<(pzic+pxic))
                    #    xancx[xc]=(xancx[xc]+1)%2
                    elseif (r>(pzic)) && (r<(pzic+pxic))
                        xancz[xc]=(xancz[xc]+1)%2
                        xancx[xc]=(xancx[xc]+1)%2
                    end
                end
            end

            #boundary data qubit wait locations

            for zq in [1:dz;]
                for xq in [1:dx;]

                    if bsch[cn,zq,xq]==1
                        r=rand(Float64)
                        if r<pzic
                            dataz[zq,xq]=(dataz[zq,xq]+1)%2
                        #elseif (r>pzic) && (r<(pzic+pxic))
                        #    datax[zq,xq]=(datax[zq,xq]+1)%2
                        elseif (r>(pzic)) && (r<(pzic+pxic))
                            dataz[zq,xq]=(dataz[zq,xq]+1)%2
                            datax[zq,xq]=(datax[zq,xq]+1)%2
                        end
                    end
                end
            end
        end

        #t=5

        #final data qubit waits
        for zq in [1:dz;]
            for xq in [1:dx;]

                r=rand(Float64)
                if r<pzim
                    dataz[zq,xq]=(dataz[zq,xq]+1)%2
                #elseif (r>pzim) && (r<(pzim+pxim))
                #    datax[zq,xq]=(datax[zq,xq]+1)%2
                elseif (r>(pzim)) && (r<(pzim+pxim))
                    dataz[zq,xq]=(dataz[zq,xq]+1)%2
                    datax[zq,xq]=(datax[zq,xq]+1)%2
                end
            end
        end

        #final ancilla qubit measurement failures
        for zc in [1:nzc;]

            r=rand(Float64)
            if r<pmz
                zancx[zc]=(zancx[zc]+1)%2
            end
        end

        for xc in [1:nxc;]

            r=rand(Float64)
            if r<pmx
                xancz[xc]=(xancz[xc]+1)%2
            end
        end

        #record measurement outcomes and reset ancillas

        zcmeas[rep,:]=zancx
        xcmeas[rep,:]=xancz

        zancx=zeros(Int,nzc)
        zancz=zeros(Int,nzc)
        xancx=zeros(Int,nxc)
        xancz=zeros(Int,nxc)
    end


    return dataz,datax,zcmeas,xcmeas
end

#returns the logical Z and X operators

function SurfLogs(dz,dx)

    Lz=zeros(Int,dz,dx)
    Lx=zeros(Int,dz,dx)

    if dx%2==0

        Lx[dz,:]=ones(Int,dx)

    else

        Lx[1,:]=ones(Int,dx)
    end

    Lz[:,dx]=ones(Int,dz)

    return Lz,Lx
end

#returns the pure error components of the accumulated errors, up to stabilizer

function getSurfPE(Ez,Ex,dz,dx)

    Lz,Lx=SurfLogs(dz,dx)
    LC=LogComp(Ez,Ex,dz,dx)

    PEZ=LC[1]*Lz + Ez
    PEX=LC[2]*Lx + Ex

    for i in [1:dz;]
        for j in [1:dx;]

            PEZ[i,j]=PEZ[i,j]%2
            PEX[i,j]=PEX[i,j]%2

        end
    end

    return PEZ,PEX
end

#returns the logical components of the accumulated errors

function LogComp(Ez,Ex,dz,dx)

    LC=[0,0]

    if dx%2==0
        println("shouldn't see this!!")
        LC[1]=Ez[dz,dx]
        LC[2]=Ex[dz,dx]

    else
        for i in [1:dx;]
            LC[1]=(LC[1]+Ez[1,i])%2
        end
        for i in [1:dz;]
            LC[2]=(LC[2]+Ex[i,dx])%2
        end
    end

    return LC
end





#the following function takes error configurations on data and ancilla qubits
#and returns an array of strings that can be used to create the appropriate
#MPS to overlap with the circuit tensor network in order to retrieve the desired
#element of the MPS

function buildstrings(dz::Int,dx::Int,nr,pez,Synz,Synx,layout)

    N=2*dz*dx-1
    outi=[]
    outz=[]

    count=0


    for i in [1:N;]
        #println(i)
        typ=layout[i][1]
#        println(layout[i])
        #println(typ)
        if typ==0
           #data qubit

            qz=layout[i][2][1]
            qx=layout[i][2][2]

            if qx==dx
                #on logical boundary
                #println(pez)
                #println(qz)
                #println(qx)
                e=pez[qz,qx]
                #println(e)
                if e==0
                    push!(outi,"Z+")
                    push!(outz,"Z-")
                else
                    push!(outi,"Z-")
                    push!(outz,"Z+")
                end
                count=count+1
            else
                #println(pez)
                #println(qz)
                #println(qx)
                e=pez[qz,qx]
                #println(e)
                #e=pez[qz,qx]
                if e==0
                    push!(outi,"Z+")
                    push!(outz,"Z+")
                else
                    push!(outi,"Z-")
                    push!(outz,"Z-")
                end
                count=count+1
            end

        elseif typ==1
            #println("hay1")
            #z check ancilla
            zc=layout[i][2][1]
            if i==1
                #println("hay2")
                if Synz[nr,zc]==0
                    #println("hay3")
                    outi=["X+"]
                    outz=["X+"]
                else
                    #println("hay3")
                    outi=["X+"]
                    outz=["X+"]
                end
                count=count+1
            else
                if Synz[nr,zc]==0
                    push!(outi,"X+")
                    push!(outz,"X+")
                else
                    push!(outi,"X+")
                    push!(outz,"X+")
                end
                count=count+1
            end

        elseif typ==2
            #x check ancilla
            xc=layout[i][2][1]
            if Synx[nr,xc]==0
                push!(outi,"Z+")
                push!(outz,"Z+")
            else
                push!(outi,"Z-")
                push!(outz,"Z-")
            end
            count=count+1
        end

    end
    #println(outi)
    #println(outz)
    return outi,outz

end


function buildinitstring(dz::Int,dx::Int,nr,layout)

    N=2*dz*dx-1
    out=[]

    #println("buildstring layout/pez/synz/synx")
    #println(layout)
    #println(pez)
    #println(Synz)
    #println(Synx)
    #println(Synx[nr,:])
    #println(typeof(layout))
    #println(typeof(pez))
    #println(typeof(Synz))
    #println(typeof(Synx))
    count=0


    for i in [1:N;]
        #println(i)
        typ=layout[i][1]
#        println(layout[i])
        #println(typ)
        if typ==0
           #data qubit
            push!(out,"Z+")
            count=count+1

        elseif typ==1
            #println("hay1")
            #z check ancilla

            if i==1
                out=["Z+"]
                #println("hay2")
                count=count+1
            else
                push!(out,"Z+")
                count=count+1
            end

        elseif typ==2
            push!(out,"Z+")
            #x check ancilla
            count=count+1
        end

    end
    #println(outi)
    #println(outz)
    return out

end


function MLError(B,dz,dx,nr,PEZ,PEX,Synz,Synx,zsch,xsch,bsch,layout,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)

    psi=productstate(B)
    #@show psi
    #println("noprime")
    MPS=noprime(*(B,psi;cutoff=1e-15))
    #@show MPS

    for re in [1:(nr-1);]
        #println("glue")
        MPS=glue(MPS,dz,dx,nr,re,PEZ,PEX,Synx,zsch,xsch,bsch,layout,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)
        #println("noprime")
        MPS=noprime(*(B,MPS;cutoff=1e-15))

    end
    #println("glue")
    MPS=gluefinal(MPS,dz,dx,nr,nr,PEZ,PEX,Synx,zsch,xsch,bsch,layout,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)

    stringi,stringz=buildstrings(dz,dx,nr,PEZ,Synz,Synx,layout)
    mpsi=productstate(siteinds(MPS),stringi)
    mpsz=productstate(siteinds(MPS),stringz)

    #compute "probability" that this pure error and meas outcomes occurred, and no logical z
    pli=inner(MPS,mpsi)

    #compute "probability" that this pure error and meas outcomes occurred, as well as logical z
    plz=inner(MPS,mpsz)
    #println(pli)
    #println(plz)

    if pli>plz
        #if it was more likely that no logical z occurred, decoder does nothing
        ml=0
    else
        #if it was more likely that logical z occurred, decoder applies logical z
        ml=1
    end
    #println(ml)

    return ml

end


function glue(psi,dz,dx,nr,re,PEZ,PEX,Synx,zsch,xsch,bsch,layout,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)

    #First initialize all hardware noise parameters and failure rates. Note that p=k1/k2

    ppax=(15/2)*p
    ppaz=0.39*exp(-4*al2)
    pzip=10*p*(1+exp(-4*al2))
    pxip=10*p*exp(-4*al2)
    pzic=0.31*((p)^(0.5))*(1+exp(-4*al2))
    pxic=0.31*((p)^(0.5))*exp(-4*al2)
    pzim=k2*p*al2*(tmeas+(10/(k2*al2)) )*(1+2*nth)*(1+exp(-4*al2))
    pxim=k2*p*al2*(tmeas+(10/(k2*al2)) )*exp(-4*al2)
    pz1=0.91*((p)^(0.5))+3*0.93*((p)^(0.5))*exp(-2*al2)
    pz2=0.15*((p)^(0.5))+3*0.28*p*exp(-2*al2)
    pz1z2=0.15*((p)^(0.5))+3*0.28*p*exp(-2*al2)
    px1=0.93*((p)^(0.5))*exp(-2*al2)
    px2=0.93*((p)^(0.5))*exp(-2*al2)
    px1x2=0.93*((p)^(0.5))*exp(-2*al2)
    pz1x2=0.93*((p)^(0.5))*exp(-2*al2)
    py1=0.93*((p)^(0.5))*exp(-2*al2)
    py1x2=0.93*((p)^(0.5))*exp(-2*al2)
    py2=0.28*p*exp(-2*al2)
    py1z2=0.28*p*exp(-2*al2)
    px1z2=0.28*p*exp(-2*al2)
    pz1y2=0.28*p*exp(-2*al2)
    py1y2=0.28*p*exp(-2*al2)
    px1y2=0.28*p*exp(-2*al2)

    gates=Tuple[]

    #build the circuit out of gates
    #each measurement round consists of four time steps


    nzc=size(zsch)[1]
    nxc=size(xsch)[1]

    nq=size(layout)[1]


    for q in [1:nq;]

        #first time step -- initialization/idling
        typ=layout[q][1]

        if typ==0
            #data qubit
            push!(gates,("pz",q,(p=0,)))

        elseif typ==1
            #z check ancilla qubit
            if re==(nr-1)
                push!(gates,("pizresfinal",q,(p=0,)))
            else
                push!(gates,("pizres",q,(p=0,)))
            end


        elseif typ==2
            #x check ancilla qubit
            xc=layout[q][2][1]
            bit=Synx[re,xc]
            if bit==0
                push!(gates,("pi0",q,(p=0,)))
            else
                push!(gates,("pi1",q,(p=0,)))
            end

        end

    end

    #istr=buildinitstring(dz,dx,nr,layout)
    #T=productstate(nq) #all zeros initial state sets initial error configuration to be trivial
    #TI=productstate(siteinds(T),istr)
    #display(T)
    psiOut=runcircuit(psi,gates,cutoff=acc,maxdim=bd)

    return psiOut #return the MPS, whose physical indices are as described at top of function
end

function gluefinal(psi,dz,dx,nr,re,PEZ,PEX,Synx,zsch,xsch,bsch,layout,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)

    #First initialize all hardware noise parameters and failure rates. Note that p=k1/k2

    ppax=(15/2)*p
    ppaz=0.39*exp(-4*al2)
    pzip=10*p*(1+exp(-4*al2))
    pxip=10*p*exp(-4*al2)
    pzic=0.31*((p)^(0.5))*(1+exp(-4*al2))
    pxic=0.31*((p)^(0.5))*exp(-4*al2)
    pzim=k2*p*al2*(tmeas+(10/(k2*al2)) )*(1+2*nth)*(1+exp(-4*al2))
    pxim=k2*p*al2*(tmeas+(10/(k2*al2)) )*exp(-4*al2)
    pz1=0.91*((p)^(0.5))+3*0.93*((p)^(0.5))*exp(-2*al2)
    pz2=0.15*((p)^(0.5))+3*0.28*p*exp(-2*al2)
    pz1z2=0.15*((p)^(0.5))+3*0.28*p*exp(-2*al2)
    px1=0.93*((p)^(0.5))*exp(-2*al2)
    px2=0.93*((p)^(0.5))*exp(-2*al2)
    px1x2=0.93*((p)^(0.5))*exp(-2*al2)
    pz1x2=0.93*((p)^(0.5))*exp(-2*al2)
    py1=0.93*((p)^(0.5))*exp(-2*al2)
    py1x2=0.93*((p)^(0.5))*exp(-2*al2)
    py2=0.28*p*exp(-2*al2)
    py1z2=0.28*p*exp(-2*al2)
    px1z2=0.28*p*exp(-2*al2)
    pz1y2=0.28*p*exp(-2*al2)
    py1y2=0.28*p*exp(-2*al2)
    px1y2=0.28*p*exp(-2*al2)

    gates=Tuple[]

    #build the circuit out of gates
    #each measurement round consists of four time steps


    nzc=size(zsch)[1]
    nxc=size(xsch)[1]

    nq=size(layout)[1]

    for rep in [1:1;]
        for q in [1:nq;]


            #first time step -- initialization/idling
            typ=layout[q][1]


            if typ==0
                #data qubit
                push!(gates,("pz",q,(p=pzip,)))

            elseif typ==1
                #z check ancilla qubit
                push!(gates,("pizresfinal",q,(p=0,)))
                #push!(gates,("pz",q,(p=0,)))

            elseif typ==2
                #x check ancilla qubit
                push!(gates,("pz",q,(p=0,)))
                #xc=layout[q][2][1]
                #bit=Synx[re,xc]
                #if bit==0
                #    push!(gates,("pi0",q,(p=0,)))
                #else
                #    push!(gates,("pi1",q,(p=0,)))
                #end

            end


        end


    end


    #istr=buildinitstring(dz,dx,nr,layout)
    #T=productstate(nq) #all zeros initial state sets initial error configuration to be trivial
    #TI=productstate(siteinds(T),istr)
    #display(T)
    psiOut=runcircuit(psi,gates,cutoff=acc,maxdim=bd)

    return psiOut #return the MPS, whose physical indices are as described at top of function
end

function buildblock(dz,dx,nr,zsch,xsch,bsch,layout,ql,zl,xl,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)

    #First initialize all hardware noise parameters and failure rates. Note that p=k1/k2

    ppax=(15/2)*p
    ppaz=0.39*exp(-4*al2)
    pzip=10*p*(1+exp(-4*al2))
    pxip=10*p*exp(-4*al2)
    pzic=0.31*((p)^(0.5))*(1+exp(-4*al2))
    pxic=0.31*((p)^(0.5))*exp(-4*al2)
    pzim=k2*p*al2*(tmeas)*(1+2*nth)*(1+exp(-4*al2))
    pxim=k2*p*al2*(tmeas)*exp(-4*al2)
    pz1=0.91*((p)^(0.5))+3*0.93*((p)^(0.5))*exp(-2*al2)
    pz2=0.15*((p)^(0.5))+3*0.28*p*exp(-2*al2)
    pz1z2=0.15*((p)^(0.5))+3*0.28*p*exp(-2*al2)
    px1=0.93*((p)^(0.5))*exp(-2*al2)
    px2=0.93*((p)^(0.5))*exp(-2*al2)
    px1x2=0.93*((p)^(0.5))*exp(-2*al2)
    pz1x2=0.93*((p)^(0.5))*exp(-2*al2)
    py1=0.93*((p)^(0.5))*exp(-2*al2)
    py1x2=0.93*((p)^(0.5))*exp(-2*al2)
    py2=0.28*p*exp(-2*al2)
    py1z2=0.28*p*exp(-2*al2)
    px1z2=0.28*p*exp(-2*al2)
    pz1y2=0.28*p*exp(-2*al2)
    py1y2=0.28*p*exp(-2*al2)
    px1y2=0.28*p*exp(-2*al2)

    gates=Tuple[]

    #build the circuit out of gates
    #each measurement round consists of four time steps


    nzc=size(zsch)[1]
    nxc=size(xsch)[1]

    nq=size(layout)[1]

    for rep in [1:1;]
        for q in [1:nq;]


            #first time step -- initialization/idling
            typ=layout[q][1]


            if typ==0
                #data qubit
                push!(gates,("pz",q,(p=pzip,)))

            elseif typ==1
                #z check ancilla qubit
                push!(gates,("pz",q,(p=ppaz,)))

            elseif typ==2
                #x check ancilla qubit
                push!(gates,("pz",q,(p=ppax,)))

            end


        end

        for cnt in [1:4;]

            for q in [1:nq;]

                typ=layout[q][1]

                if typ==0
                    #data qubit, need to see if this is a wait location
                    addr=layout[q][2]
                    wait=bsch[cnt,addr[1],addr[2]]
                    if wait==1

                        push!(gates,("pz",q,(p=pzic,)))

                    end

                elseif typ==1

                    #z check ancilla qubit, identifies x errors, data is control, ancilla is target

                    zc=layout[q][2][1] #which z check it is in the list of z check schedules
                    partner=zsch[zc][cnt]

                    if partner==[-1,-1]
                        #wait location for this ancilla
                        push!(gates,("pz",q,(p=pzic,)))
                    else
                        pzad=partner[1]
                        pxad=partner[2]
                        q2=ql[pzad,pxad]
                        push!(gates,("CNOT",(q2,q),(pz1=pz1,pz2=pz2,pz12=pz1z2,)))
                    end

                elseif typ==2

                    #x check ancilla qubit, identifies z errors, data is target, ancilla is control

                    xc=layout[q][2][1] #which x check it is in the list of x check schedules
                    partner=xsch[xc][cnt]

                    if partner[1]==-1
                        #wait location for this ancilla
                        push!(gates,("pz",q,(p=pzic,)))
                    else
                        pzad=partner[1]
                        pxad=partner[2]
                        q2=ql[pzad,pxad]
                        push!(gates,("CNOT",(q,q2),(pz1=pz1,pz2=pz2,pz12=pz1z2,)))
                    end

                end

            end
        end


        for q in [1:nq;]


            #final time step -- idling/measurement
            typ=layout[q][1]

            if typ==0
                #data qubit
                push!(gates,("pz",q,(p=pzim,)))

            elseif typ==1
                #z check ancilla qubit
                push!(gates,("pz",q,(p=0,)))

            elseif typ==2
                #x check ancilla qubit
                push!(gates,("pz",q,(p=pmx,)))

            end

        end
    end


    #istr=buildinitstring(dz,dx,nr,layout)
    #T=productstate(nq) #all zeros initial state sets initial error configuration to be trivial
    #TI=productstate(siteinds(T),istr)
    #display(T)
    TOut=runcircuit(gates,cutoff=acc,maxdim=bd; process =  true)

    return TOut #return the MPS, whose physical indices are as described at top of function
end


function SurfCirc(dz,dx,nr,PEZ,PEX,Synz,Synx,zsch,xsch,bsch,layout,ql,zl,xl,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)

    #First initialize all hardware noise parameters and failure rates. Note that p=k1/k2

    ppax=(15/2)*p
    ppaz=0.39*exp(-4*al2)
    pzip=10*p*(1+exp(-4*al2))
    pxip=10*p*exp(-4*al2)
    pzic=0.31*((p)^(0.5))*(1+exp(-4*al2))
    pxic=0.31*((p)^(0.5))*exp(-4*al2)
    pzim=k2*p*al2*(tmeas+(10/(k2*al2)) )*(1+2*nth)*(1+exp(-4*al2))
    pxim=k2*p*al2*(tmeas+(10/(k2*al2)) )*exp(-4*al2)
    pz1=0.91*((p)^(0.5))+3*0.93*((p)^(0.5))*exp(-2*al2)
    pz2=0.15*((p)^(0.5))+3*0.28*p*exp(-2*al2)
    pz1z2=0.15*((p)^(0.5))+3*0.28*p*exp(-2*al2)
    px1=0.93*((p)^(0.5))*exp(-2*al2)
    px2=0.93*((p)^(0.5))*exp(-2*al2)
    px1x2=0.93*((p)^(0.5))*exp(-2*al2)
    pz1x2=0.93*((p)^(0.5))*exp(-2*al2)
    py1=0.93*((p)^(0.5))*exp(-2*al2)
    py1x2=0.93*((p)^(0.5))*exp(-2*al2)
    py2=0.28*p*exp(-2*al2)
    py1z2=0.28*p*exp(-2*al2)
    px1z2=0.28*p*exp(-2*al2)
    pz1y2=0.28*p*exp(-2*al2)
    py1y2=0.28*p*exp(-2*al2)
    px1y2=0.28*p*exp(-2*al2)

    gates=Tuple[]

    #build the circuit out of gates
    #each measurement round consists of four time steps


    nzc=size(zsch)[1]
    nxc=size(xsch)[1]

    nq=size(layout)[1]

    for rep in [1:nr;]
        for q in [1:nq;]


            #first time step -- initialization/idling
            typ=layout[q][1]

            if rep>1

                if typ==0
                    #data qubit
                    push!(gates,("pz",q,(p=0,)))

                elseif typ==1
                    #z check ancilla qubit
                    zc=layout[q][2][1]
                    #println(zc)
                    #println(rep-1)
                    #println(size(Synz))
                    mbit=Synz[rep-1,zc]
                    if mbit==0
                        if rep<nr
                            push!(gates,("pizres",q,(p=ppaz,)))
                        else
                            push!(gates,("pizresfinal",q,(p=0,)))
                        end

                    else
                        if rep<nr
                            push!(gates,("pizres",q,(p=ppaz,)))
                        else
                            push!(gates,("pizresfinal",q,(p=0,)))
                        end
                    end


                elseif typ==2
                    #x check ancilla qubit
                    xc=layout[q][2][1]
                    mbit=Synx[rep-1,xc]
                    if mbit==0

                        push!(gates,("pi0",q,(p=ppax,)))

                    else
                        push!(gates,("pi1",q,(p=ppax,)))
                    end

                end

            else

                if typ==0
                    #data qubit
                    push!(gates,("pz",q,(p=pzip,)))

                elseif typ==1
                    #z check ancilla qubit
                    push!(gates,("pz",q,(p=ppaz,)))

                elseif typ==2
                    #x check ancilla qubit
                    push!(gates,("pz",q,(p=ppax,)))

                end

            end

        end

        for cnt in [1:4;]

            for q in [1:nq;]

                typ=layout[q][1]

                if typ==0
                    #data qubit, need to see if this is a wait location
                    addr=layout[q][2]
                    wait=bsch[cnt,addr[1],addr[2]]
                    if wait==1

                        push!(gates,("pz",q,(p=pzic,)))

                    end

                elseif typ==1

                    #z check ancilla qubit, identifies x errors, data is control, ancilla is target

                    zc=layout[q][2][1] #which z check it is in the list of z check schedules
                    partner=zsch[zc][cnt]

                    if partner==[-1,-1]
                        #wait location for this ancilla
                        push!(gates,("pz",q,(p=pzic,)))
                    else
                        pzad=partner[1]
                        pxad=partner[2]
                        q2=ql[pzad,pxad]
                        push!(gates,("CNOT",(q2,q),(pz1=pz1,pz2=pz2,pz12=pz1z2,)))
                    end

                elseif typ==2

                    #x check ancilla qubit, identifies z errors, data is target, ancilla is control

                    xc=layout[q][2][1] #which x check it is in the list of x check schedules
                    partner=xsch[xc][cnt]

                    if partner[1]==-1
                        #wait location for this ancilla
                        push!(gates,("pz",q,(p=pzic,)))
                    else
                        pzad=partner[1]
                        pxad=partner[2]
                        q2=ql[pzad,pxad]
                        push!(gates,("CNOT",(q,q2),(pz1=pz1,pz2=pz2,pz12=pz1z2,)))
                    end

                end

            end
        end


        for q in [1:nq;]


            #final time step -- idling/measurement
            typ=layout[q][1]

            if typ==0
                #data qubit
                push!(gates,("pz",q,(p=pzim,)))

            elseif typ==1
                #z check ancilla qubit
                push!(gates,("pz",q,(p=0,)))

            elseif typ==2
                #x check ancilla qubit
                push!(gates,("pz",q,(p=pmx,)))

            end

        end
    end


    istr=buildinitstring(dz,dx,nr,layout)
    T=productstate(nq) #all zeros initial state sets initial error configuration to be trivial
    TI=productstate(siteinds(T),istr)
    #display(T)
    TOut=runcircuit(TI,gates,cutoff=acc,maxdim=bd)

    return TOut #return the MPS, whose physical indices are as described at top of function
end

function SurfMC(dz,dx,nr,p,al2,tmeas,k2,nth,acc,bd,err,nt; sim_id::Int=-1)

    if sim_id < 0
      fname = "experiment_dz$(dz)_dx$(dx)_p$(p).txt"
    else
      fname = "experiment_dz$(dz)_dx$(dx)_p$(p)_id$(sim_id).txt"
    end
    Random.seed!(1234 * (sim_id+2))

    pmz=exp(-1.5-0.9*al2)
    if p==1e-5
        pmx=1.85e-5
    elseif p==2e-5
        pmx=1.9e-5
    elseif p==3e-5
        pmx=1.95e-5
    elseif p==4e-5
        pmx=2.05e-5
    elseif p==5e-5
        pmx=2.1e-5
    elseif p==6e-5
        pmx=2.15e-5
    elseif p==7e-5
        pmx=2.2e-5
    elseif p==8e-5
        pmx=2.25e-5
    elseif p==9e-5
        pmx=2.3e-5
    elseif p==1e-4
        pmx=2.4e-5
    elseif p==2e-4
        pmx=3.1e-5
    elseif p==3e-4
        pmx=3.95e-5
    elseif p==4e-4
        pmx=5.05e-5
    elseif p==5e-4
        pmx=6.6e-5
    elseif p==6e-4
        pmx=8.15e-5
    elseif p==7e-4
        pmx=9.7e-5
    elseif p==8e-4
        pmx=11.7e-5
    elseif p==9e-4
        pmx=14.25e-5
    elseif p==1e-3
        pmx=16.8e-5
    else
        println("bad mmnt error rate!")
    end

    zsch,xsch=schedulemaker(dz,dx)
    xsch=xreorder(xsch,dz,dx)
    bsch=bdrysched(dz,dx)
    nzc=size(zsch)[1]
    nxc=size(xsch)[1]
    layout=SurfLayout(dz,dx,zsch,xsch)
    layout=optlayout(layout,dz,dx)
    ql,zl,xl=linemaps(layout,dz,dx,nzc,nxc)
    #println(zsch)
    #println(xsch)
    Lz,Lx=SurfLogs(dz,dx)
    #display(Lz)
    #display(Lx)
    n=0
    f=0
    fx=0
    breakflag=0
    pct=0
    #println("buildblock")
    B=buildblock(dz,dx,nr,zsch,xsch,bsch,layout,ql,zl,xl,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)

    totime=0
    while n<nt
        #breakflag=1
        #println("pspsps")
        Ez,Ex,Synz,Synx=SurfMCError(dz,dx,nr,zsch,xsch,bsch,p,al2,tmeas,k2,nth,pmz,pmx)
        #display(Ez)
        #display(Synx)
        #display(Ex)
        #display(Synz)

        PEZ,PEX=getSurfPE(Ez,Ex,dz,dx)
        #display(PEZ)
        #display(PEX)
        L=LogComp(Ez,Ex,dz,dx)
        LZ=L[1]
        LX=L[2]
        #display(L)
        #println(PEZ)
        elapsed = @elapsed begin
            MLZ=MLError(B,dz,dx,nr,PEZ,PEX,Synz,Synx,zsch,xsch,bsch,layout,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)
	    end
	    if n>0
	        totime=totime+elapsed
        end

	#PEZ[1,1]=(PEZ[1,1]+1)%2
        #PEZ[1,2]=(PEZ[1,2]+1)%2
        #PEZ[2,1]=(PEZ[2,1]+1)%2
        #PEZ[2,2]=(PEZ[2,2]+1)%2
        #MLZ=MLError(B,dz,dx,nr,PEZ,PEX,Synz,Synx,zsch,xsch,bsch,layout,p,al2,tmeas,k2,nth,pmz,pmx)
        #println(PEZ)
        #MPS=SurfCirc(dz,dx,nr,PEZ,PEX,Synz,Synx,zsch,xsch,bsch,layout,ql,zl,xl,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)
        #MPS
        #MLZ=MLError(MPS,dz,dx,nr,PEZ,PEX,Synz,Synx,zsch,xsch,bsch,layout)
        #display(MLZ)
        #MLZ=MLL[1]
        #print(L)
        #print(MLL)
        #println(L)
        #println(MLZ)
        resz=(LZ+MLZ)%2
        if resz==1
            f=f+1
        elseif LX==1
            fx=fx+1
        end

        n=n+1
        pct=pct+1
        mu=f/n
        #println(mu)
        if (n>1)&&(f>4)
            stdev=sqrt((f*(1-mu)*(1-mu) + (n-f)*(mu)*(mu))/(n-1))
            stderr=stdev/sqrt(n)
            if stderr<(err*mu)
                breakflag=1
            end
        end
        #if pct>5
        #    breakflag=1
        #end
        if (pct>9)
        atime=totime/(pct-1)
            totime=0
        pct=0
            println("mu is")
            println(mu)
            #println(nfail/ntr)
            println("number of trials is")
            println(n)
            println("number of Z or Y failures is")
            println(f)
            println("number of X failures is")
            println(fx)
        println("average trial time")
        println(atime)
            if f>4
                println("stderr is")
                println(stderr)
                println("target is")
                println(err*mu)
            end
        end
        fout = open(fname,"w")
        writedlm(fout,[n,f])
        close(fout)
        end
        println(f)
        println(fx)
        println(n)
        println((f+fx)/n)

    return
end

dzin=parse(Int64,ARGS[1])
pin=parse(Float64,ARGS[2])
#println(pin)
ntin=parse(Int64,ARGS[3])
cutin=parse(Float64,ARGS[4])
bdin=parse(Int64,ARGS[5])
if length(ARGS)==6
    sidin=parse(Int64,ARGS[6])
    SurfMC(dzin,3,dzin,pin,8,500e-9,1e7,0,cutin,10,0.1,ntin,sim_id = sidin)
else
    SurfMC(dzin,3,dzin,pin,8,500e-9,1e7,0,cutin,10,0.1,ntin)
end
