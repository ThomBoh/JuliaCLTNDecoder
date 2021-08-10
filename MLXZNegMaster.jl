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

function gate(::GateName"idle"; pz::Number, px::Number, py::Number)
    return [1-pz-px-py px          pz          py
            px         1-px-py-pz  py          pz
            pz         py          1-px-py-pz  px
            py         pz          px          1-px-py-pz]
end

function gate(::GateName"Ystab"; pz::Number)
    return [1   0   0   0
            0   1   0   0
            0   0   0   1
            0   0   1   0]
end

#function gate(::GateName"pz"; p::Number)
#  return [1-p  p
#          p    1-p]
#end

# CNOT gate failure -- first qubit is control, second qubit is target


function bin2dec(bin)

    l=length(bin)
    out=0
    for i in [1:l;]

        out=out+(bin[i]%2)*2^(i-1)

    end

    return out

end
#(::GateName"CX"; px2::Number,pz2::Number,py2::Number,px1::Number,px1x2::Number,px1z2::Number,px1y2::Number,pz1::Number,pz1x2::Number,pz1z2::Number,pz1y2::Number,py1::Number,py1x2::Number,py1z2::Number,py1y2::Number)

function gate(::GateName"CXT";px2::Number,pz2::Number,py2::Number,px1::Number,px1x2::Number,px1z2::Number,px1y2::Number,pz1::Number,pz1x2::Number,pz1z2::Number,pz1y2::Number,py1::Number,py1x2::Number,py1z2::Number,py1y2::Number)
    cx=buildCX(px2,pz2,py2,px1,px1x2,px1z2,px1y2,pz1,pz1x2,pz1z2,pz1y2,py1,py1x2,py1z2,py1y2)
    #println(typeof(cx))
    #println(size(cx))
    #println(cx)
    return cx
end

function gate(::GateName"CYT";px2::Number,pz2::Number,py2::Number,px1::Number,px1x2::Number,px1z2::Number,px1y2::Number,pz1::Number,pz1x2::Number,pz1z2::Number,pz1y2::Number,py1::Number,py1x2::Number,py1z2::Number,py1y2::Number)
    cy=buildCY(px2,pz2,py2,px1,px1x2,px1z2,px1y2,pz1,pz1x2,pz1z2,pz1y2,py1,py1x2,py1z2,py1y2)
    #println(typeof(cy))
    #println(size(cy))
    #println(cy)
    return cy
end

function buildCX(px2::Number,pz2::Number,py2::Number,px1::Number,px1x2::Number,px1z2::Number,px1y2::Number,pz1::Number,pz1x2::Number,pz1z2::Number,pz1y2::Number,py1::Number,py1x2::Number,py1z2::Number,py1y2::Number)

    #al2=8
    #p=1e-4

    #@show px2=0.93*((p)^(0.5))*exp(-2*al2)
    #@show pz2=0.15*((p)^(0.5))
    #@show py2=0.28*p*exp(-2*al2)
    #@show px1=0.93*((p)^(0.5))*exp(-2*al2)
    #@show px1x2=0.93*((p)^(0.5))*exp(-2*al2)
    #@show px1z2=0.28*p*exp(-2*al2)
    #@show px1y2=0.28*p*exp(-2*al2)
    #@show pz1=0.91*((p)^(0.5))
    #@show pz1x2=0.93*((p)^(0.5))*exp(-2*al2)
    #@show pz1z2=0.15*((p)^(0.5))
    #@show pz1y2=0.28*p*exp(-2*al2)
    #@show py1=0.93*((p)^(0.5))*exp(-2*al2)
    #@show py1x2=0.93*((p)^(0.5))*exp(-2*al2)
    #@show py1z2=0.28*p*exp(-2*al2)
    #@show py1y2=0.28*p*exp(-2*al2)


    probs = [px2,pz2,py2,px1,px1x2,px1z2,px1y2,pz1,pz1x2,pz1z2,pz1y2,py1,py1x2,py1z2,py1y2]
    tot=sum(probs)
    dist=[1-tot,px2,pz2,py2,px1,px1x2,px1z2,px1y2,pz1,pz1x2,pz1z2,pz1y2,py1,py1x2,py1z2,py1y2]
    perm=[1,2,11,12,6,5,16,15,9,10,3,4,14,13,8,7]
    gateout=zeros(Float64,16,16)
    #println(typeof(gateout))
    for i in 1:16

        for j in 1:16
            j2=perm[j]
            bini=digits(i-1,base=2,pad=4)
            binj2=digits(j2-1,base=2,pad=4)
            binerrno=bini+binj2
            errno=bin2dec(binerrno)+1
            gateout[i,j]=dist[errno]

        end
    end
    #println("CX")
    #display(gateout)
    #println(gateout)
    #println(size(gateout))
    return gateout
end

function buildCY(px2::Number,pz2::Number,py2::Number,px1::Number,px1x2::Number,px1z2::Number,px1y2::Number,pz1::Number,pz1x2::Number,pz1z2::Number,pz1y2::Number,py1::Number,py1x2::Number,py1z2::Number,py1y2::Number)

    #al2=8
    #p=1e-4
    #@show CYpx2=0.35*exp(-2*al2)*(p^(0.5))
    #@show CYpz2=(0.32 + 1.05/al2)*(p^(0.5))
    #@show CYpy2=1.67*exp(-2*al2)*(p^(0.5))
    #@show CYpx1=1.17*exp(-2*al2)*(p^(0.5))
    #@show CYpx1x2=1.13*exp(-2*al2)*(p^(0.5))
    #@show CYpx1z2=1.32*exp(-2*al2)*(p^(0.5))
    #@show CYpx1y2=exp(-2*al2)*(p^(0.5))
    #@show CYpz1=(1.01 + 0.77/al2)*(p^(0.5))
    #@show CYpz1x2=0.42*exp(-2*al2)*(p^(0.5))
    #@show CYpz1z2=(0.32 + 1.05/al2)*(p^(0.5))
    #@show CYpz1y2=1.06*exp(-2*al2)*(p^(0.5))
    #@show CYpy1=0.86*exp(-2*al2)*(p^(0.5))
    #@show CYpy1x2=1.14*exp(-2*al2)*(p^(0.5))
    #@show CYpy1z2=1.35*exp(-2*al2)*(p^(0.5))
    #@show CYpy1y2=0.76*exp(-2*al2)*(p^(0.5))


    # +0.42*exp(-2*al2)*p +0.86*exp(-2*al2)*p +1.14*exp(-2*al2)*p
    #@show CYpz2=(0.32 + 1.05/al2)*(p^(0.5))# +1.67*exp(-2*al2)*p +1.32*exp(-2*al2)*p +exp(-2*al2)*p
    # +1.35*exp(-2*al2)*p +1.06*exp(-2*al2)*p +0.76*exp(-2*al2)*p

    probs = [px2,pz2,py2,px1,px1x2,px1z2,px1y2,pz1,pz1x2,pz1z2,pz1y2,py1,py1x2,py1z2,py1y2]
    tot=sum(probs)
    dist=[1-tot,px2,pz2,py2,px1,px1x2,px1z2,px1y2,pz1,pz1x2,pz1z2,pz1y2,py1,py1x2,py1z2,py1y2]
    perm=[1,10,11,4,8,15,14,5,9,2,3,12,16,7,6,13]
    gateout=zeros(Float64,16,16)
    #println(typeof(gateout))
    for i in 1:16

        for j in 1:16

            j2=perm[j]
            bini=digits(i-1,base=2,pad=4)
            binj2=digits(j2-1,base=2,pad=4)
            binerrno=bini+binj2
            errno=bin2dec(binerrno)+1
            gateout[i,j]=dist[errno]

        end
    end
    #println("CY")
    #display(gateout)
    #println(gateout)
    #println(size(gateout))
    return gateout
end


#  return [1-pz1-pz2-pz12  pz12            pz1             pz2
#          pz2             pz1             pz12            1-pz1-pz2-pz12
#          pz1             pz2             1-pz1-pz2-pz12  pz12
#          pz12            1-pz1-pz2-pz12  pz2             pz1]
#end

function gate(::GateName"pz"; p::Number)
  return [ 1-p p
           p   1-p]
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

function gate(::GateName"Xres"; p::Number)
  return [ 1   1
           0   0]
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

function schedulemaker(dy::Integer,dx::Integer)

    #figure out how many y checks per row

    ycounts=zeros(Int,dy+1)
    for i in [1:(dy+1);]
        if i%2==1
            ycounts[i]=floor((dx-1)/2)
        else
            ycounts[i]=floor((dx)/2)
        end
    end

    #build y ancilla qubit schedule

    ysched=[]

    #i indexes which row of checks
    for i in [1:size(ycounts)[1];]
        #j indexes which check in row i
        for j in [1:ycounts[i];]

            if i==1

                sch=[]

                push!(sch,[1,2*(j-1)+3])
                push!(sch,[1,2*(j-1)+2])
                push!(sch,[-1,-1])
                push!(sch,[-1,-1])

                push!(ysched,sch)

            elseif i%2==0

                sch=[]

                if i==(size(ycounts)[1])

                    push!(sch,[-1,-1])
                    push!(sch,[-1,-1])
                    push!(sch,[i-1,2*(j-1)+2])
                    push!(sch,[i-1,2*(j-1)+1])

                else

                    push!(sch,[i,2*(j-1)+2])
                    push!(sch,[i,2*(j-1)+1])
                    push!(sch,[i-1,2*(j-1)+2])
                    push!(sch,[i-1,2*(j-1)+1])

                end

                push!(ysched,sch)

            elseif i%2==1

                sch=[]

                if i==(size(ycounts)[1])

                    push!(sch,[-1,-1])
                    push!(sch,[-1,-1])
                    push!(sch,[i-1,2*(j-1)+3])
                    push!(sch,[i-1,2*(j-1)+2])

                else

                    push!(sch,[i,2*(j-1)+3])
                    push!(sch,[i,2*(j-1)+2])
                    push!(sch,[i-1,2*(j-1)+3])
                    push!(sch,[i-1,2*(j-1)+2])

                end

                push!(ysched,sch)
            end
        end
    end

    #figure out how many x checks per column
    xcounts=zeros(Int,dx+1)
    for i in [1:(dx+1);]
        if i%2==1
            xcounts[i]=floor(dy/2)
        else
            xcounts[i]=floor((dy-1)/2)
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
                push!(sch,[2*(j-1)+2,1])
                push!(sch,[2*(j-1)+1,1])
                push!(sch,[-1,-1])
                push!(sch,[-1,-1])

                push!(xsched,sch)

            elseif i%2==0

                sch=[]

                if i==(size(xcounts)[1])

                    push!(sch,[-1,-1])
                    push!(sch,[-1,-1])
                    push!(sch,[2*(j-1)+3,i-1])
                    push!(sch,[2*(j-1)+2,i-1])

                else

                    push!(sch,[2*(j-1)+3,i])
                    push!(sch,[2*(j-1)+2,i])
                    push!(sch,[2*(j-1)+3,i-1])
                    push!(sch,[2*(j-1)+2,i-1])

                end

                push!(xsched,sch)

            elseif i%2==1

                sch=[]

                if i==(size(xcounts)[1])

                    push!(sch,[-1,-1])
                    push!(sch,[-1,-1])
                    push!(sch,[2*(j-1)+2,i-1])
                    push!(sch,[2*(j-1)+1,i-1])

                else

                    push!(sch,[2*(j-1)+2,i])
                    push!(sch,[2*(j-1)+1,i])
                    push!(sch,[2*(j-1)+2,i-1])
                    push!(sch,[2*(j-1)+1,i-1])

                end

                push!(xsched,sch)
            end
        end
    end

    println(ysched)
    println(xsched)
    return ysched, xsched
end

#create schedule of when data qubits on the boundary have wait locations instead of CNOTs

function bdrysched(dy::Integer,dx::Integer)

    bsch=zeros(Int,4,dy,dx)

    for iy in [1:dy;]

        for ix in [1:dx;]

            if (iy==1) && (ix==1)

                bsch[1,iy,ix]=1
                bsch[3,iy,ix]=1

            elseif iy==(dy) && (ix==1)

                if dy%2==0
                    bsch[3,iy,ix]=1
                    bsch[4,iy,ix]=1
                else
                    bsch[2,iy,ix]=1
                    bsch[1,iy,ix]=1
                end

            elseif (ix==(dx)) && (iy==1)

                if dx%2==0

                    bsch[1,iy,ix]=1
                    bsch[2,iy,ix]=1

                else

                    bsch[3,iy,ix]=1
                    bsch[4,iy,ix]=1
                end

            elseif (ix==dx) && (iy==dy)

                if (dx+dy)%2==0
                    bsch[2,iy,ix]=1
                    bsch[4,iy,ix]=1
                else
                    bsch[2,iy,ix]=1
                    bsch[4,iy,ix]=1
                end

            elseif (iy==1) && ((ix>1) && (ix<dx))

                if ix%2==1
                    bsch[3,iy,ix]=1
                else
                    bsch[1,iy,ix]=1
                end

            elseif (ix==1) && ((iy>1) && (iy<dy))

                if iy%2==1
                    bsch[1,iy,ix]=1
                else
                    bsch[3,iy,ix]=1
                end


            elseif (iy==dy) && ((ix>1) && (ix<dx))

                if dy%2==0

                    if ix%2==1
                        bsch[4,iy,ix]=1
                    else
                        bsch[2,iy,ix]=1
                    end
                else
                    if ix%2==1
                        bsch[2,iy,ix]=1
                    else
                        bsch[4,iy,ix]=1
                    end
                end

            elseif (ix==dx) && ((iy>1) && (iy<dy))

                if dx%2==0

                    if iy%2==1
                        bsch[2,iy,ix]=1
                    else
                        bsch[4,iy,ix]=1
                    end
                else
                    if iy%2==1
                        bsch[4,iy,ix]=1
                    else
                        bsch[2,iy,ix]=1
                    end
                end
            end
        end
    end

    return bsch
end

#reorder the x ancilla qubits for easier 1D layout

function xreorder(xsch,dy,dx)
    #assume dy,dx both odd
    #then there are dx+1 "columns" of (dz-1)/2 x checks each
    #and dz-1 "rows" of (dx+1)/2 x checks each
    nxc=size(xsch)[1]
    cpc=(dy-1)/2
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

function SurfLayout(dy::Integer,dx::Integer,ysch,xsch)


    nzc=size(ysch)[1]
    nxc=size(xsch)[1]
    ycpr=(dx-1)/2
    xcpr=(dy+1)/2
    qtot=2*dy*dx-1


    count=0
    mode=0
    yc=1
    xc=1
    drow=1
    layout=[]
    while count<(2*dy*dx-1)

        if mode==0
            #first/last row of y checks, type 1
            for i in [1:ycpr;]

                push!(layout,[1,[yc,yc]])
                yc=yc+1
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
            if drow>dy
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
                    push!(layout,[1,[yc,yc]])
                    yc=yc+1
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

function linemaps(layout,dy,dx,nyc,nxc)

    qlinemap=zeros(Int,dy,dx)
    ylinemap=zeros(Int,nyc)
    xlinemap=zeros(Int,nxc)
    N=2*dx*dy-1
    for i in [1:N;]

        typ=layout[i][1]

        if typ==0
            #data qubit
            addr=layout[i][2]
            yad=addr[1]
            xad=addr[2]

            qlinemap[yad,xad]=i

        elseif typ==1
            #y check qubit
            addr=layout[i][2][1]

            ylinemap[addr]=i

        elseif typ==2
            #x check qubit
            addr=layout[i][2][1]

            xlinemap[addr]=i

        end

    end

    return qlinemap,ylinemap,xlinemap
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


function CYfail(r,p,al2,tmeas,k2,nth,pmz,pmx)

    #pz1=0.91*((p)^(0.5))
    #pz2=0.15*((p)^(0.5))
    #pz1z2=0.15*((p)^(0.5))
    #px1=0.93*((p)^(0.5))*exp(-2*al2)
    #px2=0.93*((p)^(0.5))*exp(-2*al2)
    #px1x2=0.93*((p)^(0.5))*exp(-2*al2)
    #pz1x2=0.93*((p)^(0.5))*exp(-2*al2)
    #py1=0.93*((p)^(0.5))*exp(-2*al2)
    #py1x2=0.93*((p)^(0.5))*exp(-2*al2)
    #py2=0.28*p*exp(-2*al2)
    #py1z2=0.28*p*exp(-2*al2)
    #px1z2=0.28*p*exp(-2*al2)
    #pz1y2=0.28*p*exp(-2*al2)
    #py1y2=0.28*p*exp(-2*al2)
    #px1y2=0.28*p*exp(-2*al2)


    #CYpz1=(1.01 + 0.77/al2)*p
    #CYpz2=(0.32 + 1.05/al2)*p
    #CYpz1z2=(0.32 + 1.05/al2)*p
    CYpz1=(1.01 + 0.77/al2)*(p^(0.5))# +0.42*exp(-2*al2)*p +0.86*exp(-2*al2)*p +1.14*exp(-2*al2)*p
    CYpz2=(0.32 + 1.05/al2)*(p^(0.5))# +1.67*exp(-2*al2)*p +1.32*exp(-2*al2)*p +exp(-2*al2)*p
    CYpz1z2=(0.32 + 1.05/al2)*(p^(0.5))# +1.35*exp(-2*al2)*p +1.06*exp(-2*al2)*p +0.76*exp(-2*al2)*p

    CYpx1=1.17*exp(-2*al2)*(p^(0.5))
    CYpx2=0.35*exp(-2*al2)*(p^(0.5))
    CYpx1x2=1.13*exp(-2*al2)*(p^(0.5))
    CYpz1x2=0.42*exp(-2*al2)*(p^(0.5))
    CYpy1=0.86*exp(-2*al2)*(p^(0.5))
    CYpy1x2=1.14*exp(-2*al2)*(p^(0.5))
    CYpy2=1.67*exp(-2*al2)*(p^(0.5))
    CYpy1z2=1.35*exp(-2*al2)*(p^(0.5))
    CYpx1z2=1.32*exp(-2*al2)*(p^(0.5))
    CYpz1y2=1.06*exp(-2*al2)*(p^(0.5))
    CYpy1y2=0.76*exp(-2*al2)*(p^(0.5))
    CYpx1y2=exp(-2*al2)*(p^(0.5))


    if r < CYpz1
        zout=[1,0]
        xout=[0,0]
    elseif (r>CYpz1) && (r<(CYpz1+CYpz2))
        zout=[0,1]
        xout=[0,0]
    elseif (r>(CYpz1+CYpz2)) && (r<(CYpz1+CYpz2+CYpz1z2))
        zout=[1,1]
        xout=[0,0]
    elseif (r>(CYpz1+CYpz2+CYpz1z2)) && (r<(CYpz1+CYpz2+CYpz1z2+CYpx1))
        zout=[0,0]
        xout=[1,0]
    elseif (r>(CYpz1+CYpz2+CYpz1z2+CYpx1)) && (r<(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2))
        zout=[0,0]
        xout=[0,1]
    elseif (r>(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2)) && (r<(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2))
        zout=[0,0]
        xout=[1,1]
    elseif (r>(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2)) && (r<(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2))
        zout=[1,0]
        xout=[0,1]
    elseif (r>(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2)) && (r<(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2+CYpy1))
        zout=[1,0]
        xout=[1,0]
    elseif (r>(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2+CYpy1)) && (r<(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2+CYpy1+CYpy1x2))
        zout=[1,0]
        xout=[1,1]
    elseif (r>(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2+CYpy1+CYpy1x2)) && (r<(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2+CYpy1+CYpy1x2+CYpy2))
        zout=[0,1]
        xout=[0,1]
    elseif (r>(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2+CYpy1+CYpy1x2+CYpy2)) && (r<(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2+CYpy1+CYpy1x2+CYpy2+CYpy1z2))
        zout=[1,1]
        xout=[1,0]
    elseif (r>(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2+CYpy1+CYpy1x2+CYpy2+CYpy1z2)) && (r<(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2+CYpy1+CYpy1x2+CYpy2+CYpy1z2+CYpx1z2))
        zout=[0,1]
        xout=[1,0]
    elseif (r>(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2+CYpy1+CYpy1x2+CYpy2+CYpy1z2+CYpx1z2)) && (r<(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2+CYpy1+CYpy1x2+CYpy2+CYpy1z2+CYpx1z2+CYpz1y2))
        zout=[1,1]
        xout=[0,1]
    elseif (r>(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2+CYpy1+CYpy1x2+CYpy2+CYpy1z2+CYpx1z2+CYpz1y2)) && (r<(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2+CYpy1+CYpy1x2+CYpy2+CYpy1z2+CYpx1z2+CYpz1y2+CYpy1y2))
        zout=[1,1]
        xout=[1,1]
    elseif (r>(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2+CYpy1+CYpy1x2+CYpy2+CYpy1z2+CYpx1z2+CYpz1y2+CYpy1y2)) && (r<(CYpz1+CYpz2+CYpz1z2+CYpx1+CYpx2+CYpx1x2+CYpz1x2+CYpy1+CYpy1x2+CYpy2+CYpy1z2+CYpx1z2+CYpz1y2+CYpy1y2+CYpx1y2))
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


function SurfMCErrorOld(dz,dx,nr,zsched,xsched,bsch,p,al2,tmeas,k2,nth,pmz,pmx)

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

function SurfMCError(dz,dx,nr,zsched,xsched,bsch,p,al2,tmeas,k2,nth,pmz,pmx)

    ppaxz=(15/2)*p
    ppaxy=0#(15/2)*p

    ppazx=0.39*exp(-4*al2)
    ppazy=0#.39*exp(-4*al2)

    pzip=10*p#*(1+exp(-4*al2))
    pxip=0
    pyip=10*p*exp(-4*al2)

    pzic=0.31*((p)^(0.5))#*(1+exp(-4*al2))
    pxic=0
    pyic=0.31*((p)^(0.5))*exp(-4*al2)

    pzim=k2*p*al2*(tmeas+(10/(k2*al2)) )*(1+2*nth)#*(1+exp(-4*al2))
    pxim=0
    pyim=k2*p*al2*(tmeas+(10/(k2*al2)) )*exp(-4*al2)

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
        #println("ROUND 1 FAULTS")
        #t=0

        for zc in [1:nzc;]
            r=rand(Float64)
            if r<ppazx
                #Z check ancilla gets X error
                zancx[zc]=(zancx[zc]+1)%2
            end
        end
        #println("Init X Ancillae")
        for xc in [1:nxc;]
            r=rand(Float64)
            if r<ppaxz
                #println("X ancilla $(xc) gets Z fault")
                #x check ancilla gets z error
                xancz[xc]=(xancz[xc]+1)%2
            end
        end
        #println("Data Init Idle")
        for zq in [1:dz;]
            for xq in [1:dx;]
                if rep==1
                    r=rand(Float64)
                    qub=[zq,xq]
                    if r<pzip
                        #println("Data qubit $(qub) gets Z fault")
                        #data qubit gets z error
                        dataz[zq,xq]=(dataz[zq,xq]+1)%2
                    #elseif (r>pzip) && (r<(pzip+pxip))
                        #data qubit gets x error
                    #    datax[zq,xq]=(datax[zq,xq]+1)%2
                    elseif (r>(pzip)) && (r<(pzip+pxip))
                        #println("Data qubit $(qub) gets X fault")
                        #data qubit gets x error
                        datax[zq,xq]=(datax[zq,xq]+1)%2
                    elseif (r>(pzip+pxip)) && (r<pzip+pxip+pyip)
                        #println("Data qubit $(qub) gets Y fault")
                        #data qubit gets y error
                        dataz[zq,xq]=(dataz[zq,xq]+1)%2
                        datax[zq,xq]=(datax[zq,xq]+1)%2
                    end
                end
            end
        end

        #t=1-4, CNOTs and CYs

        for cn in [1:4;]
            #println("CXY ROUND $(cn)")
            for zc in [1:nzc;]
                #z-type check is to identify x errors on data
                #data is CNOT control, ancilla is target
                qub=zsched[zc][cn]
                if qub[1]>-1

                    r=rand(Float64)
                    zf,xf=CNOTfail(r,p,al2,tmeas,k2,nth,pmz,pmx)
                    ei=zf[1]+zf[2]+xf[1]+xf[2]
                    if ei>0
                        #println("CY FAILURE between ancilla $(yc) and qubit $(qub)")
                        #println(zf)
                        #println(xf)
                    end

		            datazold=dataz[qub[1],qub[2]]
		            dataxold=datax[qub[1],qub[2]]
		            zancxold=zancx[zc]
                    zanczold=zancz[zc]

                    #z error of the data qubit (control) inherits prior z error from ancilla (target)
                    #as well as the control z error from the CNOT operation
                    dataz[qub[1],qub[2]]=(zanczold+datazold+zf[1])%2
                    #z error of the ancilla qubit (target) inherits target z error from CNOT operation
                    zancz[zc]=(zanczold+zf[2])%2


                    #x error of the ancilla qubit (target) inherits prior x error from data (control)
                    #as well as the target x error from the CNOT operation
                    zancx[zc]=(dataxold+zancxold+xf[2])%2
                    #x error of the data qubit (control) inherits control x error from CNOT operation
                    datax[qub[1],qub[2]]=(dataxold+xf[1])%2

                else

                    #ancilla wait
                    r=rand(Float64)

                    if r<pzic
                        #println("Y ancilla $(yc) idle Z FAULT")
                        zancz[zc]=(zancz[zc]+1)%2
                    #elseif (r>pzic) && (r<(pzic+pxic))
                    #    zancx[zc]=(zancx[zc]+1)%2
                    elseif (r>(pzic)) && (r<(pzic+pyic))
                        #println("Y ancilla $(yc) idle Y FAULT")
                        zancz[zc]=(zancz[zc]+1)%2
                        zancx[zc]=(zancx[zc]+1)%2
                    elseif (r>(pzic+pyic)) && (r<(pzic+pyic+pxic))
                        #println("Y ancilla $(yc) idle X FAULT")
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
                    ei=zf[1]+zf[2]+xf[1]+xf[2]
                    #if ei>0
                        #println("CNOT FAILURE between ancilla $(xc) and qubit $(qub)")
                        #println(zf)
                        #println(xf)
                    #end
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
                    #println(dataz)
                    #println(datax)
                    #println(xc)
                    #println(xancz)
                else

                    #ancilla wait

                    r=rand(Float64)

                    if r<pzic
                        #println("X ancilla $(xc) idle Z FAULT")
                        xancz[xc]=(xancz[xc]+1)%2
                    #elseif (r>pzic) && (r<(pzic+pxic))
                    #    zancx[zc]=(zancx[zc]+1)%2
                    elseif (r>(pzic)) && (r<(pzic+pyic))
                        #println("X ancilla $(xc) idle Y FAULT")
                        xancz[xc]=(xancz[xc]+1)%2
                        xancx[xc]=(xancx[xc]+1)%2
                    elseif (r>(pzic+pyic)) && (r<(pzic+pyic+pxic))
                        #println("X ancilla $(xc) idle X FAULT")
                        xancx[xc]=(xancx[xc]+1)%2
                    end
                end
            end

            #boundary data qubit wait locations

            for zq in [1:dz;]
                for xq in [1:dx;]
                    qub=[zq,xq]
                    if bsch[cn,zq,xq]==1
                        r=rand(Float64)
                        if r<pzic
                            #println("DATA qubit $(qub) idle Z FAULT")
                            dataz[zq,xq]=(dataz[zq,xq]+1)%2
                        #elseif (r>pzic) && (r<(pzic+pxic))
                        #    datax[zq,xq]=(datax[zq,xq]+1)%2
                        elseif (r>(pzic)) && (r<(pzic+pyic))
                            #println("DATA qubit $(qub) idle Y FAULT")
                            dataz[zq,xq]=(dataz[zq,xq]+1)%2
                            datax[zq,xq]=(datax[zq,xq]+1)%2
                        elseif (r>(pzic+pyic)) && (r<(pzic+pyic+pxic))
                            #println("DATA qubit $(qub) idle X FAULT")
                            datax[zq,xq]=(datax[zq,xq]+1)%2
                        end
                    end
                end
            end
        end

        #t=5

        #final data qubit waits
        #println("MEASUREMENT TIME")
        for zq in [1:dz;]
            for xq in [1:dx;]
                qub=[zq,xq]
                r=rand(Float64)
                if r<pzim
                    #println("DATA qubit $(qub) idle Z FAULT")
                    dataz[zq,xq]=(dataz[zq,xq]+1)%2
                elseif (r>(pzim)) && (r<(pzim+pyim))
                    #println("DATA qubit $(qub) idle Y FAULT")
                    dataz[zq,xq]=(dataz[zq,xq]+1)%2
                    datax[zq,xq]=(datax[zq,xq]+1)%2
                elseif (r>(pzim+pyim)) && (r<(pzim+pyim+pxim))
                    #println("DATA qubit $(qub) idle X FAULT")
                    datax[zq,xq]=(datax[zq,xq]+1)%2
                end
            end
        end

        #final ancilla qubit measurement failures
        for zc in [1:nzc;]

            r=rand(Float64)
            if r<pmz
                #println("Y ancilla $(yc) MEASURMENT FAULT")
                zancz[zc]=(zancz[zc]+1)%2
            end
        end

        for xc in [1:nxc;]

            r=rand(Float64)
            if r<pmx
                #println("X ancilla $(xc) MEASURMENT FAULT")
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



#the following function takes error configurations on data and ancilla qubits
#and returns an array of strings that can be used to create the appropriate
#MPS to overlap with the circuit tensor network in order to retrieve the desired
#element of the MPS

function buildstrings(dz::Int,dx::Int,nr,pez,pex,Synz,Synx,layout)

    N=2*dz*dx-1
    outi=[]
    outx=[]
    outy=[]
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

            if qx==dx && qz>1
                #on Z logical boundary
                #println(pez)
                #println(qz)
                #println(qx)
                ez=pez[qz,qx]%2
                ex=pex[qz,qx]%2

                if ez==0
                    #is there a Z component presenton Z logical boundary
                    push!(outi,"Z+")
                    push!(outz,"Z-")
                    push!(outy,"Z-")
                    push!(outx,"Z+")
                elseif ez==1
                    push!(outi,"Z-")
                    push!(outz,"Z+")
                    push!(outy,"Z+")
                    push!(outx,"Z-")
                else
                    println("bad news!")
                end
                if ex==0
                    #is there an X component presenton Z logical boundary
                    push!(outi,"Z+")
                    push!(outz,"Z+")
                    push!(outy,"Z+")
                    push!(outx,"Z+")
                elseif ex==1
                    push!(outi,"Z-")
                    push!(outz,"Z-")
                    push!(outy,"Z-")
                    push!(outx,"Z-")
                else
                    println("bad news!")
                end
                count=count+1
            elseif qx==dx && qz==1
                #on X and Z logical boundaries (corner)
                ez=pez[qz,qx]%2
                ex=pex[qz,qx]%2
                #println(e)
                if ez==0
                    #is there a Z component present on Z+X logical boundary
                    push!(outi,"Z+")
                    push!(outz,"Z-")
                    push!(outy,"Z-")
                    push!(outx,"Z+")
                elseif ez==1
                    push!(outi,"Z-")
                    push!(outz,"Z+")
                    push!(outy,"Z+")
                    push!(outx,"Z-")
                else
                    println("bad news!")
                end
                if ex==0
                    #is there an X component present on Z+X logical boundary
                    push!(outi,"Z+")
                    push!(outz,"Z+")
                    push!(outy,"Z-")
                    push!(outx,"Z-")
                elseif ex==1
                    push!(outi,"Z-")
                    push!(outz,"Z-")
                    push!(outy,"Z+")
                    push!(outx,"Z+")
                else
                    println("bad news!")
                end
            elseif qz==1 && qx<dx
                #on X logical boundary
                ez=pez[qz,qx]%2
                ex=pex[qz,qx]%2
                #println(e)
                if ez==0
                    #is there a Z component present on X logical boundary
                    push!(outi,"Z+")
                    push!(outz,"Z+")
                    push!(outy,"Z+")
                    push!(outx,"Z+")
                elseif ez==1
                    push!(outi,"Z-")
                    push!(outz,"Z-")
                    push!(outy,"Z-")
                    push!(outx,"Z-")
                else
                    println("bad news!")
                end
                if ex==0
                    #is there an X component present on X logical boundary
                    push!(outi,"Z+")
                    push!(outz,"Z+")
                    push!(outy,"Z-")
                    push!(outx,"Z-")
                elseif ex==1
                    push!(outi,"Z-")
                    push!(outz,"Z-")
                    push!(outy,"Z+")
                    push!(outx,"Z+")
                else
                    println("bad news!")
                end
            else
                #println(pez)
                #println(qz)
                #println(qx)
                ez=pez[qz,qx]%2
                ex=pex[qz,qx]%2
                #println(e)
                #e=pez[qz,qx]
                if ez==0
                    push!(outi,"Z+")
                    push!(outz,"Z+")
                    push!(outy,"Z+")
                    push!(outx,"Z+")
                elseif ez==1
                    push!(outi,"Z-")
                    push!(outz,"Z-")
                    push!(outy,"Z-")
                    push!(outx,"Z-")
                else
                    println("bad news!")
                end
                if ex==0
                    push!(outi,"Z+")
                    push!(outz,"Z+")
                    push!(outy,"Z+")
                    push!(outx,"Z+")
                elseif ex==1
                    push!(outi,"Z-")
                    push!(outz,"Z-")
                    push!(outy,"Z-")
                    push!(outx,"Z-")
                else
                    println("bad news!")
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
                    outi=["Z+"]
                    outz=["Z+"]
                    outy=["Z+"]
                    outx=["Z+"]
                else
                    #println("hay3")
                    outi=["Z-"]
                    outz=["Z-"]
                    outy=["Z-"]
                    outx=["Z-"]
                end


                #push!(outi,"X+")
                #push!(outz,"X+")
                #push!(outy,"X+")
                #push!(outx,"X+")
                count=count+1
            else
                if Synz[nr,zc]==0
                    push!(outi,"Z+")
                    push!(outz,"Z+")
                    push!(outy,"Z+")
                    push!(outx,"Z+")
                else
                    push!(outi,"Z-")
                    push!(outz,"Z-")
                    push!(outy,"Z-")
                    push!(outx,"Z-")
                end
                #push!(outi,"X+")
                #push!(outz,"X+")
                #push!(outy,"X+")
                #push!(outx,"X+")
                count=count+1
            end

        elseif typ==2
            #x check ancilla
            xc=layout[i][2][1]
            if Synx[nr,xc]==0
                push!(outi,"Z+")
                push!(outz,"Z+")
                push!(outy,"Z+")
                push!(outx,"Z+")
            else
                push!(outi,"Z-")
                push!(outz,"Z-")
                push!(outy,"Z-")
                push!(outx,"Z-")
            end
            #push!(outi,"X+")
            #push!(outz,"X+")
            #push!(outy,"X+")
            #push!(outx,"X+")
            count=count+1
        end

    end
    #println("IZYX")
    #println(outi)
    #println(outz)
    #println(outy)
    #println(outx)
    return outi,outz,outy,outx

end


function buildinitstring(dz::Int,dx::Int,nr,layout)

    N=2*dz*dx-1
    out=[]

    count=0


    for i in [1:N;]
        #println(i)
        typ=layout[i][1]
#        println(layout[i])
        #println(typ)
        if typ==0
           #data qubit
            push!(out,"Z+")
            push!(out,"Z+")
            count=count+1

        elseif typ==1
            #println("hay1")
            #y check ancilla

            if i==1
                out=["Z+"]
                #push!(out,"X+")
                #println("hay2")
                count=count+1
            else
                push!(out,"Z+")
                #push!(out,"X+")
                count=count+1
            end

        elseif typ==2
            push!(out,"Z+")
            #push!(out,"X+")
            #x check ancilla
            count=count+1
        end

    end
    #println(outi)
    #println(outz)
    return out

end

function YStabs(dy,dx)

    stabgens=[]
    #stabgens = Array{Array{Int,(dy,dx)}}(dy+1)
    for i in [0:dy;]

         stab=zeros(Int,dy,dx)
         if i==0
             stab[1,2]=1
             stab[1,3]=1
         elseif i==dy
             stab[i,1]=1
             stab[i,2]=1
         elseif i%2==1
             stab[i,1]=1
             stab[i,2]=1
             stab[i+1,2]=1
             stab[i+1,1]=1
         elseif i%2==0
             stab[i,2]=1
             stab[i,3]=1
             stab[i+1,2]=1
             stab[i+1,3]=1
         end
         push!(stabgens,stab)
    end
    stabsout=[]
    #stabsout = Array{Array{Int,(dy,dx)}}(2^(dy+1))
    for i in [0:(2^(dy+1)-1);]
        bi=digits(i,base=2,pad=dy+1)
        stab=zeros(Int,dy,dx)
        for j in [1:(dy+1);]
            stab=stab+bi[j]*stabgens[j]
        end
        display(stab)
        println("\n")
        push!(stabsout,stab)
        #stabsout[i]=stab
    end
    return stabsout
end

function XStabs(dy,dx)

    stabgens=[]
    #stabgens = Array{Array{Int,(dy,dx)}}(dy+1)
    for i in [0:dy;]

         stab=zeros(Int,dy,dx)
         if i==0
             stab[1,1]=1
             stab[2,1]=1
         elseif i==dx
             stab[2,i]=1
             stab[3,i]=1
         elseif i%2==1
             stab[i,2]=1
             stab[i,3]=1
             stab[i+1,2]=1
             stab[i+1,3]=1
         elseif i%2==0
             stab[i,1]=1
             stab[i,2]=1
             stab[i+1,1]=1
             stab[i+1,2]=1
         end
         push!(stabgens,stab)
    end
    stabsout=[]
    #stabsout = Array{Array{Int,(dy,dx)}}(2^(dy+1))
    for i in [0:(2^(dy+1)-1);]
        bi=digits(i,base=2,pad=dy+1)
        stab=zeros(Int,dy,dx)
        for j in [1:(dy+1);]
            stab=stab+bi[j]*stabgens[j]
        end
        display(stab)
        println("\n")
        push!(stabsout,stab)
        #stabsout[i]=stab
    end
    return stabsout
end

function MLError(B1,B2,dz,dx,nr,PEZ,PEX,Synz,Synx,zsch,xsch,bsch,layout,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)
    #@show B
    phi=productstate(B1)
    instring=buildinitstring(dz,dx,nr,layout)
    psi=productstate(siteinds(phi),instring)
    #psi=productstate(B)
    #@show psi
    #println("noprime")
    #psi=yfilter(psi,dz,dx,zsch,xsch,layout,acc,bd)8
    A=noprime(*(B1,psi;cutoff=acc,maxdim=bd))
    #@show MPS
    #println(norm(A))
    for re in [1:(nr-1);]
        #println("glue")
        if nr==1
            println("WARNING YOU SHOULDNT SEE THIS!!!!")
        end
        A=glue(A,dz,dx,nr,re,PEZ,PEX,Synz,Synx,zsch,xsch,bsch,layout,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)
        #println(norm(A))
        #@show MPS
        #println("noprime")
        A=noprime(*(B2,A;cutoff=acc,maxdim=bd))
        #println(norm(A))
        #@show MPS
    end
    A=gluefinal(A,dz,dx,nr,nr,PEZ,PEX,Synz,Synx,zsch,xsch,bsch,layout,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)
    #println("glue")
    #psi=glue(psi,dz,dx,nr,nr,PEZ,PEX,Synz,Synx,zsch,xsch,bsch,layout,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)
    #@show MPS
    #display(PEZ)
    #display(PEX)
    #display(Synz)
    #display(Synx)
    stringi,stringz,stringy,stringx=buildstrings(dz,dx,nr,PEZ,PEX,Synz,Synx,layout)
    mpsi=productstate(siteinds(A),stringi)
    mpsz=productstate(siteinds(A),stringz)
    mpsy=productstate(siteinds(A),stringy)
    mpsx=productstate(siteinds(A),stringx)

    #compute "probability" that this pure error and meas outcomes occurred, and no logical z
    pli=inner(A,mpsi)

    #compute "probability" that this pure error and meas outcomes occurred, as well as logical z
    plz=inner(A,mpsz)
    ply=inner(A,mpsy)
    plx=inner(A,mpsx)
    #println(pli)
    #println(plz)
    #println(ply)
    #println(plx)
    #println(ml)

    return pli,plz,ply,plx

end

function glue(psi,dz,dx,nr,re,PEZ,PEX,Synz,Synx,zsch,xsch,bsch,layout,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)

    #First initialize all hardware noise parameters and failure rates. Note that p=k1/k2
    j=0
    ppax=(15/2)*p
    ppaz=0.39*exp(-4*al2)
    pzip=10*p#*(1+exp(-4*al2))
    pyip=10*p*exp(-4*al2)
    pzic=0.31*((p)^(0.5))#*(1+exp(-4*al2))
    pyic=0.31*((p)^(0.5))*exp(-4*al2)
    pzim=k2*p*al2*(tmeas)*(1+2*nth)#*(1+exp(-4*al2))
    pyim=k2*p*al2*(tmeas)*exp(-4*al2) #tmeas+(10/(k2*al2))
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


    nzc=size(zsch)[1]
    nxc=size(xsch)[1]

    nq=size(layout)[1]


    for q in [1:nq;]

        #first time step -- initialization/idling
        typ=layout[q][1]

        if typ==0
            #data qubit
            push!(gates,("idle",(2*q-1-j,2*q-j),(pz=0,px=0,py=0)))
            #push!(gates("pz",2*q-1-j,(p=0,)))
            #push!(gates("pz",2*q-j,(p=0,)))


        elseif typ==1
            #y check ancilla qubit
            zc=layout[q][2][1]
            bit=Synz[re,zc]
            if bit==0
                push!(gates,("pi0",2*q-1-j,(p=0,)))
                #push!(gates,("Xres",2*q,(p=0,)))
            else
                push!(gates,("pi1",2*q-1-j,(p=0,)))
                #push!(gates,("Xres",2*q,(p=0,)))
            end
            j=j+1

        elseif typ==2
            #x check ancilla qubit
            xc=layout[q][2][1]
            bit=Synx[re,xc]
            if bit==0
                push!(gates,("pi0",2*q-1-j,(p=0,)))
                #push!(gates,("Xres",2*q,(p=0,)))
            else
                push!(gates,("pi1",2*q-1-j,(p=0,)))
                #push!(gates,("Xres",2*q,(p=0,)))
            end
            j=j+1
        end

    end

    #istr=buildinitstring(dz,dx,nr,layout)
    #T=productstate(nq) #all zeros initial state sets initial error configuration to be trivial
    #TI=productstate(siteinds(T),istr)
    #display(T)
    psiOut=runcircuit(psi,gates,cutoff=acc,maxdim=bd)
    #display(psiOut)
    return psiOut #return the MPS, whose physical indices are as described at top of function
end


function gluefinal(psi,dz,dx,nr,re,PEZ,PEX,Synz,Synx,zsch,xsch,bsch,layout,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)

    #First initialize all hardware noise parameters and failure rates. Note that p=k1/k2
    j=0
    ppax=(15/2)*p
    ppaz=0.39*exp(-4*al2)
    pzip=10*p#*(1+exp(-4*al2))
    pyip=10*p*exp(-4*al2)
    pzic=0.31*((p)^(0.5))#*(1+exp(-4*al2))
    pyic=0.31*((p)^(0.5))*exp(-4*al2)
    pzim=k2*p*al2*(tmeas)*(1+2*nth)#*(1+exp(-4*al2))
    pyim=k2*p*al2*(tmeas)*exp(-4*al2) #tmeas+(10/(k2*al2))
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


    nzc=size(zsch)[1]
    nxc=size(xsch)[1]

    nq=size(layout)[1]


    for q in [1:nq;]

        #first time step -- initialization/idling
        typ=layout[q][1]

        if typ==0
            #data qubit
            push!(gates,("idle",(2*q-1-j,2*q-j),(pz=pzip,px=0,py=pyip)))
            #push!(gates("pz",2*q-1-j,(p=0,)))
            #push!(gates("pz",2*q-j,(p=0,)))


        elseif typ==1
            #y check ancilla qubit
            zc=layout[q][2][1]
            bit=Synz[re,zc]
            if bit==0
                push!(gates,("pz",2*q-1-j,(p=0,)))
                #push!(gates,("Xres",2*q,(p=0,)))
            else
                push!(gates,("pz",2*q-1-j,(p=0,)))
                #push!(gates,("Xres",2*q,(p=0,)))
            end
            j=j+1

        elseif typ==2
            #x check ancilla qubit
            xc=layout[q][2][1]
            bit=Synx[re,xc]
            if bit==0
                push!(gates,("pz",2*q-1-j,(p=0,)))
                #push!(gates,("Xres",2*q,(p=0,)))
            else
                push!(gates,("pz",2*q-1-j,(p=0,)))
                #push!(gates,("Xres",2*q,(p=0,)))
            end
            j=j+1
        end

    end

    #istr=buildinitstring(dz,dx,nr,layout)
    #T=productstate(nq) #all zeros initial state sets initial error configuration to be trivial
    #TI=productstate(siteinds(T),istr)
    #display(T)
    psiOut=runcircuit(psi,gates,cutoff=acc,maxdim=bd)
    #display(psiOut)
    return psiOut #return the MPS, whose physical indices are as described at top of function
end

function buildblock(b,dz,dx,nr,zsch,xsch,bsch,layout,ql,zl,xl,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)

    #First initialize all hardware noise parameters and failure rates. Note that p=k1/k2

    ppaxz=(15/2)*p
    ppaxy=0#(15/2)*p

    ppazx=0.39*exp(-4*al2)
    ppazy=0#.39*exp(-4*al2)

    if b==2

        pzip=10*p#*(1+exp(-4*al2))
        pxip=0
        pyip=10*p*exp(-4*al2)

    else
        pzip=0
        pxip=0
        pyip=0

    end

    pzic=0.31*((p)^(0.5))#*(1+exp(-4*al2))
    pxic=0
    pyic=0.31*((p)^(0.5))*exp(-4*al2)

    pzim=k2*p*al2*(tmeas)*(1+2*nth)#*(1+exp(-4*al2))#+(10/(k2*al2))
    pxim=0
    pyim=k2*p*al2*(tmeas )*exp(-4*al2)#+(10/(k2*al2))

    CXpz1=0.91*((p)^(0.5))#+3*0.93*((p)^(0.5))*exp(-2*al2)
    CXpz2=0.15*((p)^(0.5))#+3*0.28*p*exp(-2*al2)
    CXpz1z2=0.15*((p)^(0.5))#+3*0.28*p*exp(-2*al2)
    CXpx1=0.93*((p)^(0.5))*exp(-2*al2)
    CXpx2=0.93*((p)^(0.5))*exp(-2*al2)
    CXpx1x2=0.93*((p)^(0.5))*exp(-2*al2)
    CXpz1x2=0.93*((p)^(0.5))*exp(-2*al2)
    CXpy1=0.93*((p)^(0.5))*exp(-2*al2)
    CXpy1x2=0.93*((p)^(0.5))*exp(-2*al2)
    CXpy2=0.28*p*exp(-2*al2)
    CXpy1z2=0.28*p*exp(-2*al2)
    CXpx1z2=0.28*p*exp(-2*al2)
    CXpz1y2=0.28*p*exp(-2*al2)
    CXpy1y2=0.28*p*exp(-2*al2)
    CXpx1y2=0.28*p*exp(-2*al2)

    #cxm=buildCX(px2,pz2,py2,px1,px1x2,px1z2,px1y2,pz1,pz1x2,pz1z2,pz1y2,py1,py1x2,py1z2,py1y2)

    CYpz1=(1.01 + 0.77/al2)*(p^(0.5)) #+0.42*exp(-2*al2)*(p^(0.5)) +0.86*exp(-2*al2)*(p^(0.5)) +1.14*exp(-2*al2)*(p^(0.5))
    CYpz2=(0.32 + 1.05/al2)*(p^(0.5)) #+1.67*exp(-2*al2)*(p^(0.5)) +1.32*exp(-2*al2)*(p^(0.5)) +exp(-2*al2)*(p^(0.5))
    CYpz1z2=(0.32 + 1.05/al2)*(p^(0.5)) #+1.35*exp(-2*al2)*(p^(0.5)) +1.06*exp(-2*al2)*(p^(0.5)) +0.76*exp(-2*al2)*(p^(0.5))

    CYpx1=1.17*exp(-2*al2)*(p^(0.5))
    CYpx2=0.35*exp(-2*al2)*(p^(0.5))
    CYpx1x2=1.13*exp(-2*al2)*(p^(0.5))
    CYpz1x2=0.42*exp(-2*al2)*(p^(0.5))
    CYpy1=0.86*exp(-2*al2)*(p^(0.5))
    CYpy1x2=1.14*exp(-2*al2)*(p^(0.5))
    CYpy2=1.67*exp(-2*al2)*(p^(0.5))
    CYpy1z2=1.35*exp(-2*al2)*(p^(0.5))
    CYpx1z2=1.32*exp(-2*al2)*(p^(0.5))
    CYpz1y2=1.06*exp(-2*al2)*(p^(0.5))
    CYpy1y2=0.76*exp(-2*al2)*(p^(0.5))
    CYpx1y2=exp(-2*al2)*(p^(0.5))

    #cym=buildCY(CYpx2,CYpz2,CYpy2,CYpx1,CYpx1x2,CYpx1z2,CYpx1y2,CYpz1,CYpz1x2,CYpz1z2,CYpz1y2,CYpy1,CYpy1x2,CYpy1z2,CYpy1y2)

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
                push!(gates,("idle",(2*q-1,2*q),(pz=pzip,px=pxip,py=pyip)))

            elseif typ==1
                #z check ancilla qubit
                push!(gates,("idle",(2*q-1,2*q),(pz=0,px=ppazx,py=ppazy)))

            elseif typ==2
                #x check ancilla qubit
                push!(gates,("idle",(2*q-1,2*q),(pz=ppaxz,px=0,py=ppaxy)))

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

                        push!(gates,("idle",(2*q-1,2*q),(pz=pzic,px=pxic,py=pyic)))

                    end

                elseif typ==1

                    #Z check ancilla qubit, identifies X errors, data is target, ancilla is control

                    zc=layout[q][2][1] #which z check it is in the list of z check schedules
                    partner=zsch[zc][cnt]

                    if partner[1]==-1
                        #wait location for this ancilla
                        push!(gates,("idle",(2*q-1,2*q),(pz=pzic,px=pxic,py=pyic)))
                    else
                        pzad=partner[1]
                        pxad=partner[2]
                        q2=ql[pzad,pxad]
                        push!(gates,("CXT",(2*q2-1,2*q2,2*q-1,2*q),(px2=CXpx2,pz2=CXpz2,py2=CXpy2,px1=CXpx1,px1x2=CXpx1x2,px1z2=CXpx1z2,px1y2=CXpx1y2,pz1=CXpz1,pz1x2=CXpz1x2,pz1z2=CXpz1z2,pz1y2=CXpz1y2,py1=CXpy1,py1x2=CXpy1x2,py1z2=CXpy1z2,py1y2=CXpy1y2)))
                    end

                elseif typ==2

                    #x check ancilla qubit, identifies z errors, data is target, ancilla is control

                    xc=layout[q][2][1] #which x check it is in the list of x check schedules
                    partner=xsch[xc][cnt]

                    if partner[1]==-1
                        #wait location for this ancilla
                        push!(gates,("idle",(2*q-1,2*q),(pz=pzic,px=pxic,py=pyic)))
                    else
                        pzad=partner[1]
                        pxad=partner[2]
                        q2=ql[pzad,pxad]
                        push!(gates,("CXT",(2*q-1,2*q,2*q2-1,2*q2),(px2=CXpx2,pz2=CXpz2,py2=CXpy2,px1=CXpx1,px1x2=CXpx1x2,px1z2=CXpx1z2,px1y2=CXpx1y2,pz1=CXpz1,pz1x2=CXpz1x2,pz1z2=CXpz1z2,pz1y2=CXpz1y2,py1=CXpy1,py1x2=CXpy1x2,py1z2=CXpy1z2,py1y2=CXpy1y2)))
                    end

                end

            end
        end


        for q in [1:nq;]


            #final time step -- idling/measurement
            typ=layout[q][1]

            if typ==0
                #data qubit
                push!(gates,("idle",(2*q-1,2*q),(pz=pzim,px=pxim,py=pyim)))

            elseif typ==1
                #z check ancilla qubit
                push!(gates,("idle",(2*q-1,2*q),(pz=0,px=pmz,py=0)))

            elseif typ==2
                #x check ancilla qubit
                push!(gates,("idle",(2*q-1,2*q),(pz=pmx,px=0,py=0)))

            end

        end
    end

    TOut=runcircuit(gates,cutoff=acc,maxdim=bd; process =  true)

    return TOut #return the MPS, whose physical indices are as described at top of function
end

function ZXtoYX(dy,dx,EZzxin,EXzxin)

    EYyxout=zeros(Int,dy,dx)
    EXyxout=zeros(Int,dy,dx)

    for i in [1:dy;]
        for j in [1:dx;]

            Zzx=EZzxin[i,j]
            Xzx=EXzxin[i,j]

            if Zzx==1 && Xzx==1

                EYyxout[i,j]=1

            elseif Zzx==0 && Xzx==1

                EXyxout[i,j]=1

            elseif Zzx==1 && Xzx==0

                EYyxout[i,j]=1
                EXyxout[i,j]=1

            end

        end
    end

    return EYyxout,EXyxout

end


#returns the logical Y and X operators

function SurfLogs(dy,dx)

    Ly=zeros(Int,dy,dx)
    Lx=zeros(Int,dy,dx)

    if dx%2==0

        Lx[dy,:]=ones(Int,dx)

    else

        Lx[1,:]=ones(Int,dx)
    end

    Ly[:,dx]=ones(Int,dy)

    return Ly,Lx
end

#returns the pure error components of the accumulated errors, up to stabilizer

function getSurfPE(Ey,Ex,dy,dx)

    Ly,Lx=SurfLogs(dy,dx)
    LC=LogComp(Ey,Ex,dy,dx)

    PEY=LC[1]*Ly + Ey
    PEX=LC[2]*Lx + Ex
    #PEY=zeros(Int,dz,dx)
    for i in [1:dy;]
        for j in [1:dx;]

            PEY[i,j]=PEY[i,j]%2
            PEX[i,j]=PEX[i,j]%2
#            if PEZ[i,j]==1 && PEX[i,j]==1
#                PEY[i,j]=1
#            end
        end
    end

    return PEY,PEX
end

#returns the logical components of the accumulated errors

function LogComp(Ey,Ex,dy,dx)

    LC=[0,0]

    if dx%2==0
        println("shouldn't see this!!")
        LC[1]=Ey[dy,dx]
        LC[2]=Ex[dy,dx]

    else
        for i in [1:dx;]
            LC[1]=(LC[1]+Ey[1,i])%2
        end
        for i in [1:dy;]
            LC[2]=(LC[2]+Ex[i,dx])%2
        end
    end

    return LC
end


function traceblock2(U,layout)

    sites = firstsiteinds(U)

    M = ITensor[]
    j=1
    for ix in 1:length(U)
        if ix%2==1
            i=Int((ix+1)/2)
        else
            i=Int(ix/2)
        end
        typ=layout[i][1]
        if (typ ==2 ) && (ix%2==0)
            tracetensor = ITensor([1,1], sites[ix])
            initensor = ITensor([1,0], sites[ix])
            T = U[ix] * initensor
            T = T * tracetensor'
            M[ix-j] = M[ix-j] * T
            j=j+1
        elseif (typ==1) && (ix%2==1)
            tracetensor = ITensor([1,1], sites[ix])
            initensor = ITensor([1,0], sites[ix])
            T = U[ix] * initensor
            T = T * tracetensor'
            if ix==1
                push!(M,T)
                j=j+1
            else
                M[ix-j] = M[ix-j] * T
                j=j+1
            end
        elseif ix==2
            M[1]=M[1]*U[2]
        else
            push!(M, U[ix])
        end
    end
    out=MPO(M)

    #println(typeof(B))
    #println(typeof(M))
    #println(typeof(out))
    #M=0
    #GC.gc()
    return out

end

function traceblock1(U,layout)

    sites = firstsiteinds(U)

    M = ITensor[]
    j=1
    for ix in 1:length(U)
        if ix%2==1
            i=Int((ix+1)/2)
        else
            i=Int(ix/2)
        end
        typ=layout[i][1]
        if (typ == 2) && (ix%2==0)
            tracetensor = ITensor([1,1], sites[ix])
            T = U[ix] * tracetensor
            T = T * tracetensor'
            M[ix-j] = M[ix-j] * T
            j=j+1
        elseif (typ==1) && (ix%2==1)
            tracetensor = ITensor([1,1], sites[ix])
            initensor = ITensor([1,1], sites[ix])
            T = U[ix] * initensor
            T = T * tracetensor'
            if ix==1
                push!(M,T)
                j=j+1
            else
                M[ix-j] = M[ix-j] * T
                j=j+1
            end
        elseif ix==2
            M[1]=M[1]*U[2]
        else
            push!(M, U[ix])
        end
    end
    out=MPO(M)

    #println(typeof(B))
    #println(typeof(M))
    #println(typeof(out))
    #M=0
    #GC.gc()
    return out

end

function simpleblocks(dz,dx,nr,zsch,xsch,bsch,layout,ql,zl,xl,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)

    #B1=buildblock(1,dz,dx,nr,zsch,xsch,bsch,layout,ql,zl,xl,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)
    B2=buildblock(2,dz,dx,nr,zsch,xsch,bsch,layout,ql,zl,xl,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)
    C=traceblock1(B2,layout)
    D=traceblock2(B2,layout)
    println(C)
    println(D)
    return C,D

end

function MLEs(pli,plz,ply,plx)
    plin,plzn,plyn,plxn=pli,plz,ply,plx
    nffl=0
    if pli < 0
        plin=-1*pli
        nffl=1
    end
    if plz < 0
        plzn=-1*plz
        nffl=1
    end
    if ply < 0
        plyn=-1*ply
        nffl=1
    end
    if plx < 0
        plxn=-1*plx
        nffl=1
    end

    if (pli>plz) && ((pli>ply) && (pli>plx))
        #if it was more likely that no logical error occurred, decoder does nothing
        ml=[0,0]
    elseif (plz>pli) && ((plz>ply) && (plz>plx))
        #if it was more likely that logical z occurred, decoder applies logical z
        ml=[1,0]
    elseif (plx>pli) && ((plx>ply) && (plx>plz))
        #if it was more likely that logical x occurred, decoder applies logical x
        ml=[0,1]
    else
        #else it's a logical Y error
        ml=[1,1]
    end

    if (plin>plzn) && ((plin>plyn) && (plin>plxn))
        #if it was more likely that no logical error occurred, decoder does nothing
        mln=[0,0]
    elseif (plzn>plin) && ((plzn>plyn) && (plzn>plxn))
        #if it was more likely that logical z occurred, decoder applies logical z
        mln=[1,0]
    elseif (plxn>plin) && ((plxn>plyn) && (plxn>plzn))
        #if it was more likely that logical x occurred, decoder applies logical x
        mln=[0,1]
    else
        #else it's a logical Y error
        mln=[1,1]
    end

    return ml, mln, nffl

end

function SurfMC(dz,dx,nr,p,al2,tmeas,k2,nth,acc,bd,err,nt; sim_id::Int=-1)

    if sim_id < 0
      fname = "XZ_dy$(dz)_dx$(dx)_p$(p)_bd$(bd)_cutoff$(acc).txt"
    else
      fname = "XZ_dy$(dz)_dx$(dx)_p$(p)_bd$(bd)_cutoff$(acc)_id$(sim_id).txt"
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
    #pmx=0
    #ystab=YStabs(dz,dx)
    zsch,xsch=schedulemaker(dz,dx)
    #println(zsch)
    #println(xsch)
    xsch2=xreorder(xsch,dz,dx)
    println(xsch2)
    bsch=bdrysched(dz,dx)
    #println(bsch[1,:,:])
    #println(bsch[2,:,:])
    #println(bsch[3,:,:])
    #println(bsch[4,:,:])
    nzc=size(zsch)[1]
    nxc=size(xsch)[1]
    layout=SurfLayout(dz,dx,zsch,xsch2)
    #println(layout)
    layout=optlayout(layout,dz,dx)
    #println(layout)
    ql,zl,xl=linemaps(layout,dz,dx,nzc,nxc)
    #println(ql)
    #println(zl)
    #println(xl)
    #println(zsch)
    #println(xsch)
    Lz,Lx=SurfLogs(dz,dx)
    #display(Lz)
    #display(Lx)
    n=0
    f=0
    fx=0
    fy=0
    fz=0
    fn=0
    fxn=0
    fyn=0
    fzn=0
    pct=0
    totime=0
    difn=0
    diff=0
    difdm=0
    distot=0
    B1,B2=simpleblocks(dz,dx,nr,zsch,xsch2,bsch,layout,ql,zl,xl,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)
    #C=traceblock(B,layout)
    breakflag=0
    while n<nt
        #breakflag=1
        #println("pspsps")
        EZzx,EXzx,Synz,Synx=SurfMCError(dz,dx,nr,zsch,xsch2,bsch,p,al2,tmeas,k2,nth,pmz,pmx)
        #EYyx,EXyx=ZXtoYX(dz,dx,EZzx,EXzx)
        #display(EZzx)
        #display(Syny)
        #display(EXzx)
        #display(Synx)
        #display(EYyx)
        #display(EZzx)
        #display(Ey)
        PEZzx,PEXzx=getSurfPE(EZzx,EXzx,dz,dx)
        #display(PEYyx)
        #display(PEXyx)
        #display(PEY)
        L=LogComp(EZzx,EXzx,dz,dx)
        #display(L)
        #println(PEZ)

        #PEZzx,PEXzx=ZXtoYX(dz,dx,PEYyx,PEXyx)
        #display(PEZzx)
        #display(PEXzx)
        pli,plz,ply,plx=MLError(B1,B2,dz,dx,nr,PEZzx,PEXzx,Synz,Synx,zsch,xsch2,bsch,layout,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)
        ML,MLN,nfl=MLEs(pli,plz,ply,plx)

        if (ML[1]==MLN[1]) && (ML[2]==MLN[2])
            dis=0
        else
            dis=1
            distot=distot+1
        end


        resz=(ML[1]+L[1])%2
        resx=(ML[2]+L[2])%2
        reszn=(MLN[1]+L[1])%2
        resxn=(MLN[2]+L[2])%2
        #display(resy)
        #display(resx)
        le=0
        if resx==1 && resz==0
            le=1
            #println("!!!!!!!!!!!!LOGICAL X FAILURE!!!!!!!!!!!!!!")
            fx=fx+1
            f=f+1
        elseif resx==0 && resz==1
            le=1
            #println("!!!!!!!!!!!!LOGICAL Z FAILURE!!!!!!!!!!!!!!")
            fz=fz+1
            f=f+1
        elseif resx==1 && resz==1
            le=1
            #println("!!!!!!!!!!!!LOGICAL Z FAILURE!!!!!!!!!!!!!!")
            fy=fy+1
            f=f+1
        end
        len=0
        if resxn==1 && reszn==0
            len=1
            #println("!!!!!!!!!!!!LOGICAL X FAILURE!!!!!!!!!!!!!!")
            fxn=fxn+1
            fn=fn+1
        elseif resxn==0 && reszn==1
            len=1
            #println("!!!!!!!!!!!!LOGICAL Y FAILURE!!!!!!!!!!!!!!")
            fzn=fzn+1
            fn=fn+1
        elseif resxn==1 && reszn==1
            len=1
            #println("!!!!!!!!!!!!LOGICAL Z FAILURE!!!!!!!!!!!!!!")
            fyn=fyn+1
            fn=fn+1
        end

        if dis==1
            if le==1 && len==0
                difn=difn+1
            elseif le==0 && len==1
                diff=diff+1
            else
                difdm=difdm+1
            end
        end


        if len==1
            println("Here's logical Y and X error components")
            println(L)
            println("Here's Ez")
            println(EZzx)
            println("Here's Ex")
            println(EXzx)
            println("Here's PEz")
            println(PEZzx)
            println("Here's PEx")
            println(PEXzx)
            println("Here's Synz")
            println(Synz)
            println("Here's Synx")
            println(Synx)
            println("Here's Lz")
            println(Lz)
            println("Here's Lx")
            println(Lx)
            println("probability of no error")
            println(pli)
            println("probability of Z error")
            println(plz)
            println("probability of Y error")
            println(ply)
            println("probability of X error")
            println(plx)
        end

        n=n+1
        pct=pct+1
        mu=f/n
        mun=fn/n
        #println(mu)
        if (n>1)&&(f>4)
            stdev=sqrt((f*(1-mu)*(1-mu) + (n-f)*(mu)*(mu))/(n-1))
            stderr=stdev/sqrt(n)
            if stderr<(err*mu)
                breakflag=1
            end
        end
        if (pct>99)
	        atime=totime/pct
            totime=0
	        pct=0
            println("Total logical failure rate is")
            println(mu)
            println(mun)
            #println(nfail/ntr)
            println("number of trials is")
            println(n)
            println("number of logical Z failures is")
            println(fz)
            println(fzn)
            println("number of logical Y failures is")
            println(fy)
            println(fyn)
            println("number of logical X failures is")
            println(fx)
            println(fxn)
            println("total number of logical failures is")
            println(f)
            println(fn)
            println("number of disagreements")
            println(distot)
            println("number of disagreements where negating gave right answer")
            println(difn)
            println("number of disagreements where negating gave wrong answer")
            println(diff)
            println("number of disagreements where negating didn't change answer")
            println(difdm)
	        #println("average trial time")
	        #println(atime)
            if f>4
                println("stderr is")
                println(stderr)
            end
        end
        fout = open(fname,"w")
        writedlm(fout,[n,f])
        close(fout)
    end
    println(f)
    println(n)
    println(f/n)

    return
end


function SurfMCStabInv(dz,dx,nr,p,al2,tmeas,k2,nth,acc,bd,err,nt; sim_id::Int=-1)

    if sim_id < 0
      fname = "XY_dy$(dz)_dx$(dx)_p$(p)_bd$(bd)_cutoff$(acc).txt"
    else
      fname = "XY_dy$(dz)_dx$(dx)_p$(p)_bd$(bd)_cutoff$(acc)_id$(sim_id).txt"
    end
    #Random.seed!(1234 * (sim_id+2))

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
    #pmx=0
    #ystab=YStabs(dz,dx)
    zsch,xsch=schedulemaker(dz,dx)
    #println(zsch)
    #println(xsch)
    xsch2=xreorder(xsch,dz,dx)
    #println(xsch2)
    bsch=bdrysched(dz,dx)
    #println(bsch[1,:,:])
    #println(bsch[2,:,:])
    #println(bsch[3,:,:])
    #println(bsch[4,:,:])
    nzc=size(zsch)[1]
    nxc=size(xsch)[1]
    layout=SurfLayout(dz,dx,zsch,xsch2)
    #println(layout)
    layout=optlayout(layout,dz,dx)
    #println(layout)
    ql,zl,xl=linemaps(layout,dz,dx,nzc,nxc)
    #println(ql)
    #println(zl)
    #println(xl)
    #println(zsch)
    #println(xsch)
    Lz,Lx=SurfLogs(dz,dx)
    #display(Lz)
    #display(Lx)
    n=0
    f=0
    fx=0
    fy=0
    fz=0
    pct=0
    totime=0
    B1,B2=simpleblocks(dz,dx,nr,zsch,xsch2,bsch,layout,ql,zl,xl,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)
    #C=traceblock(B,layout)
    breakflag=0
    plis=zeros(Float64,16,16)
    plzs=zeros(Float64,16,16)
    plxs=zeros(Float64,16,16)
    plys=zeros(Float64,16,16)
    outcomes=zeros(Int64,16,16)
    ystab=YStabs(dz,dx)
    xstab=XStabs(dz,dx)
    PEYyxstab=zeros(Int,dz,dx)
    PEXyxstab=zeros(Int,dz,dx)
    #println("WE BEGIN")
    #while n<nt
        #breakflag=1
        #println("pspsps")
        EZzx,EXzx,Syny,Synx=SurfMCError(dz,dx,nr,zsch,xsch2,bsch,p,al2,tmeas,k2,nth,pmz,pmx)
        EYyx,EXyx=ZXtoYX(dz,dx,EZzx,EXzx)
        #display(EZzx)
        #display(Syny)
        #display(EXzx)
        #display(Synx)
        #display(EYyx)
        #display(EZzx)
        #display(Ey)
        PEYyx,PEXyx=getSurfPE(EYyx,EXyx,dz,dx)
        #display(PEYyx)
        #display(PEXyx)
        #display(PEY)
        L=LogComp(EYyx,EXyx,dz,dx)
        LY=L[1]
        LX=L[2]
        #display(L)
        #println(PEZ)
        ML=[0,0]
        println("PRE STAB LOOP")
        pli=0
        plz=0
        plx=0
        ply=0
        println("here's the errors")
        println(EYyx)
        println(EXyx)
        println("here's the pure errors")
        println(PEYyx)
        println(PEXyx)
        for sy in [1:16;]
            for sx in [1:16;]

                yst=ystab[sy]
                xst=xstab[sx]

                for iy in [1:dz;]
                    for ix in [1:dx;]
                        PEYyxstab[iy,ix]=(PEYyx[iy,ix]+yst[iy,ix])%2
                        PEXyxstab[iy,ix]=(PEXyx[iy,ix]+xst[iy,ix])%2
                    end
                end

                PEZzx,PEXzx=ZXtoYX(dz,dx,PEYyxstab,PEXyxstab)
                pli,plz,ply,plx=MLError(B1,B2,dz,dx,nr,PEZzx,PEXzx,Syny,Synx,zsch,xsch2,bsch,layout,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)
                ML,MLN,nfl=MLEs(pli,plz,ply,plx)
                MLNI=MLN[1]+2*MLN[2]
                plis[sy,sx]=pli
                plzs[sy,sx]=plz
                plxs[sy,sx]=plx
                plys[sy,sx]=ply
                outcomes[sy,sx]=MLNI
            end
        end
        println("POST STAB LOOP")
        println("Here are the Identity probabilities")
        for sy in [1:16;]
            for sx in [1:16;]
                yst=ystab[sy]
                xst=xstab[sx]
                println([sy,sx])
                #println("ystab")
                #println(ystab[sy])
                #println("xstab")
                #println(xstab[sx])
                for iy in [1:dz;]
                    for ix in [1:dx;]
                        PEYyxstab[iy,ix]=(PEYyx[iy,ix]+yst[iy,ix])%2
                        PEXyxstab[iy,ix]=(PEXyx[iy,ix]+xst[iy,ix])%2
                    end
                end
                #println("Ystab plus pure error")
                #println(PEYyxstab)
                #println("Xstab plus pure error")
                #println(PEXyxstab)
                PEZzx,PEXzx=ZXtoYX(dz,dx,PEYyxstab,PEXyxstab)
                #println("coverted pure errors")
                #println(PEZzx)
                #println(PEXzx)
                println(plis[sy,sx])
            end
        end
        println("Here are the Z probabilities")
        for sy in [1:16;]
            for sx in [1:16;]
                println([sy,sx])
                println(plzs[sy,sx])

            end
        end
        println("Here are the Y probabilities")
        for sy in [1:16;]
            for sx in [1:16;]
                println([sy,sx])
                println(plys[sy,sx])

            end
        end
        println("Here are the X probabilities")
        for sy in [1:16;]
            for sx in [1:16;]
                println([sy,sx])
                println(plxs[sy,sx])

            end
        end

        println("Here are the decoder outcomes")
        for sy in [1:16;]
            for sx in [1:16;]
                println([sy,sx])
                println(outcomes[sy,sx])

            end
        end



    return
end

dzin=parse(Int64,ARGS[1])
pin=parse(Float64,ARGS[2])
#println(pin)
ntin=parse(Int64,ARGS[3])
cutin=parse(Float64,ARGS[4])
bdin=parse(Int64,ARGS[5])
#compin=parse(Int64,ARGS[6])
repin=parse(Int64,ARGS[6])
if length(ARGS)==7
    sidin=parse(Int64,ARGS[7])
    SurfMC(dzin,3,repin,pin,8,550e-9,1e7,0,cutin,bdin,0.1,ntin,sim_id = sidin)
    #SurfMCStabInv(dzin,3,repin,pin,8,550e-9,1e7,0,cutin,bdin,0.1,ntin,sim_id = sidin)
else
    #fourgatecheck(pin,8,500e-9,1e7,0,0.001,0.001,1e-15,10000)
    SurfMC(dzin,3,repin,pin,8,550e-9,1e7,0,cutin,bdin,0.1,ntin)
    #SurfMCStabInv(dzin,3,repin,pin,8,550e-9,1e7,0,cutin,bdin,0.1,ntin)
end
