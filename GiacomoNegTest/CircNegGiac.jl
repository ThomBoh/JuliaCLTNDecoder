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

    #println(ysched)
    #println(xsched)
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



function buildstrings(dy::Int,dx::Int,nr,pez,pex,Syny,Synx,layout)

    N=2*dy*dx-1
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

            qy=layout[i][2][1]
            qx=layout[i][2][2]

            if qx==dx && qy>1
                #on Y logical boundary
                #println(pez)
                #println(qz)
                #println(qx)
                ez=pez[qy,qx]%2
                ex=pex[qy,qx]%2

                if ez==0
                    #is there a Z component presenton Y logical boundary
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
                    #is there an X component presenton Y logical boundary
                    push!(outi,"Z+")
                    push!(outz,"Z-")
                    push!(outy,"Z-")
                    push!(outx,"Z+")
                elseif ex==1
                    push!(outi,"Z-")
                    push!(outz,"Z+")
                    push!(outy,"Z+")
                    push!(outx,"Z-")
                else
                    println("bad news!")
                end
                count=count+1
            elseif qx==dx && qy==1
                #on X and Y logical boundaries (corner)
                ez=pez[qy,qx]%2
                ex=pex[qy,qx]%2
                #println(e)
                if ez==0
                    #is there a Z component present on Y+X logical boundary
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
                    #is there an X component present on Y+X logical boundary
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
            elseif qy==1 && qx<dx
                #on X logical boundary
                ez=pez[qy,qx]%2
                ex=pex[qy,qx]%2
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
                    push!(outz,"Z-")
                    push!(outy,"Z+")
                    push!(outx,"Z-")
                elseif ex==1
                    push!(outi,"Z-")
                    push!(outz,"Z+")
                    push!(outy,"Z-")
                    push!(outx,"Z+")
                else
                    println("bad news!")
                end
            else
                #println(pez)
                #println(qz)
                #println(qx)
                ez=pez[qy,qx]%2
                ex=pex[qy,qx]%2
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
            #y check ancilla
            yc=layout[i][2][1]
            if i==1
                #println("hay2")
                if Syny[nr,yc]==0
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


                push!(outi,"X+")
                push!(outz,"X+")
                push!(outy,"X+")
                push!(outx,"X+")
                count=count+1
            else
                if Syny[nr,yc]==0
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
                push!(outi,"X+")
                push!(outz,"X+")
                push!(outy,"X+")
                push!(outx,"X+")
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
            push!(outi,"X+")
            push!(outz,"X+")
            push!(outy,"X+")
            push!(outx,"X+")
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
                push!(out,"X+")
                #println("hay2")
                count=count+1
            else
                push!(out,"Z+")
                push!(out,"X+")
                count=count+1
            end

        elseif typ==2
            push!(out,"Z+")
            push!(out,"X+")
            #x check ancilla
            count=count+1
        end

    end
    #println(outi)
    #println(outz)
    return out

end

function buildCirc(b,dz,dx,nr,zsch,xsch,bsch,layout,ql,zl,xl,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)

    #First initialize all hardware noise parameters and failure rates. Note that p=k1/k2

    ppaxz=(15/2)*p
    ppaxy=0#(15/2)*p

    ppazx=0.39*exp(-4*al2)
    ppazy=0.39*exp(-4*al2)

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
                #y check ancilla qubit
                push!(gates,("idle",(2*q-1,2*q),(pz=ppaxz,px=0,py=ppaxy)))

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

                    #Y check ancilla qubit, identifies Z errors, data is target, ancilla is control

                    zc=layout[q][2][1] #which z check it is in the list of z check schedules
                    partner=zsch[zc][cnt]

                    if partner==[-1,-1]
                        #wait location for this ancilla
                        push!(gates,("idle",(2*q-1,2*q),(pz=pzic,px=pxic,py=pyic)))
                    else
                        pzad=partner[1]
                        pxad=partner[2]
                        q2=ql[pzad,pxad]
                        push!(gates,("CYT",(2*q-1,2*q,2*q2-1,2*q2),(px2=CYpx2,pz2=CYpz2,py2=CYpy2,px1=CYpx1,px1x2=CYpx1x2,px1z2=CYpx1z2,px1y2=CYpx1y2,pz1=CYpz1,pz1x2=CYpz1x2,pz1z2=CYpz1z2,pz1y2=CYpz1y2,py1=CYpy1,py1x2=CYpy1x2,py1z2=CYpy1z2,py1y2=CYpy1y2)))
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
                #y check ancilla qubit
                push!(gates,("idle",(2*q-1,2*q),(pz=pmx,px=0,py=0)))

            elseif typ==2
                #x check ancilla qubit
                push!(gates,("idle",(2*q-1,2*q),(pz=pmx,px=0,py=0)))

            end

        end
    end
    istr=buildinitstring(dz,dx,nr,layout)
    #println(istr)
    T=productstate(2*nq) #all zeros initial state sets initial error configuration to be trivial
    TI=productstate(siteinds(T),istr)
    #display(T)
    TOut=runcircuit(TI,gates,cutoff=acc,maxdim=bd)

    return TOut,gates #return the MPS, whose physical indices are as described at top of function
end

function NegCirc(dz,dx,nr,p,al2,tmeas,k2,nth,acc,bd,err,nt; sim_id::Int=-1)


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

    xsch2=xreorder(xsch,dz,dx)

    bsch=bdrysched(dz,dx)

    nzc=size(zsch)[1]
    nxc=size(xsch)[1]
    layout=SurfLayout(dz,dx,zsch,xsch2)

    layout=optlayout(layout,dz,dx)

    ql,zl,xl=linemaps(layout,dz,dx,nzc,nxc)
    CircMPS,gates=buildCirc(1,dz,dx,nr,zsch,xsch,bsch,layout,ql,zl,xl,p,al2,tmeas,k2,nth,pmz,pmx,acc,bd)

    return CircMPS,gates
end

#Calling this function returns an MPS created by runcircuit and the circuit used by
#runcirctuit to do it. The MPS is the ML decoder
#for the XY surface code with a single measurement round. Different logical cosets can
#be computed by taking the inner product of this MPS with certain product state MPS

#This MPS has many negative entries whose positive versions will hopefully tell us
#what's going on

CircMPS,gates=NegCirc(3,3,1,1e-4,8,550e-9,1e7,0,1e-15,40,0.1,0,sim_id =1)
