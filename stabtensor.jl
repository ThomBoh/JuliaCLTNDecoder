using Printf
using PastaQ
using ITensors
using Random
import PastaQ: gate
import PastaQ.ITensors: array
using DelimitedFiles
using LinearAlgebra

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

Id4 =1* Matrix(I, 4, 4)
