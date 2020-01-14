using VariationalGraphs
using PyPlot
using Printf
using DelimitedFiles

#A=[1 2; 3 4]

#function mat(A::Array{T, 2},i::Int64)where T <: Real
#    n,m=size(A)
#    ret=zeros(n+2*i,m+2*i)
    #Mitte
#    ret[(i+1):(i+n),(i+1):(i+m)]=A
    #Ecken
#    ret[1:i,1:i]=A[1,1]*ones(i,i)
#    ret[(n+i+1):(n+2*i),(m+i+1):(m+2*i)]=A[n,m]*ones(i,i)
#    ret[(n+i+1):(n+2*i),1:i]=A[n,1]*ones(i,i)
#    ret[1:i,(m+i+1):(m+2*i)]=A[1,m]*ones(i,i)
    #Kanten
#    ret[1:i,(i+1):(m+i)]=transpose(repeat(A[1,:],1,i))
#    ret[(n+i+1):(n+2*i),(i+1):(m+i)]=transpose(repeat(A[n,:],1,i))
#    ret[(i+1):(n+i),1:m]=repeat(A[:,1],1,i)
#    ret[(i+1):(n+i),(m+i+1):(m+2*i)]=repeat(A[:,m],1,i)

#    return ret
#end
#A=mat(A,2)
#show(A)

function rand_normal(mean, stdev)
    if stdev <= 0.0
        error("standard deviation must be positive")
    end
    u1 = rand()
    u2 = rand()
    r = sqrt( -2.0*log(u1) )
    theta = 2.0*pi*u2
    mean + stdev*r*sin(theta)
end


#printstyled(@sprintf("Hi\n"); color =:reverse)
#PyPlot.close("all")
#PyPlot.figure()
#PyPlot.imshow(I, cmap="gray")

#Ineu=VariationalGraphs.diffusion(I)
#P=zeros(2,1*)
#P[:,1:2]=I
#P[:,3:4]=Ineu
#PyPlot.figure()
#PyPlot.imshow(P, cmap="gray")

# Set data
I = readdlm("../data/cameraman.txt", Float64)
#### adds nois
n,m=size(I)
#I[1,1]=-50
nois=10*randn(n,m)
#for j=1:n
#    for k=1:m
#        nois[n,m]=rand_normal(0,300)
#    end
#end
#nois=40*rand(n,m)-20*ones(n,m)
I=I+nois
### shows pictur
PyPlot.close("all")
PyPlot.figure()
PyPlot.imshow(I, cmap="gray")

Ineu=VariationalGraphs.diffusion(1.0,1.0,I,kn())#fivepointweighted(), nonlocal()
#Ineu[1,1]=50
#show(Ineu[1:20,1:20])
PyPlot.figure()
PyPlot.imshow(Ineu, cmap="gray")
P=zeros(n,2*m)
P[:,1:m]=I
P[:,(m+1):(2*m)]=Ineu
PyPlot.figure()
PyPlot.imshow(P, cmap="gray")
#VariationalGraph([1 2;3 4], fivepoint())
