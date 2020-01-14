using VariationalGraphs
using PyPlot
using Printf
using DelimitedFiles

# Set data
predata=readdlm("../data/cameraman.txt", Float64)
I = predata[50:200,50:200]#rand(30,30)#readdlm("../data/cameraman.txt", Float64)
#I=I./200
show(maximum(I))
show(minimum(I))
n,m=size(I)
PyPlot.close("all")
PyPlot.figure()
PyPlot.imshow(I, cmap="gray")

#maybethis
#Im=zeros(n,m)        #This is a workaround
#Im=Im+I'              #for a reshape problem
#f0=reshape(Im,n*m)

##
#n, m = size(I)
#f = reshape(I, n * m)
#g = VariationalGraph(n, m, f, forward())
##
#w=zeros(n,m)        #This is a workaround
#w=w+I'              #for a reshape problem

#f0=reshape(w, n*m)

f0=reshape(I,n*m)#this reshape isprobably wrong
#g = VariationalGraph(n, m, f0, forward())
g=VariationalGraph(I, fivepoint())
Image=cutpursuit2(1, f0, g, g.weights_mat, 100.0, 2, cut_aniso())
printstyled(@sprintf("Maximum:\n"); color =:reverse)
show(maximum(Image))
#
#Im=zeros(n,m)        #This is a workaround
#Im=Im+reshape(Image,n,m)' #for a reshape problem

Im=reshape(Image,n,m)

#PyPlot.close("all")
PyPlot.figure()
PyPlot.imshow(Im, cmap="gray")
