using DelimitedFiles
using VariationalGraphs
using LightGraphs
using PyPlot
using Printf
printstyled(@sprintf("Finished loading!\n"); color =:reverse)
printstyled(@sprintf("Assembling matrices and building the graph!\n"); color =:reverse)
# Loading +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set data
I = readdlm("../data/cameraman.txt", Float64)
#I = I[40:70, 60:100]
n, m = size(I)
f = reshape(I, n * m)
g = VariationalGraph(n, m, f, forward())
# Set PD Configuration ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Store parameters in img_params struct

# Set parameters for pd algorithm

setSize=ones(length(f))
par = red_params(0.35, 0.35, 1.0, 1e-07, 100.0, 1000, 500, f,setSize, g)
# Apply PD algorithm ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
x = zeros(size(f))
y = VariationalGraphs.gradient(x, g, g.config)
printstyled(@sprintf("Solving the Problem via PrimalDual!\n"); color =:reverse)

u, hist = red_primaldual( par); #@time primaldual(x, y, par);
w = reshape(u, n, m)

 printstyled(@sprintf("Finished calculating!\n"); color =:reverse)
 printstyled(@sprintf("Started visualization!\n"); color =:reverse)
 # Visualization +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 PyPlot.close("all")
 PyPlot.figure()
 PyPlot.imshow(I, cmap="gray")
 PyPlot.figure()
 PyPlot.imshow(w, cmap="gray")
