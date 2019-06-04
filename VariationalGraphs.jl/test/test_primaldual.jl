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
par = img_params(0.35, 0.35, 1.0, 1e-07, 100.0, 1000, 500, f, g)
# Apply PD algorithm ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
x = zeros(size(f))
y = VariationalGraphs.gradient(x, g, g.config)
printstyled(@sprintf("Solving the Problem via PrimalDual!\n"); color =:reverse)
#u, hist = 
@time primaldual(x, y, par);
# w = reshape(u, n, m)

# printstyled(@sprintf("Finished calculating!\n"); color =:reverse)
# printstyled(@sprintf("Started visualization!\n"); color =:reverse)
# # Visualization +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PyPlot.close("all")
# PyPlot.figure()
# PyPlot.imshow(I, cmap="gray")
# PyPlot.figure()
# PyPlot.imshow(w, cmap="gray")
# PyPlot.figure()
# PyPlot.plot(hist.energy)



# Visualization +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# function pos(n::Int64, m::Int64, u::Int64)
#     y = n - div(u - 1, m)
#     x = 1 + mod(u - 1, m) 
#     return (x, y)
# end
# n, m = size(I)
# maxi = maximum(I)
# for u = 1:g.num_verts
#     i = 1 + div(u - 1, m)
#     j = 1 + mod(u - 1, m)
#     c = string(I[i,j]/maxi)
#     (x,y) = pos(n, m, u)
#     plot(x, y, color = c, marker =:s)
#     #annotate(string(u), (points[1,u], points[2,u]))
#     for i = 1:length(g.edges_list[u])
#         v = g.edges_list[u][i]
#         (x2, y2) = pos(n, m, v)
#         plot((x, x2), (y, y2), color = "b")
#     end
# end
# Visualization +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

