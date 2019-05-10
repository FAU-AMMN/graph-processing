using DelimitedFiles
using VariationalGraphs
using LightGraphs
using PyPlot
using Printf
printstyled(@sprintf("Finished Loading\n"); color =:reverse)
# Loading +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
I = readdlm("../data/cameraman.txt", Float64)
#I = I[1:5,1:6]
g = VariationalGraph(I)
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

