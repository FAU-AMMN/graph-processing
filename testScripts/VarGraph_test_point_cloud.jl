include("VarGraph.jl")
using.VarGraph

# using Plots
# using GraphRecipes

points = rand(2,30)
F = rand(30)


inner_norm(x) = x.^2
outer_norm(x) = sqrt(x)
w_fct(x) = 1 ./ (x.^2)
f() = "dummy"

num_nghs = 4

ngh = Dict([("nType", "knn"), ("direction", "full"), ("searchRadius", 1), ("num_nghs", num_nghs)])
distFct = Dict([("inner_norm", inner_norm), ("outer_norm", outer_norm)])
wFct = Dict([("sigma", 0), ("fct", w_fct)])
data = Dict([("type", "point_cloud"), ("dimDomain", 2), ("dimRange", 1), ("f", F),("points", points)])

display("Starting to calculate")
G = constructGraph(data,ngh,distFct,wFct)

H = undirect(G)
J = getadjMat(H)


# A = [3 0 0 0; 0 5 0 0; 2 2 0 0; 0 0 2 0]
# display("Starting to plot")
# graphplot(A, x=points[1,1:4], y=points[2,1:4])
# 
# # plt = plot([points[1, :]],[points[2, :]], seriestype=:scatter);
# # for u = 1:G.num_edges
# #     p1 = [points[1,G.edges[u,1]], points[1,G.edges[u,2]]]
# #     p2 = [points[2,G.edges[u,1]], points[2,G.edges[u,2]]]
# #     plot!(plt, p1, p2)
# # end
# # display("Starting GUI")
# # gui(plt)
