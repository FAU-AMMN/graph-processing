include("VarGraph.jl")
using.VarGraph
#------------------------------------------------
# This takes some time, include only for visualization
include("GraphVis.jl")
using.GraphVis
#------------------------------------------------
# Visualization works for d = 2,3 
# Calculation should work for "arbitrary" positive integers
dim_domain = 2
num_nghs = 4
#------------------------------------------------
const points = rand(dim_domain,30)
const F = rand(30)
#------------------------------------------------
inner_norm(x) = x.^2
outer_norm(x) = sqrt(x)
w_fct(x) = 1 ./ (x.^2)
f() = "dummy"
#------------------------------------------------
ngh = Dict([("nType", "knn"), ("direction", "full"), ("searchRadius", 1), ("num_nghs", num_nghs)])
distFct = Dict([("inner_norm", inner_norm), ("outer_norm", outer_norm)])
wFct = Dict([("sigma", 0), ("fct", w_fct)])
data = Dict([("type", "point_cloud"), ("dimDomain", dim_domain), ("dimRange", 1), ("f", F),("points", points)])
#------------------------------------------------
display("Starting to calculate!")
G = constructGraph(data,ngh,distFct,wFct)
#------------------------------------------------
H = undirect(G)
display("Starting to plot! Have a coffee in the meantime ;)")
vis(H.edges, points) 
#------------------------------------------------
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
