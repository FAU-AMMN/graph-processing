include("VarGraph.jl")
using.VarGraph

points = rand(2,30)
F = rand(30)


f() = "dummy"


ngh = Dict([("nType", "knn"), ("direction", "full"), ("searchRadius", 1), ("num_nghs", 5)])
distFct = Dict([("innerNorm", f), ("outerNorm", f)])
wFct = Dict([("sigma", 0), ("fct", f)])
data = Dict([("type", "point_cloud"), ("dimDomain", 2), ("dimRange", 1), ("f", F),("points", points)])


b = constructGraph(data,ngh,distFct,wFct)