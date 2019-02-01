using FileIO
using Colors
include("AbstractVarGraph.jl")
include("VarGraph.jl")

using.AbstractVarGraph
using .VarGraph



img = load("FAU_SIG_3.jpg")
img_2 = cat(red.(img), green.(img), blue.(img);dims = 3)
f() = "dummy"


ngh = Neighborhood("grid","full",1,0)
distFct = DistanceFunction(f,f,"dummy")
wFct = WeightFunction(0,f)
data = ProblemData("grid",2,3,img_2)


b = constructGraph(data,ngh,distFct,wFct)




