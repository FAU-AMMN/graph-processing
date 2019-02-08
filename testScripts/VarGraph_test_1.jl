using FileIO
using Colors
include("VarGraph.jl")

using.VarGraph



img = load("FAU_SIG_3.jpg")
img_2 = cat(red.(img), green.(img), blue.(img);dims = 3)
f() = "dummy"


ngh = Dict([("nType","grid"),("direction","full"),("searchRadius",1),("patchRadius",0)])
distFct = Dict([("innerNorm",f),("outerNorm",f),("padMethod","dummy")])
wFct = Dict([("sigma",0),("fct",f)])
data = Dict([("nType","grid"),("dimDomain",2),("dimRange",3),("f",img_2)])


b = constructGraph(data,ngh,distFct,wFct)
