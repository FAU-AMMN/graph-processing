using FileIO
using Colors
include("GridGraph.jl")

using.GridGraph: computePatchDistance2D



img = load("FAU_SIG_3.jpg")
img_2 = cat(red.(img), green.(img), blue.(img);dims = 3)
const img_3 = img_2[1:4,1:4,3]

function innerNorm(A::AbstractArray)
  A.^2
end

function outerNorm(A::AbstractArray)
  A.^(0.5)
end

display("Hello there!")

ngh = Dict([("nType","grid"),("direction","full"),("searchRadius",1),("patchRadius",2),("padMethod","replicate")])
dFct = Dict([("innerNorm",innerNorm),("outerNorm",outerNorm)])
wFct = Dict([("sigma",0),("fct",innerNorm)])
data = Dict([("nType","grid"),("dimDomain",2),("dimRange",3),("f",img_3), 
             ("nx", size(img_3,1)), ("ny", size(img_3,2))])

#@time computePatchDistance2D(data, data, ngh, wFct, dFct)