include("GraphVis.jl")
using.GraphVis
display("finished loading stuff")
#------------------------------------------------
num_points = 30
#------------------------------------------------
const points = rand(2,num_points)
const edges = Int64[rand(1:num_points) for i in 1:(2 * num_points), j in 1:2]
vis(edges,points)