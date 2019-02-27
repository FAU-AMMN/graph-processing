module VariationalGraphs

using SparseArrays
using NearestNeighbors
using LightGraphsFlows

import LightGraphsFlows: 
mincut
import LightGraphs:
DiGraph, add_edge!

import NearestNeighbors: 
KDTree, BallTree, BruteTree, knn, inrange

export AbstractVariationalGraph, constructGraph, getadjMat, undirect,
generate_flowgraph, cutpursuit, get_alt_edgerep

abstract type AbstractVariationalGraph{T<:Integer, U<:Real} end

include("./knn_graph.jl")
include("./epsilon_graph.jl")
include("./variationalgraph.jl")
include("./generate_flowgraph.jl")
include("./cutpursuit.jl")

end # module
