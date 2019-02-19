module VariationalGraphs

using SparseArrays
using NearestNeighbors


import NearestNeighbors: 
KDTree, BallTree, BruteTree, knn, inrange

export AbstractVariationalGraph, constructGraph, getadjMat, undirect

abstract type AbstractVariationalGraph{T<:Integer, U<:Real} end

include("./knn_graph.jl")
include("./epsilon_graph.jl")
include("./variationalgraph.jl")

end # module
