module VariationalGraphs

using SparseArrays
using NearestNeighbors
using LightGraphsFlows
using Printf

import LightGraphsFlows: 
mincut
import LightGraphs:
DiGraph, AbstractGraph, add_edge!, rem_edge!, edges,
strongly_connected_components, connected_components

import NearestNeighbors: 
KDTree, BallTree, BruteTree, knn, inrange

export AbstractVariationalGraph, constructGraph, getadjMat, undirect,
generate_flowgraph, cutpursuit, get_alt_edgerep,
cut_aniso, reg_aniso

###############################################################################

abstract type AbstractVariationalGraph{T<:Integer, U<:Real} end

abstract type cut_type end

abstract type reg_type end

struct cut_aniso <: cut_type end
struct cut_iso <: cut_type end

struct reg_aniso <: cut_type end

###############################################################################

include("./knn_graph.jl")
include("./epsilon_graph.jl")
include("./variationalgraph.jl")
include("./generate_flowgraph.jl")
include("./primaldual.jl")
include("./cutpursuit.jl")

end # module
