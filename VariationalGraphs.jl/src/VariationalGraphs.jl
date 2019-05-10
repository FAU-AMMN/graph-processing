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

export AbstractVariationalGraph, constructGraph, VariationalGraph, getadjMat, undirect,
generate_flowgraph, cutpursuit,
cut_aniso, reg_aniso

###############################################################################

abstract type AbstractVariationalGraph{T<:Integer, U<:Real} end

# structs to dispatch the cut type
abstract type cut_type end
struct cut_aniso <: cut_type end
struct cut_iso <: cut_type end

# structs to dispatch the regularization type
abstract type reg_type end
struct reg_aniso <: cut_type end
struct reg_iso <: cut_type end
struct reg_aniso_aniso <: cut_type end
struct reg_aniso_iso <: cut_type end
struct reg_iso_aniso <: cut_type end
struct reg_iso_iso <: cut_type end

###############################################################################

include("./variationalgraph.jl")
include("./constructgraph.jl")
include("./generate_flowgraph.jl")
include("./primaldual.jl")
include("./cutpursuit.jl")

end # module
