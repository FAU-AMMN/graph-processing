module VariationalGraphs

using Distances
using LightGraphsFlows
using NearestNeighbors
using Printf
using SparseArrays

import LightGraphsFlows: 
mincut
import LightGraphs:
DiGraph, AbstractGraph, add_edge!, rem_edge!, edges,
strongly_connected_components, connected_components

import NearestNeighbors: 
KDTree, BallTree, BruteTree, knn, inrange

import Distances:
euclidean

export AbstractVariationalGraph, constructGraph, VariationalGraph, getadjMat, undirect,
generate_flowgraph, cutpursuit, primaldual,
cut_aniso, reg_aniso, 
pd_params, img_params,
fivepoint, forward

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
