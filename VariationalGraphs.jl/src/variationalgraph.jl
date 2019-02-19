###############################################################################
# The struct representing a Graph and outer constructor functions.
"""
    VariationalGraph{T, U}
A type representing an directed weighted Graph graph.
"""
struct VariationalGraph{T <: Integer, U <: Real} <: AbstractVariationalGraph{T, U}
    num_verts::Int
    num_edges::Int
    edges::Array{Array{T, 1}, 1}
    weights::Array{Array{U, 1}, 1}
end
#------------------------------------------------
VariationalGraph(num_edges, edges::Array{Array{T, 1}, 1}, 
weights::Array{Array{U, 1}, 1}) where T <: Integer where U <: Real = 
VariationalGraph(length(edges), num_edges, edges, weights)
#------------------------------------------------
###############################################################################
"""
    rm_edge!
Remove edge from edgeset.
"""
function rm_edge!(edges::Array{Array{T, 1}, 1}, u::T, v::T) where T <: Integer
    setdiff!(edges[u], v)
end

rm_edge!(g::VariationalGraph, u::Int, v::Int) =  
rm_edge!(g.edges, u::Int, v::Int)

###############################################################################

"""
    add_edge!
Add edge to an edgeset.
"""
function add_edge!(edges::Array{Array{T, 1}, 1}, u::T, v::T) where T <: Integer
    if !(v in edges[u])
        push!(edges[u], v)
    end
end

add_edge!(g::VariationalGraph, u::Int, v::Int) =  
add_edge!(g.edges, u::Int, v::Int)

############################################################################### 

"""
    get_adjMat
Get adjacency Matrix.
"""
function getadjMat(g::VariationalGraph) 
    # A = Array{typeof(g.weights[1])}(undef,g.num_verts, g.num_verts)
    A = spzeros(g.num_verts, g.num_verts)
    e = g.edges
    for u = 1:g.num_verts
        for v = 1:length(e[u])
            A[u,e[u][v]] = g.weights[u][v]
        end
    end
    return A
end

###############################################################################

"""
    undirect
Undirect a directed Graph in a *very* inefficent fashion.
"""
function undirect(g::VariationalGraph)
    A = getadjMat(g)
    A = 0.5 * (A + transpose(A))
    num_edges = nnz(A)
    edges = [Vector{Int64}() for _ in 1:g.num_verts]
    weights = [Vector{Float64}() for _ in 1:g.num_verts]
    cols = A.colptr
    u = 1
    for i = 1:(length(cols) - 1)
        for j = cols[i]:(cols[i + 1] - 1)
            push!(edges[u],A.rowval[j])
            push!(weights[u],A.nzval[j]) 
        end
        u += 1
    end
    VariationalGraph(g.num_verts, num_edges, edges, weights)
end


###############################################################################
# constructGraph

"""
    constructGraph(data::Dict, ngh::Dict, distFct::Dict, wFct::Dict)
    
Construct a graph from given data and specification. 
...
# Arguments
- `data::Dict`: Dictionary specifing the data, the graph should be constructed from. 
- `ngh::Dict`: Dictionary specifing the neighborhood of each vertex in the graph. 
...
"""
function constructGraph(data::Dict, ngh::Dict, dist_fcts::Dict, weight_fcts::Dict)
    # Grid-Grid
    if lowercase(ngh["nType"]) == "grid" && lowercase(data["type"]) == "grid"
        error("Currently unsupported!")
        #----------------------------------------------------------------
        #Get the indices that correspond to the neighbors specified by ngh and data in the stencil
        relNghInd = getRelNghInd(lowercase(ngh["direction"]), data["dimDomain"], ngh["searchRadius"])
        #----------------------------------------------------------------
        #Define useful values
        num_verts = div(length(data["f"]),data["dimRange"])
        num_nghs = length(relNghInd)
        num_edges = numOfVerts * numOfNgh
        g = SimpleWeightedDiGraph(numOfVerts)
    # knn - PointCloud
    elseif lowercase(ngh["nType"]) == "knn"
        knn_graph(data, ngh, dist_fcts, weight_fcts)
    elseif lowercase(ngh["nType"]) == "epsilon"
        epsilon_graph(data, ngh, dist_fcts, weight_fcts)
    else
        error("Unsupported neighborhood type!")
    end
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
