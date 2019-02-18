module VarGraph
using SparseArrays
using NearestNeighbors
import NearestNeighbors: KDTree, BallTree, BruteTree, knn, inrange

export constructGraph, getadjMat, undirect


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
    VariationalGraph{T, U}
A type representing an directed weighted Graph graph.
"""
struct VariationalGraph{T<:Integer, U<:Real}
    num_verts::Int
    num_edges::Int
    edges::Array{Array{T, 1}, 1}
    weights::Array{Array{U, 1}, 1}
#     function VariationalGraph{T, U}(num_verts, edges, weights) where T<:Integer where U<:Real
#         num_edges = size(edges, 1)
#         new{T, U}(num_verts, num_edges, edges, weights)
#     end
end
#------------------------------------------------
VariationalGraph(nv::Int, edges::Array{T,2}, weights::Array{U,1}) where T<:Integer where U<:Real = 
VariationalGraph{T, U}(nv, edges, weights)
#------------------------------------------------
"""
    get_adjMat
Get adjacency Matrix.
"""
function getadjMat(g::VariationalGraph) 
    #A = Array{typeof(g.weights[1])}(undef,g.num_verts, g.num_verts)
    A = spzeros(g.num_verts, g.num_verts)
    e = g.edges
    for u = 1:g.num_verts
        for v = 1:length(e[u])
            A[u,e[u][v]] = g.weights[u][v]
        end
    end
    return A
end
#------------------------------------------------
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####
#### constructGraph
####
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
        # Extract given input
        # num_verts = div(length(data["f"]),data["dimRange"])
        num_verts = size(data["points"], 2)
        num_nghs = ngh["num_nghs"]
        inner_norm = dist_fcts["inner_norm"]
        outer_norm = dist_fcts["outer_norm"]
        weight_fct = weight_fcts["fct"]
        # Set number of edges
        num_edges = num_verts * num_nghs;
        # Initialize Array for Edgeweights
        weights = [Vector{Float64}(undef, num_nghs) for _ in 1:num_verts]
        # Check the given DataType
        if lowercase(data["type"]) == "grid"
            error("Currently unsupported!")
        elseif lowercase(data["type"]) == "point_cloud"
            # Generate kd-Tree from given Data
            points = data["points"]
            kdtree = KDTree(points);
            # Extraxt neighbors for each vertex via knn algorithm
            # One can obtain a second output "dists", which might be useful
            edges, = NearestNeighbors.knn(kdtree, points, num_nghs + 1, false)
            # Modify the obtained inds, store them into the edge array and 
            # calculate the dists
            for u = 1:num_verts
                # delete the index of u itself in the inds array
                setdiff!(edges[u],u)
                pt = points[:,u]
                # set edges and weights accordingly
                for v = 1:num_nghs
                    # Set weights via outer and inner Norm and weight function
                    inner_dist = inner_norm(pt - points[:,edges[u][v]])
                    weights[u][v] = weight_fct(outer_norm(inner_dist[1,1] + inner_dist[2,1]))
                end
            end
        end
        g = VariationalGraph(num_verts, num_edges, edges, weights)
    elseif lowercase(ngh["nType"]) == "epsilon"
        # Extract given input
        # num_verts = div(length(data["f"]),data["dimRange"])
        num_verts = size(data["points"], 2)
        epsilon = ngh["epsilon"]
        inner_norm = dist_fcts["inner_norm"]
        outer_norm = dist_fcts["outer_norm"]
        weight_fct = weight_fcts["fct"]
        # Initialize Array for Edgeweights
        weights = [Vector{Float64}() for _ in 1:num_verts]
        # Check the given DataType
        if lowercase(data["type"]) == "grid"
            error("Currently unsupported!")
        elseif lowercase(data["type"]) == "point_cloud"
            # Generate ball-Tree from given Data
            points = data["points"]
            kdtree = BallTree(points);
            # Extraxt neighbors for each vertex via a range algorithm
            edges = NearestNeighbors.inrange(kdtree, points, epsilon, false)
            # Modify the obtained inds, store them into the edge array and 
            # calculate the dists
            num_edges = 0
            for u = 1:num_verts
                # delete the index of u itself in the inds array
                setdiff!(edges[u],u)
                pt = points[:,u]
                # set edges and weights accordingly
                for v = 1:length(edges[u])
                    # Set weights via outer and inner Norm and weight function
                    inner_dist = inner_norm(pt - points[:,edges[u][v]])
                    push!(weights[u], weight_fct(outer_norm(inner_dist[1,1] + inner_dist[2,1])))
                    num_edges += 1
                end
            end
        end
        g = VariationalGraph(num_verts, num_edges, edges, weights)
    else
        error("Unsupported neighborhood type!")
    end
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
