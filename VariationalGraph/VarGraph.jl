module VarGraph
using SparseArrays
using NearestNeighbors
import NearestNeighbors: KDTree, BallTree, BruteTree, knn

export constructGraph, getadjMat, undirect


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
    VariationalGraph{T, U}
A type representing an directed weighted Graph graph.
"""
struct VariationalGraph{T<:Integer, U<:Real}
    num_verts::Int
    num_edges::Int
    edges::Array{T,2}
    weights::Array{U,1}
    function VariationalGraph{T, U}(num_verts, edges, weights) where T<:Integer where U<:Real
        num_edges = size(edges, 1)
        new{T, U}(num_verts, num_edges, edges, weights)
    end
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
    for u = 1:g.num_edges
        A[e[u,1],e[u,2]] = g.weights[u]
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
    edges = Array{Int64}(undef, num_edges, 2)
    weights = Array{Float64}(undef, num_edges)
    cols = A.colptr
    u = 1
    for i = 1:(length(cols) - 1)
        for j = cols[i]:(cols[i + 1] - 1)
            edges[j, 1] = u
            edges[j, 2] = A.rowval[j]
            weights[j] = A.nzval[j]
        end
        u += 1
    end
    VariationalGraph(g.num_verts, edges, weights)
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
        num_verts = div(length(data["f"]),data["dimRange"])
        num_nghs = ngh["num_nghs"]
        inner_norm = dist_fcts["inner_norm"]
        outer_norm = dist_fcts["outer_norm"]
        weight_fct = weight_fcts["fct"]
        # Initialize the neighbors of nodes
        node_nghs = Array{Int64}(undef, num_verts, 3) 
        for i = 1:num_verts
            node_nghs[i, 1] = i
            node_nghs[i, 2] = (i - 1) * num_nghs + 1
            node_nghs[i, 3] = num_nghs
        end
        # Set number of edges
        num_edges = num_verts * num_nghs;
        # Initialize Array for Edges and Edgeweights
        edges = Array{Int64}(undef,num_edges,2)
        weights = Array{Float64}(undef, num_edges)
        # Set first column of edges to node indices
        for i = 0:(num_edges-1)
            edges[i + 1, 1] = div(i + num_nghs,num_nghs)
        end
        # Check the given DataType
        if lowercase(data["type"]) == "grid"
            error("Currently unsupported!")
        elseif lowercase(data["type"]) == "point_cloud"
            # Generate kd-Tree from given Data
            points = data["points"]
            kdtree = KDTree(points);
            # Extraxt neighbors for each vertex via knn algorithm
            # One can obtain a second output "dists", which might be useful
            inds, = NearestNeighbors.knn(kdtree, points, num_nghs + 1, false)
            # Modify the obtained inds, store them into the edge array and 
            # calculate the dists
            for u = 1:num_verts
                # delete the index of u itself in the inds array
                setdiff!(inds[u],u)
                pt = points[:,u]
                # set edges and weights accordingly
                for v = 1:num_nghs
                    edges[(u - 1) * num_nghs + v, 2] = inds[u][v]
                    # Set weights via outer and inner Norm and weight function
                    inner_dist = inner_norm(pt - points[:,inds[u][v]])
                    weights[(u - 1) * num_nghs + v] = weight_fct(outer_norm(inner_dist[1,1] + inner_dist[2,1]))
                end

            end
        end
        # Construct the corresponding Graph
        g = VariationalGraph(num_verts, edges, weights)
    else
        error("Unsupported neighborhood type!")
    end
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
