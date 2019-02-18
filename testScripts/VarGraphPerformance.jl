using SparseArrays
using NearestNeighbors
using BenchmarkTools
import NearestNeighbors: KDTree, BallTree, BruteTree, knn, inrange

function constructknnGraph_1(data::Dict, ngh::Dict, dist_fcts::Dict, weight_fcts::Dict)
    # Extract given input
    # num_verts = div(length(data["f"]),data["dimRange"])
    num_verts = size(data["points"], 2)
    num_nghs = ngh["num_nghs"]
    inner_norm = dist_fcts["inner_norm"]
    outer_norm = dist_fcts["outer_norm"]
    weight_fct = weight_fcts["fct"]
    # Initialize the neighbors of nodes
    node_nghs = Array{Int64}(undef, 3, num_verts) 
    for i = 1:num_verts
        node_nghs[1, i] = i
        node_nghs[2, i] = (i - 1) * num_nghs + 1
        node_nghs[3, i] = num_nghs
    end
    # Set number of edges
    num_edges = num_verts * num_nghs;
    # Initialize Array for Edges and Edgeweights
    edges = Array{Int64}(undef,num_edges,2)
    weights = Array{Float64}(undef, num_edges)
    # Set first column of edges to node indices
    for i = 0:(num_edges - 1)
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
    return edges, weights, node_nghs
end

function constructknnGraph_2(data::Dict, ngh::Dict, dist_fcts::Dict, weight_fcts::Dict)
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
    return edges, weights
end



dim_domain = 3
num_pts = 10000
const points = rand(dim_domain, num_pts)
const F = rand(num_pts)
num_nghs = 4
#------------------------------------------------
inner_norm(x) = x.^2
outer_norm(x) = sqrt(x)
w_fct(x) = 1 ./ (x.^2)
f() = "dummy"
#------------------------------------------------
ngh = Dict([("nType", "knn"), ("direction", "full"), ("searchRadius", 1), ("num_nghs", num_nghs)])
dist_fcts = Dict([("inner_norm", inner_norm), ("outer_norm", outer_norm)])
weight_fcts = Dict([("sigma", 0), ("fct", w_fct)])
data = Dict([("type", "point_cloud"), ("dimDomain", dim_domain), ("dimRange", 1), ("f", F),("points", points)])



display("doing sum")

@benchmark constructknnGraph_1(data, ngh, dist_fcts, weight_fcts)

