function knn_graph(data::Dict, ngh::Dict, dist_fcts::Dict, weight_fcts::Dict)
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
end