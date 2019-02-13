module VarGraph
using SimpleWeightedGraphs
using NearestNeighbors
import NearestNeighbors: KDTree, BallTree, BruteTree, knn

export constructGraph

# Currently unsupported
# include("GridGraph.jl")
# using .GridGraph

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####
#### constructGraph
####
"""
    constructGraph(data::Dict, ngh::Dict, distFct::Dict, wFct::Dict)
    
Construct a graph from given data and specification. 
...
# Arguments
- `data::Dict`: Dictionary specifing the data, the graph should be constructed from. 
"type" -> String
"dimDomain" -> Integer
"dimRange" -> Integer
"f" -> AbstractArray
...
- `ngh::Dict`: the dimensions along which to perform the computation.
- "direction" -> String
- "searchRadius" -> Integer
-
...
"""
function constructGraph(data::Dict, ngh::Dict, distFct::Dict, wFct::Dict)
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
        num_verts = div(length(data["f"]),data["dimRange"])
        num_nghs = ngh["num_nghs"]
        # Array of verts
        verts = collect(1:num_verts)
        # Initialize the neighbors of nodes
        node_nghs = Array{Integer}(undef, num_verts, 3) 
        for i = 1:num_verts
            node_nghs[i, 1] = i
            node_nghs[i, 2] = (i - 1) * num_nghs + 1
            node_nghs[i, 3] = num_nghs
        end
        # Set number of edges
        num_edges = num_verts * num_nghs;
        # Initialize Array for Edges and Edgeweights
        edges = Array{Integer}(undef,num_edges,2)
        weights = Array{Float64}(undef, num_edges,1)
        # Set first column of edges to node indices
        for i = 1:num_edges
            edges[i, 1] = i%num_nghs
        end
        # Check the given DataType
        if lowercase(data["type"]) == "grid"
            error("Currently unsupported!")
        elseif lowercase(data["type"]) == "point_cloud"
            # Generate kd-Tree from given Data
            points = data["points"]
            kdtree = KDTree(points);
            # Extraxt neighbors for each vertex via knn algorithm
            for u = 1:num_verts
                indices = NearestNeighbors.knn(kdtree, points[:, u], num_nghs + 1);
            end
        end
    else
        error("Unsupported neighborhood type!")
    end
end
end
