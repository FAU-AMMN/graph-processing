###################################################################################################
# constructGraph

"""
```julia
    constructGraph(data::Dict, ngh::Dict, distFct::Dict, wFct::Dict)
```
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

###################################################################################################
struct fivepoint<:GraphConfig end
struct forward<:GraphConfig end
struct nonlocal<:GraphConfig end            # Ich hoffe richtig exportiert
struct fivepointweighted<:GraphConfig end   #
struct kn<:GraphConfig end

function VariationalGraph(f::Array{T, 2},conf::kn)where T <: Real
    loc_n=3
    num_nghs=10
    n, m = size(f)
    num_verts = n * m
    num_edges= num_verts*num_nghs
    f=mat(f,loc_n)
    points = zeros((2*loc_n+1)^2,num_verts)

    for v = 1:num_verts
        i = 1 + div(v - 1, m)
        j = 1 + mod(v - 1, m)
        points[:,v]=reshape(f[i:(2*loc_n+i),j:(2*loc_n+j)],(2*loc_n+1)^2,1)
    end
    kdtree = KDTree(points);
    edges_list, weights_list = NearestNeighbors.knn(kdtree, points, num_nghs + 1, true)
    for v = 1:num_verts
        deleteat!(edges_list[v],1)
        deleteat!(weights_list[v],1)
        for k=1:num_nghs
            weights_list[v][k]=exp(-weights_list[v][k]/50^2)
        end
    end
    return VariationalGraph(num_verts, num_edges, edges_list, weights_list, conf)
end

function mat(A::Array{T, 2},i::Int64)where T <: Real
    n,m=size(A)
    ret=zeros(n+2*i,m+2*i)
    #Mitte
    ret[(i+1):(i+n),(i+1):(i+m)]=A
    #Ecken
    ret[1:i,1:i]=A[1,1]*ones(i,i)
    ret[(n+i+1):(n+2*i),(m+i+1):(m+2*i)]=A[n,m]*ones(i,i)
    ret[(n+i+1):(n+2*i),1:i]=A[n,1]*ones(i,i)
    ret[1:i,(m+i+1):(m+2*i)]=A[1,m]*ones(i,i)
    #Kanten
    ret[1:i,(i+1):(m+i)]=transpose(repeat(A[1,:],1,i))
    ret[(n+i+1):(n+2*i),(i+1):(m+i)]=transpose(repeat(A[n,:],1,i))
    ret[(i+1):(n+i),1:i]=repeat(A[:,1],1,i)
    ret[(i+1):(n+i),(m+i+1):(m+2*i)]=repeat(A[:,m],1,i)

    return ret
end

function VariationalGraph(f::Array{T, 2}, conf::nonlocal) where T <: Real
    loc_n=3;#size of local neighborhood
    glob_n=10;#size of global neighborhood
    n, m = size(f)
    num_verts = n * m
    num_edges = 0
    edges_list = [Vector{Int64}() for _ in 1:num_verts]
    weights_list = [Vector{Float64}() for _ in 1:num_verts]
    loc_neighborhood=[Array{Float64,2}(undef,2*loc_n+1,2*loc_n+1) for _ in 1:num_verts]

    for v = 1:num_verts
        i = 1 + div(v - 1, m)
        j = 1 + mod(v - 1, m)
        neighborhood=NaN*zeros(2*loc_n+1,2*loc_n+1)#vielleicht mit not a number
        for k = -loc_n:loc_n
            for l = -loc_n:loc_n
                if  (1 <= i + k <= n) && (1 <= j + l <= m) #node itself is not included in neigbours
                    neighborhood[1+loc_n+k,1+loc_n+l]=f[i + k,j + l]
                    #push!(edges_list[v], (i + k - 1) * m + (j + l))
                    #push!(weights_list[v], 1)
                    #num_edges += 1
                end
            end
        end
        loc_neighborhood[v]=neighborhood
        #push!(loc_neighborhood[v],neighborhood)
    end

    for v = 1:num_verts
        i = 1 + div(v - 1, m)
        j = 1 + mod(v - 1, m)
        for k = -glob_n:glob_n
            for l = -glob_n:glob_n
                if (1 <= i + k <= n) && (1 <= j + l <= m)&&(j!=0||k!=0) #node itself is not included in neigbours
                    push!(edges_list[v], (i + k - 1) * m + (j + l))
                    push!(weights_list[v],
                     exp(-helpernorm(loc_neighborhood[v],loc_neighborhood[(i + k - 1) * m + (j + l)])/50^2))
                    num_edges += 1
                end
            end
        end
    end
    #show(edges_list[1])
    return VariationalGraph(num_verts, num_edges, edges_list, weights_list, conf)
end



function helpernorm(x::Array{T, 2},y::Array{T, 2})where T <: Real
    vec=x-y
    n,m=size(vec)
    k=sum(vec.==vec)

    vec=replace!(vec,NaN=>0)
    ret=(n*m/(k))*sum((vec.^2))
    #show(ret)
    #if ret<0
    #    show(1)#prüft ret
    #end

    return ret #mit oder ohne sqrt?
    #was passiert am Rand vom Bild?
    #(n*m/(n*m-k)) Korrekturfaktor für nicht enthaltenen Elemente
end

function Laplace(f::Array{T, 1}, g::VariationalGraph, conf::nonlocal) where T <: Real
    #show(6)
    #show(g.edges_list[1])
    #show(length(g.edges_list[1]))
    lap = zeros(g.num_verts)#Array{T, 1}(undef, g.num_verts)
    #show(g.num_verts)
    #show(" ")
    #show(1:(g.num_verts))
    #show(" ")
    nummer=g.num_verts
    for v=1:nummer
        divider=1  #oder 0
        #show(length(g.edges_list[v]))
        #show(" ")
        for n=1:length(g.edges_list[v])
            #if g.weights_list[v][n]>200
            #    show(2)
            #end

            divider=divider+g.weights_list[v][n]

            lap[v]=lap[v]+g.weights_list[v][n]*(f[g.edges_list[v][n]]-f[v])
        end
        #if divider>200||divider<1
        #    show(9)
        #end
        #if v==1
        #    show(" ")
        #    show(divider)
        #end

        lap[v]=lap[v]/divider
        ### vectorizes
        #lap[v]=lap[v]+sum(g.weights_list[v].*(f[g.edges_list[v]].-f[v])
    end
    return lap
end

function VariationalGraph(f::Array{T, 2}, conf::fivepointweighted) where T <: Real
    n, m = size(f)
    num_verts = n * m
    num_edges = 0
    edges_list = [Vector{Int64}() for _ in 1:num_verts]
    weights_list = [Vector{Float64}() for _ in 1:num_verts]
    for v = 1:num_verts
        i = 1 + div(v - 1, m)
        j = 1 + mod(v - 1, m)
        for k = -1:1
            for l = -1:1
                if (abs(k + l) == 1) && (1 <= i + k <= n) && (1 <= j + l <= m) #node itself is not included in neigbours
                    push!(edges_list[v], (i + k - 1) * m + (j + l))
                    push!(weights_list[v], exp(-abs(f[i+k,j+l]-f[i,j])/20^2) )#stimmen die weights oder mit norm?
#1/exp(abs((f[i+k,j+l]-f[i,j])-40)\20)
                    num_edges += 1
                end
            end
        end
    end
    return VariationalGraph(num_verts, num_edges, edges_list, weights_list, conf)
end


function VariationalGraph(f::Array{T, 2}, conf::fivepoint) where T <: Real
    n, m = size(f)
    num_verts = n * m
    num_edges = 0
    edges_list = [Vector{Int64}() for _ in 1:num_verts]
    weights_list = [Vector{Float64}() for _ in 1:num_verts]
    for v = 1:num_verts
        i = 1 + div(v - 1, m)
        j = 1 + mod(v - 1, m)
        for k = -1:1
            for l = -1:1
                if (abs(k + l) == 1) && (1 <= i + k <= n) && (1 <= j + l <= m) #node itself is not included in neigbours
                    push!(edges_list[v], (i + k - 1) * m + (j + l))
                    push!(weights_list[v], 1)
                    num_edges += 1
                end
            end
        end
    end
    return VariationalGraph(num_verts, num_edges, edges_list, weights_list, conf)
end
#brauch ich hier überhaupt fivepoint?
function Laplace(f::Array{T, 1}, g::VariationalGraph, conf::fivepoint) where T <: Real
    lap = zeros(g.num_verts)
    for v=1:(g.num_verts)
        for n=1:length(g.edges_list[v])
            lap[v]=lap[v]+g.weights_list[v][n]*(f[g.edges_list[v][n]]-f[v])
        end
        ### vectorizes
        #lap[v]=lap[v]+sum(g.weights_list[v].*(f[g.edges_list[v]].-f[v])
    end
    return lap
end

###################################################################################################

function VariationalGraph(n::Int64, m::Int64, f::Array{T, 1}, conf::forward) where T <: Real
    num_verts = n * m
    num_edges = 0
    edges_list = [Vector{Int64}() for _ in 1:num_verts]
    weights_list = [Vector{Float64}() for _ in 1:num_verts]
    for v = 1:num_verts
        i = 1 + div(v - 1, n)
        j = 1 + mod(v - 1, n)
        for k = 0:1
            for l = 0:1
                if (abs(k + l) == 1) && (1 <= i + k <= n) && (1 <= j + l <= m) #node itself is not included in neigbours
                    push!(edges_list[v], (i + k - 1) * n + (j + l))
                    push!(weights_list[v], 1)
                    num_edges += 1
                end
            end
        end
        # Neumann boundary condition
        while length(edges_list[v]) < 2
           push!(edges_list[v], v)
           push!(weights_list[v], 1)
           num_edges += 1
        end
    end
    return VariationalGraph(num_verts, num_edges, edges_list, weights_list, conf)
end

function VariationalGraph(f::Array{T, 2}, conf::forward) where T <: Real
    n, m = size(f)
    return VariationalGraph(n, m, reshape(f, n*m), conf)
end

function gradient(x::Array{T, 1}, g::VariationalGraph, conf::forward) where T <: Real
    gradx = Array{T, 1}(undef, g.num_edges)
    for i = 1:g.num_edges
        @inbounds gradx[i] = g.weights_mat[i] * (x[g.edges_mat[2, i]] - x[g.edges_mat[1, i]])
    end
    return gradx
end

function divergence(y::Array{T, 1}, g::VariationalGraph, conf::forward) where T <: Real
    divy = zeros(g.num_verts)
    for i = 1:g.num_edges
        @inbounds divy[g.edges_mat[1, i]] += g.weights_mat[i] * y[i]
        @inbounds divy[g.edges_mat[2, i]] -= g.weights_mat[i] * y[i]
    end
    return divy
end
###################################################################################################
"""
```julia
    epsilon_graph()
```
ToDo: Write DocString!
"""
function epsilon_graph(data::Dict, ngh::Dict, dist_fcts::Dict, weight_fcts::Dict)
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
end

###################################################################################################
"""
```julia
    knn_graph()
```
ToDo: Write DocString!
"""
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

###################################################################################################
