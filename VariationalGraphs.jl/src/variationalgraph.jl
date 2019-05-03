###################################################################################################
"""
```julia
    GraphConfig
```
The struct used to dispatch different graph configurations.
"""
abstract type GraphConfig end
struct vargraphdef <: GraphConfig end
###################################################################################################
# The struct representing a Graph and outer constructor functions.
"""
    VariationalGraph{T, U}
A type representing an directed weighted Graph graph.
"""
struct VariationalGraph{T <: Integer, U <: Real} <: AbstractVariationalGraph{T, U}
    num_verts::Int
    num_edges::Int
    #
    edges_list::Array{Array{T, 1}, 1}
    edges_mat::Array{T, 2}
    weights_list::Array{Array{U, 1}, 1}
    weights_mat::Array{U, 1}
    #
    config::GraphConfig
end
#------------------------------------------------
VariationalGraph(num_verts::T, num_edges::T, edge_list::Array{Array{T, 1}, 1}, edge_mat::Array{T, 2},
                          weight_list::Array{Array{U, 1}, 1}, weights_mat::Array{T, 2}) where {T <: Integer, U <: Real} =
VariationalGraph(num_verts, num_edges, edge_list, edge_mat, weight_list, weight_mat, vargraphdef())
#------------------------------------------------
function VariationalGraph(num_verts::T, num_edges::T, edge_list::Array{Array{T, 1}, 1}, 
                          weight_list::Array{Array{U, 1}, 1}, config::GraphConfig) where {T <: Integer, U <: Real}
    edge_mat, weight_mat = list2mat(num_edges, edge_list, weight_list)
    return VariationalGraph(num_verts, num_edges, edge_list, edge_mat, weight_list, weight_mat, config)
end
#------------------------------------------------
VariationalGraph(num_verts::T, num_edges::T, edge_list::Array{Array{T, 1}, 1}, 
                 weight_list::Array{Array{U, 1}, 1}) where {T <: Integer, U <: Real} =
VariationalGraph(num_verts::T, num_edges::T, edge_list::Array{Array{T, 1}, 1}, 
                 weight_list::Array{Array{U, 1}, 1}, vargraphdef())                 
#------------------------------------------------
VariationalGraph(num_edges::T, edges::Array{Array{T, 1}, 1}, weights::Array{Array{U, 1}, 1}) where {T <: Integer, U <: Real} = 
VariationalGraph(length(edges), num_edges, edges, weights)
#------------------------------------------------
###################################################################################################
"""
    rm_edge!
Remove edge from edgeset.
"""
function rm_edge!(edges::Array{Array{T, 1}, 1}, u::T, v::T) where T <: Integer
    setdiff!(edges[u], v)
end

#rm_edge!(g::VariationalGraph, u::Int, v::Int) =  
#rm_edge!(g.edges, u::Int, v::Int)

###################################################################################################

"""
    add_edge!
Add edge to an edgeset.
"""
function add_edge!(edges::Array{Array{T, 1}, 1}, u::T, v::T) where T <: Integer
    if !(v in edges[u])
        push!(edges[u], v)
    end
end

#add_edge!(g::VariationalGraph, u::Int, v::Int) =  
#add_edge!(g.edges, u::Int, v::Int)

################################################################################################### 

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

###################################################################################################

"""
```julia
    list2mat()
```
Get alternative edge representation.
"""
function list2mat(num_edges::T, e::Array{Array{T, 1}}, w::Array{Array{U, 1}, 1}) where{T<:Int64, U<:Float64} 
    edges = Array{Int64, 2}(undef, 2, num_edges)
    weights = Array{Int64, 1}(undef, 1, num_edges)
    i = 1
    for u = 1:length(e)
        for v = 1:length(e[u])
            weights[i] = w[u][v]
            e[1, i] = u
            e[2, i] = e[u][v]
            i += 1
        end
    end
    return edges, weights
end

###################################################################################################

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

###################################################################################################

grad(d::U, f::Array{T, 1}, g::VariationalGraph{T, U}) where {T<: Real, U<:Int} = 
grad(d, f, g, g.config)

###################################################################################################

function grad(d::U, f::Array{T, 1}, g::VariationalGraph{U, T}, config::vargraphdef) where {T<: Real, U<:Integer}
    m = div(size(f, 1), d)
    gradf = Array{T, 2}(undef, d, size(g.edges_mat, 2))
    for l = 1:d
        for i = 1:size(g.edges_mat, 2)
            gradf[l, i] = g.weights_mat[i] .* (f[g.edges_mat[2, i] +  (l - 1) * m] - f[g.edges_mat[1, i] + (l - 1) * m]);
        end
    end
    return gradf
end

###################################################################################################
divergence(d::U, m::U, gradf::Array{T, 2}, g::VariationalGraph{U, T}) where {T<: Real, U<:Integer} = 
diveregence(d, m, gradf, g, g.config)
###################################################################################################

function divergence(d::U, m::U, gradf::Array{T, 2}, g::VariationalGraph{U, T}, config::vargraphdef) where {T<: Real, U<:Integer}
    divy = Array{T, 1}(undef, d * m)
    for l = 1:d
        for j = 1:size(E, 2)
            divy[E[1, j] + (l - 1) * m] += 2 * w[j] * gradf[l, j]
        end
    end
    return divy
end

###############################################################################


