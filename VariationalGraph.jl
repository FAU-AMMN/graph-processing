module VariationalGraphs

using SparseArrays
import Base.length
import Core.typeof
export VertexSet, EdgeSet, DataFunction, VertexFunction, EdgeFunction
export datatype, length
# -- General Data structures -- #

# DataFunction: Data function that can correspond to an EdgeSet or a VertexSet
struct DataFunction{T<:Real}
    data::Array{T}
end

#VertexSet: The structure of vertices given by an array with integers
struct VertexSet 
    vertices::Array{Int}
end

# EdgeSet: The structure of edges is given by a sparse matrix 
struct EdgeSet
    edges::SparseMatrixCSC
    numOfNBRs::Array{Int}
    numOfEdges::Int
    function EdgeSet(edges)
        numOfNBRs = computeNumOfNBRs(edges)
        numOfEdges = nnz(edges)
        new(edges, numOfNBRs, numOfEdges)
    end
end

# EdgeFunction: A structure given by an DataFunctions that corresponds to the given EdgeSet
struct EdgeFunction{T<:Real}
    datafunction::DataFunction{T}
    edgeset::EdgeSet
    function EdgeFunction{T}(datafunction::DataFunction{T}, edgeset::EdgeSet) where {T<:Real} 
        if length(datafunction) == edgeset.numOfEdges
            new(datafunction, edgeset)
        else
            error("Number of edges and elements in data are not the same!")
        end
    end
end

# VertexFunction: Stores a datafunction corresponding to a vertex set
struct VertexFunction{T<:Real}
    vertexset::VertexSet
    datafunction::DataFunction{T} 
    function VertexFunction{T}(datafunction, vertexset) where {T<:Real}
        if length(datafunction) == length(vertexset)
            new(datafunction, vertexset)
        else
            error("Number of vertices and elements of data are not the same!")
        end
    end
end
# -- Methods for structures -- #

# - VertexSet - #
function length(vertexset::VertexSet)
    return length(vertexset.vertices)
end

# - EdgeSet - # 
# Get number of edges
function numOfEdges(ES::EdgeSet)
    return nnz(ES.edges)
end
# Get number of neighbors of each vertex 
function numOfNeighbors(ES::EdgeSet)
    return ES.numOfNeighbors
end

# Compute number of neighbors 
function computeNumOfNBRs(edges::SparseMatrixCSC)
    v = zeros(Int,edges.n,1)
    for row in edges.rowval 
        v[row] += 1
    end
    return v
end

# - DataFunction - #

# Compute length of the stored data
function length(datafunction::DataFunction)
    return length(datafunction.data)
end

# Get the data type of the stored data
function datatype(datafunction::DataFunction)
    return typeof(datafunction.data[1])
end

# - EdgeFunction - #

end # End module