using LightGraphs
using SimpleWeightedGraphs
using SparseArrays

function grad(data::Array, g::SimpleWeightedGraph)
    if (nv(g) != length(data))
        throw(ArgumentError("Data doesn't fit!"))
    end
    rows,cols = size(g.weights)
    gradient = spzeros(rows,cols)
    for col = 1:cols
        for ind = g.weights.colptr[col] : (g.weights.colptr[col + 1] - 1)
            gradient[g.weights.rowval[ind],col] =
                sqrt(g.weights[g.weights.rowval[ind],col]) * (data[col] - data[g.weights.rowval[ind]]);
        end
    end
    return gradient
end



function divergence(data::SparseMatrixCSC, g::SimpleWeightedGraph)
    rows,cols = size(g.weights)
    rows_2,cols_2 = size(data)
    if ((rows != rows_2) || (cols != cols_2))
        throw(ArgumentError("Data doesn't fit!"))
    end
    diver = Array{Float64}(undef, nv(g),1)
    for col = 1:cols
        for row = 1:rows
            diver[col,1] +=sqrt(g.weights[col,row]) *
                (data[row,col] - data[col,row])
        end
    end
    return diver
end
g = SimpleWeightedGraph(4)
add_edge!(g, 2, 3,0.5);
add_edge!(g, 1, 3,2.4);
DF = grad([4,5,6,6],g)
display(Matrix(DF))

EF = [[4 0 3 0];[2 0 4 2];[3 0 6 0];[2 0 0 0]]
display(EF)
DEF = divergence(sparse(EF),g)
display(DEF)
