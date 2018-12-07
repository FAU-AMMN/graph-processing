include("./VariationalGraph.jl")
using Main.VariationalGraphs
using SparseArrays
edges = sprand(9,9,0.5)
data = rand(nnz(edges))

ES = EdgeSet(edges)
DF = DataFunction(data)
EF = EdgeFunction{datatype(DF)}(DF,ES)

# Test that the Data stored in Edgefunction is really only a pointer

if data == EF.datafunction.data 
    println("The pointer of the original data and the data in the data function are equal.")
else
    println("The pointers are not equal. There might be a problem.")
end

