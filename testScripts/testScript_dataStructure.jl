using LightGraphs

#Simple function to switch the represantion of a graph from
#(A) fadjlist, ne
#to
#(B) edgeList, colOffs.
#(A) uses Vector{Vectors{T}} as datastrucure, this might be ineffecient
#for numerical usage(check this!), however (B)
#interprets an undirected Graph as a directed Graph, hence we double
#the amount of storage(This could be omitted).
function getRepr(g::SimpleGraph)
    edgeList = Array{eltype(g)}(undef,2 * g.ne)
    colOffs = Array{eltype(g)}(undef,nv(g))
    tmp = 0
    for vert = 1 : nv(g)
        for nvert = 1 : length(g.fadjlist[vert])
            edgeList[vert + nvert - 1] = g.fadjlist[vert][nvert]
        end
        tmp += length(g.fadjlist[vert])
        colOffs[vert] = tmp;
    end
    return edgeList,colOffs
end

G = SimpleGraph(3)
add_edge!(G, 1, 2)
add_edge!(G, 2, 3)
add_edge!(G, 2, 2)
e, c = getRepr(G)
