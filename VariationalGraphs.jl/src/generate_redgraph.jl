function generate_redgraph(num_vert_red::Int64,bins::Array{Int64,1},cut_edges::Array{Int64,2},
    cut_weights::Array{T,1}) where {T <: Real}

    num_edges_red=0
    weight_table=zeros(num_vert_red,num_vert_red) #eventuell als sparse

    # computes the edges and weights of the reduced graph
    #iterates over all cut edges
    for k=1:length(cut_weights)
        n=bins[cut_edges[1,k]]
        m=bins[cut_edges[2,k]]
        #if a new edge is found increase num_edges_red
        if weight_table[n,m]==0
            num_edges_red=num_edges_red+1
        end
        # adds the weight of the single edge to the reduced edge
        weight_table[n,m]=weight_table[n,m]+cut_weights[k]
    end

    #Assembles the edges_list, weights_list for reduced graph
    edges_list = [Vector{Int64}() for _ in 1:num_vert_red]
    weights_list = [Vector{Float64}() for _ in 1:num_vert_red]

    for i=1:num_vert_red            #eventuell i und j tauschen was ist mit ungerichtetem graph
        for j=1:num_vert_red
            if weight_table[i,j]!=0
                push!(edges_list[i],j )
                push!(weights_list[i],weight_table[i,j])#/edge_table[i,j]
            end
        end
    end

    conf=forward()
    return VariationalGraph(num_vert_red, num_edges_red, edges_list, weights_list, conf)
end
