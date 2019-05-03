###################################################################################################

function set_red_graph!(E_red::Array{U, 2}, w_red::Array{T, 1},
                                         E::Array{U, 2}, w::Array{T, 1}, 
                                         bins::Array{Int64, 1}) where {U <:Int, T <: Real}
    w_red .*= 0
    p = sortperm(bins[E[1, :]])
    mem = zeros(Int64, maximum(E))
    ct = i = 1
    while i <= length(p)
        mem .*= 0
        u = bins[E[1, p[i]]]
        while (i <= length(p)) && (bins[E[1, p[i]]] == u)
            if (mem[bins[E[2, p[i]]]] == 0)
                E_red[1, ct] = u
                E_red[2, ct] = bins[E[2, p[i]]]
                mem[bins[E[2, p[i]]]] = ct
                ct += 1
            end
            w_red[mem[bins[E[2, p[i]]]]] += w[p[i]]
            i += 1
        end
    end
    return ct - 1
end

###################################################################################################

function set_bins_f0_red!(bins::Array{U, 1}, f0_red::Array{T, 1}, f::Array{T, 1},
                    g::AbstractGraph{U}) where {U <:Int, T <: Real}
    f0_red .*= 0
    ccmps = strongly_connected_components(g)
    for i = 1:length(ccmps)
        for j = 1:length(ccmps[i])
           bins[ccmps[i][j]] = i
           f0_red[i] += f[ccmps[i][j]]
        end
    end
    return length(ccmps)
end 

###################################################################################################   

"""
```julia
    cutpursuit(d, f0, g, w, cutAlpha, nIter, cut)
```
Procedure for a cut pursuit algorithm with a primal dual approach

# Details
ToDO

# Arguments
This section describes the arguments that can be passed to this function.

## Dimension d
A parameter spcifiying the dimension of the given data. This is due to the 
development process of this project, as the data f0 has to be passed as a 
one dimenionsl array, and hence the information is lost during the reshape 
process. TR: We probably drop this procedure.
- typeof(d) = Int64.

## Data f0
The Original input data. At this stage 
multidimensional input has to be rehaped into this format, which might change 
later.
- typeof(f0) = Array{T, 1} where T <: Real.

## Variational Graph g
This variable stores the underlying graph, the algorithm should be 
performed on. 
- typeof(g) = VariationalGraph{T} where T <: Real.

## Weights w
Weights corresponding to the edges in the variationalGraph g. This is 
currently a bit unclear, as g itself stores weights. It should be 
discussed why g needs to store weights and how the correspond to w.
- typeof(w) = Array{T, 1} where T <: Real.

## cutAlpha
Regularisation parameter to generate the flow graph. 
- typeof(cutAlpha) = T where T <: Real. 

## nIter
The Number of outer iterations for the cut pursuit algorithm.
- typeof(nIter) = Int64

## cut
This parameter allows to dispatch to different choices of cut types for 
the algorithm. Below the possible options are listed.
- typeof(cut) <: cut_type

### aniso()
Perform an anisotropic cut for the cut pusruit algorithm.

### iso()
Currently unsupported

 Syntax:  [result] = cutpursuit(f0, E, w, cutType, cutAlpha, nIter, pdType, pdAlpha)

## Some basic requirements
- length(w) = g.num_edges
- length(f0) = g.num_verts
 
# Output
The result of the algorithm as a Array{T, 1} where T <: Real.

# Examples
ToDO

See also: [`primal_dual`] (@ref),  [`generate_flowgraph`](@ref).
"""
cutpursuit

###################################################################################################

function cutpursuit(d::Int64, f0::Array{T, 1}, g::VariationalGraph, w::Array{T, 1}, 
                    cutAlpha::T, nIter::Int64, cut::cut_aniso) where {T <: Real}
    N = g.num_verts
    # Initialize f with the mean value in each channel
    f = similar(f0)
    for i = 1:d
        mean = sum(f0[((i - 1) * N + 1):(i* N)])/N
        for j = ((i - 1) * N + 1):(i* N)
            f[j] = mean
        end
    end
    E = g.edges_mat
    E_red = similar(E)
    w_red = similar(w)
    f0_red = zeros(N)
    # Initialize bins and binsizes
    bins = ones(Int64, N)
    binsize = N
    # Initialize a directed LightGraph
    G = DiGraph(N)
    for i = 1:g.num_edges
        add_edge!(G, E[1, i], E[2, i])
    end
    
    ## Main Routine
    # Array for cut_inds
    B = falses(N)
    iter = 1
    # Dict for saving history
    history = Dict([("f", zeros(N, d, nIter)), ("bins", zeros(N, nIter)), 
                    ("energy_red", zeros(nIter, 1)), ("energy_full", zeros(nIter, 1))])
    # Treat each dimension independently
    for i = 1:d
        # Object to save cutted graph
        G_cut = G
        # Gathers the edges that were cut
        cutedge = falses(g.num_edges)
        iter = 1
        while iter <= 1 #&& ~isempty(find(cutedge==0,1))
            # Cut graph
            # Go through decoupled problems
            # Generate the corresponding flow-graph
            g_flow, c_matrix, gradJS, Sc = generate_flowgraph(1, f0[((i - 1) * N + 1):(i * N)], 
                                                                w, f[((i - 1) * N + 1):(i * N)], 
                                                                g, cutAlpha, reg_aniso())
            cs, = mincut(g_flow, N + 1, N + 2, c_matrix, BoykovKolmogorovAlgorithm())
            setdiff!(cs, N + 1)
            # Put together for new reduced problem
            # Get the binary partition
            B .*= false
            energy_cut = 0
            for j = 1:length(cs)
                B[cs[j]] = true
            end
            # flag for convergence and value for weight sum
            conv = true
            w_sum = 0
            # Iteration over all edges
            for j = 1:g.num_edges
                if !cutedge[j] && xor(B[E[1, j]], B[E[2, j]]) 
                    conv = false
                    cutedge[j] = true
                    # remove this edge from the graph
                    rem_edge!(G_cut, E[1, j], E[2, j])
                    if Sc[j]
                        w_sum += w[j]
                    end
                end
            end
            #display(collect(edges(G_cut)))
            energy_cut += (1/2) * cutAlpha * w_sum
            # Convergence?
            if conv
                display("Algorithm converged, minimal cut found.")
                break
            end
            #update bins, data and Graph for reduced Problem
            num_verts_red = set_bins_f0_red!(bins, f0_red, f0[((i - 1) * N + 1):(i * N)], G_cut)
            num_edges_red = set_red_graph!(E_red, w_red, E[:, cutedge], w[cutedge], bins)

            #Compute primal dual on reduced problem
            #f = primal_dual(f)
            
            iter += 1
        end
    end
end

###################################################################################################

function cutpursuit(d::Int64, f0::Array{T, 1}, g::VariationalGraph, w::Array{T, 1}, 
                    cutAlpha, nIter::Int64, cut::cut_iso) where {T <: Real}
    error("Currently unsupported!")
end
