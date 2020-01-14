###############################################################################

function generate_flowgraph2(d::Int64, f0::Array{T, 1}, w0::Array{T, 1}, f::Array{T, 1},
                                  g::VariationalGraph, alpha::T, reg::reg_aniso) where {T <: Real}
    N = d * g.num_verts
    #show(N)
    gradJS = zeros(d * g.num_verts)
    Sc = falses(d * g.num_edges)
    flowgraph = DiGraph(N + 2)
    capacity_matrix = spzeros(N + 2, N + 2)
    ect = 1
    #why not sqrt(w0[u])
    for i = 1:d
        off  = (i - 1) * g.num_verts
        for u = 1:g.num_verts
            acc = 0
            for v = 1:length(g.edges_list[u])
                tmp = f[u + off] - f[g.edges_list[u][v] + off]
                if abs(tmp) >= 1e-10
                    acc += sign(tmp) * w0[u]
                else
                    add_edge!(flowgraph, u, g.edges_list[u][v])
                    capacity_matrix[u + off, g.edges_list[u][v] + off] =
                         alpha * w0[u]# eventuel ohne (1/2)* sqrt
                    Sc[ect] = true
                end
                ect += 1
            end
            grad = f[u + off] - f0[u + off] + alpha * acc
            gradJS[u] = grad
            if grad >= 0
                add_edge!(flowgraph, N + 1, u + off)
                capacity_matrix[N + 1, u + off] = grad
            else
                add_edge!(flowgraph, u + off, N + 2)
                capacity_matrix[u + off, N + 2] = -grad
            end
        end
    end
    return flowgraph, capacity_matrix, gradJS, Sc
end

###############################################################################

"""
```julia
    generate_flowgraph(d, f0, w0, f, g, alpha, reg)
```
This function generates the flow graph for anisotropic and isotropic TV
regularisation as described in the paper.

# Details
ToDO

# Arguments
This section describes the arguments that can be passed to this function.

## Data f0
The initial input data as Array{T, 1} where T <: Real. At this stage
multidimensional input has to be reshaped into this format, which might change
later.
- typeof(f0) = Array{T, 1} where T <: Real

## Weights w0
The initial weights.
- typeof(w0) = Array{T, 1} where T <: Real

## Data f
Current iterated data f^k as . At this stage
multidimensional input has to be reshaped into this format, which might change
later.
- typeof(f) = Array{T, 1} where T <: Real

## Variational Graph g
This variable stores the underlying graph, the algorithm should be
performed on.
- typeof(g) =  VariationalGraph{T} where T <: Real.

## alpha
Regularization parameter.
- typeof(alpha) = T where T <: Real.

## reg
This parameter allows to dispatch to different choices of regularization types for
the algorithm. The possible options are listed below.
- typeof(reg) <: reg_type

## reg_aniso()
ToDO

## reg_iso
Curently unsupported.

## Some basic requirements
- length(w0) = g.num_edges
- length(f0) = length(f) = g.num_verts

# Output
This section describes the output of the algorithm.

## flowgraph
The resulting graph of the algorithm.
- typeof(flowgraph) = LightGraphs.DiGraph{U} where {U <: Int}

## capacity_matrix
Sparse matrix that storse the capacities corrensponding to the edges of
the flowgraph.
- typeof(capacity_matrix) = SparseMatrixCSC{U, T} where {U <: Int, T <: Real}
## gradJS
ToDO
- typeof(gradJS) = Array{T, 1} where T <: Real

## Sc
Array that stores wheter a corresponding edge produces a zero value.
- typeof(Sc) = Array{Bool, 1}

## Some basic properties
- length(Sc) = g.num_edges = nnz(capacity_matrix)
- length(gradJS) = g.num_verts

# Examples
ToDo
"""
generate_flowgraph2
