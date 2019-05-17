using DelimitedFiles
using VariationalGraphs
using LightGraphs
using PyPlot
using Printf
printstyled(@sprintf("Finished Loading\n"); color =:reverse)
# Loading +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set data
I = readdlm("../data/cameraman.txt", Float64)
n, m = size(I)
f = reshape(I, n * m)
g = VariationalGraph(I)
# Set PD Configuration ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Store parameters in img_params struct
struct img_params{T<:Real, U<:Integer} <: pd_params{T, U} 
    sigma::T
    tau::T
    theta::T
    threshold::T
    #
    maxiter::U
end
# Set parameters for pd algorithm
par = img_params(0.35, 0.35, 1.0, 1e-07, 4)
#prox operator for FStar
function prox_Fstar(y::Array{T, 1}, par::img_params) where {T<:Real}
    return y ./ maximum(abs(y))
end
#prox operator for G
lambda = 100
function prox_G(x::Array{T, 1}, par::img_params) where {T<:Real}
    return (x + par.tau * f ./ lambda) ./ (1 + par.tau/lambda)
end
# Operator for the problem
function pd_op(x::Array{T, 1}, par::img_params) where {T<:Real}
    return gradient(x, g, g.config)
end
# Transposed Operator for the problem
function pd_trans_op(y::Array{T, 1}, par::img_params) where {T<:Real}
    return y
end
# Apply norm to x
function compute_norm(x::Array{T, 1}, par::img_params) where {T<:Real}
    sqrt(sum((x) .^ 2))
end
# Apply PD algorithm ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
x = zeros(size(f))
y = pd_op(x, par)
u = primaldual(x, y, par);




# Visualization +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# function pos(n::Int64, m::Int64, u::Int64)
#     y = n - div(u - 1, m)
#     x = 1 + mod(u - 1, m) 
#     return (x, y)
# end
# n, m = size(I)
# maxi = maximum(I)
# for u = 1:g.num_verts
#     i = 1 + div(u - 1, m)
#     j = 1 + mod(u - 1, m)
#     c = string(I[i,j]/maxi)
#     (x,y) = pos(n, m, u)
#     plot(x, y, color = c, marker =:s)
#     #annotate(string(u), (points[1,u], points[2,u]))
#     for i = 1:length(g.edges_list[u])
#         v = g.edges_list[u][i]
#         (x2, y2) = pos(n, m, v)
#         plot((x, x2), (y, y2), color = "b")
#     end
# end
# Visualization +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

