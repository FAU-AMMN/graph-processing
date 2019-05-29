abstract type pd_params{T<:Real, U<:Integer} end

###############################################################################
# The following setting of paramter struct and operators is used for the case
# of an image-Graph and is hence dispatched by the strut img_params
struct img_params{T<:Real, U<:Integer} <: pd_params{T, U} 
    sigma::T
    tau::T
    theta::T
    threshold::T
    lambda::T
    #
    maxiter::U
    #
    f::Array{T, 1}
    #
    g::VariationalGraph{U, T}
end
#prox operator for FStar
function prox_Fstar(y::Array{T, 1}, par::img_params) where {T<:Real}
    return y ./ max.(abs.(y), 1)
end
#prox operator for G
function prox_G(x::Array{T, 1}, par::img_params) where {T<:Real}
    return (x + par.tau * par.f ./ par.lambda) ./ (1 + par.tau/par.lambda)
end
# Operator for the problem
function pd_op(x::Array{T, 1}, par::img_params) where {T<:Real}
    return VariationalGraphs.gradient(x, par.g, par.g.config)
end
# Transposed Operator for the problem
function pd_trans_op(y::Array{T, 1}, par::img_params) where {T<:Real}
    return VariationalGraphs.divergence(y, par.g, par.g.config)
end
# Apply norm to x
function compute_norm(x::Array{T, 1}, par::img_params) where {T<:Real}
    sqrt(sum((x) .^ 2))
end
###############################################################################



###############################################################################
"""
```julia
    primaldual(f::Array{T, 1}, par::pd_params) where T <: Real
```
ToDo: Write DocString!
"""
function primaldual(x::Array{T, 1}, y::Array{T, 1}, par::P) where{T<:Real, U<:Integer, P<:pd_params{T, U}}
    history = Dict{Integer, Array{T, 1}}()
    x_old = x; x_bar = x
    iter = 1
    rel_change = par.threshold + 1
    while (iter <= par.maxiter) && (rel_change > par.threshold)
        y = prox_Fstar(y + par.sigma * pd_op(x_bar, par), par)
        x = prox_G(x + par.tau * pd_trans_op(y, par), par)
        history[iter] = x
        x_bar = x + par.theta * (x - x_old)
        x_old = x
        rel_change = compute_norm(x - x_old, par) / compute_norm(x, par)
        iter += 1
    end
    return x, history
end

###############################################################################
