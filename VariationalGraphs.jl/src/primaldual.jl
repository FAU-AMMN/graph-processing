abstract type pd_params{T<:Real, U<:Integer} end
abstract type pd_history{T<:Real, U<:Integer} end

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
# history for the image type
struct img_history{T<:Real, U<:Integer} <: pd_history{T, U}
    x::Array{T, 2}
    x_bar::Array{T, 2}
    y::Array{T, 2}
    rel_change::Array{T, 1}
    energy::Array{T, 1}
    #
    maxiter::U
    iter::Vector{U}
end
# Constructor for initialization
function img_history(maxiter::U, x::Array{T, 1}, y::Array{T, 1}) where{T<:Real, U<:Integer}
    X = zeros(maxiter + 1, length(x)); X[1, :] = x;
    X2 = zeros(maxiter + 1, length(x)); X[1, :] = x;
    Y = zeros(maxiter + 1, length(y)); Y[1, :] = y;
    return img_history(X, X2, Y, zeros(maxiter + 1), zeros(maxiter + 1), maxiter, [1])
end
# Add a iteration step to the history
function img_histiory_add!(hist::img_history, iter::U, rel_change::T, energy::T, 
                          x::Array{T, 1}, x_bar::Array{T, 1}, y::Array{T, 1}) where{T<:Real, U<:Integer}
    if iter <= hist.maxiter
        hist.x[iter, :] = x
        hist.x_bar[iter, :] = x_bar
        hist.y[iter, :] = y
        hist.iter[1] = iter
        hist.rel_change[iter] = rel_change
        hist.energy[iter] = energy
    end
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
    # Initilaize history
    history = img_history(par.maxiter, x, y)
    # Initialization 
    x_old = x; x_bar = x
    iter = 1
    rel_change = par.threshold + 1
    while (iter <= par.maxiter) && (rel_change > par.threshold)
        y = prox_Fstar(y + par.sigma * pd_op(x_bar, par), par)
        x = prox_G(x + par.tau * pd_trans_op(y, par), par)
        x_bar = x + par.theta * (x - x_old)
        rel_change = compute_norm(x - x_old, par) / compute_norm(x, par)
        # Update iter and history
        iter += 1
        img_histiory_add!(history, iter, rel_change, 
                          0.5 * sum((x - par.f) .^ 2) + par.lambda * sum(pd_op(x, par) .^ 2), 
                          x, x_bar, y)
        # Save old x TR: We could use hist here
        x_old = x
    end
    return x, history
end

###############################################################################
