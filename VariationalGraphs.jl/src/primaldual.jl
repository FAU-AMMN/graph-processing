abstract type pd_params{T<:Real, U<:Integer} end

###############################################################################

"""
```julia
    primaldual(f::Array{T, 1}, par::pd_params) where T <: Real
```
ToDo: Write DocString!
"""
function primaldual(x::Array{T, 1}, y::Array{T, 1}, par::P) where{T<:Real, U<:Integer, P<:pd_params{T, U}}
    x_old = x; x_bar = x
    iter = 1
    rel_change = 0
    while iter <= par.maxiter && rel_change > par.threshold
        y = prox_Fstar(y + sigma * pd_op(x_bar, par), par)
        x = prox_G(x + par.tau * pd_trans_op(y, par), par)
        x_bar = x + par.theta * (x - x_old)
        x_old = x
        rel_change = compute_norm(x - x_old, par) / compute_norm(x, par)
    end
    return x
end

###############################################################################
