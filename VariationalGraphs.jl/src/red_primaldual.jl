abstract type pd_params{T<:Real, U<:Integer} end
abstract type pd_history{T<:Real, U<:Integer} end

###############################################################################
# The following setting of paramter struct and operators is used for the case
# of an reduced graph of the cutpursuit-algorithm and is hence dispatched by
# the strut red_params
struct red_params{T<:Real, U<:Integer} <: pd_params{T, U}
   # parameters of the primal-dual-algorithm
   sigma::T
   tau::T
   theta::T
   threshold::T
   # regularization parameter
   lambda::T
   #
   maxiter::U
   history_step::U
   # the reduced vertex function
   f::Array{T, 1}
   # number of vertices of each component of the reduced graph
   setSize::Array{T, 1}
   # the reduced graph
   g::VariationalGraph{U, T}
end
# history for the image type
struct red_history{T<:Real, U<:Integer} <: pd_history{T, U}
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
function red_history(maxiter::U, history_step::U, x::Array{T, 1}, y::Array{T, 1}) where{T<:Real, U<:Integer}
   size = div(maxiter, history_step) + 1
   X = zeros(size, length(x)); X[1, :] = x;
   X2 = zeros(size, length(x)); X[1, :] = x;
   Y = zeros(size, length(y)); Y[1, :] = y;
   return img_history(X, X2, Y, zeros(size), zeros(size), size, [1])
end
# Add a iteration step to the history
function red_histiory_add!(hist::red_history, rel_change::T, energy::T,
                         x::Array{T, 1}, x_bar::Array{T, 1}, y::Array{T, 1}) where{T<:Real, U<:Integer}
   iter =  hist.iter[1] + 1
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
#function prox_Fstar(y::Array{T, 1}, par::red_params) where {T<:Real}
#   return y ./ max.(abs.(y), 1)
#end

#prox operator for G
#function prox_G(x::Array{T, 1}, par::red_params) where {T<:Real}
#   return (x + par.tau * par.f ./ par.lambda) ./ (1 + par.tau/par.lambda)
#end

# Operator for the problem
function pd_op(x::Array{T, 1}, par::red_params) where {T<:Real}
   return VariationalGraphs.gradient(x, par.g, par.g.config)
end
# Transposed Operator for the problem
function pd_trans_op(y::Array{T, 1}, par::red_params) where {T<:Real}
   return VariationalGraphs.divergence(y, par.g, par.g.config)
end
# Apply norm to x
function compute_norm(x::Array{T, 1}, par::red_params) where {T<:Real}
   sqrt(sum((x) .^ 2))
end




# projection onto the scaled unit ball in the dual norm
function proj(z::Array{T, 1}, par::red_params) where {T<:Real}
   return par.lambda*z ./ max.(abs.(z), par.lambda)
   #d = size(z,2);
   #N = max(E(:,1));

   #% Each element in absolute value should be smaller or equal alpha
   #nz = ComputeNorm(z,N,E,type,1,1);
   #z = alpha*z./(max(alpha,reshape(nz,[],d)));
end
###############################################################################



###############################################################################
"""
```julia
   red_primaldual(f::Array{T, 1}, par::pd_params) where T <: Real
```
Procedure for the reduced primal dual algorithe for the cut-pursuit-algorithem
from the paper

# Details
ToDO

# Arguments
This section describes the arguments that can be passed to this function.

## Parameters par
Struct containing the parameters for the reduced primal-dual-algorithm.
- typeof(par) = red_params{T<:Real, U<:Integer} <: pd_params{T, U}.

ToDo: Write DocString!
"""
function red_primaldual( par::P) where{T<:Real, U<:Integer, P<:pd_params{T, U}}
   # initializes
   iter = 1
   rel_change = par.threshold + 1
   f=par.f
   f0=par.f
   y =pd_op(f,par) #grad(f,E,w);
   divy = pd_trans_op(y,par)#div(y,E,w,m);
   # Initilaize history.
   history = red_history(par.maxiter, par.history_step, f, y)
   while (iter <= par.maxiter) && (rel_change > par.threshold)#iter <= P.nIter && err > P.tol
       # Set f^k
       f_old = f
       # - Primal-Dual algorithm
       # Primal Update f^k+1
       f = 1 ./(1 .+ par.tau .*par.setSize).*(f + par.tau .*(f0 + divy))
       # Dual Update y^k+1
       y = proj(y+par.sigma .*pd_op(f+par.theta*(f-f_old),par),par)

       divy = pd_trans_op(y,par)

       rel_change = compute_norm(f - f_old, par) / compute_norm(f, par)
       #history update (not implemented)
       if iter%par.history_step == 0
          # img_histiory_add!(history, rel_change,
         #                    0.5 * sum((x - par.f) .^ 2) + par.lambda * sqrt(sum(pd_op(x, par) .^ 2)),
         #                    x, x_bar, y)
       end
       iter += 1
    end
    return f, history
end

###############################################################################
