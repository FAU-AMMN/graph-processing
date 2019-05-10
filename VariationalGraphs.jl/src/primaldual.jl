###############################################################################

struct pd_params{T<:Real, U<:Int64}
    gamma::T
    delta::T
    epsilon::T
    theta::T
    normK::T
    preconConstant::T
    tol::T
    p::T
    q::T
    #
    niter::U
    niter_power::U
    verbose::U
    setSize::U
    update_interval::U
    # Acceleration parameter
    # 2: StepsizeUpdate
    # 1: Preconditioning
    # else: None 
    acceleration::U
    # Error Type
    # 1: energy
    # 2: indifferent energy
    # 3: gap
    # 4: RMSE
    # 5: indifferent gap
    error::U
    #
    powerit::Bool
    #
    reg::reg_type
end

###############################################################################
"""
```julia
    compute_norm()
```
ToDo: Write DocString!
"""
function compute_norm(d::Int64, m::Int64, y::Array{T, 2}, E::Array{Int64, 2}, q::T, p::T, norm::reg_type) where{T<:Real}
# This function computes the different norms for different types of
# minimization problems.
#
# Syntax:  normY = ComputeNorm(y,N,E,normType)
#
# Inputs:
#    y          - Input data of the size of a gradient (length(E)xd)
#    N          - Number of vertices 
#    E          - Edge set corresponding to N
#    normType   - Type of the minimization problem. 
#                   ['aniso-aniso','aniso-iso','iso-aniso','iso-iso'] 
#   q           - The q of the pq-norm
#   p           - The p of the pq-norm
#
# Outputs:
#    normY     - Norm of the input y
#
if q == 1 # Every norm is the same as the aniso-aniso norm!
    normY = abs(y[:])
else 
    normY = compute_norm(d, m, q, y, E, norm)
end
# When p is infinity the pq-norm is the max over the q norm
if p == Inf
    normY = max(normY[:])
elseif p != 1
    normY = sum(normY[:] .^ p) .^ (1/p)
end
return normY
# If p is 1 then the pq-norm becomes a sum over inner q-norms
# TR: Is sth missing here??
end
# Dispatch the different reg types
#----------------------------------------------------------------------------------------
function compute_norm(d::U, m::U, q::Float64, y::Array{Float64, 2}, E::Array{U, 2}, norm::reg_aniso_aniso) where U <: Int64
    abs(y[:])
end
#----------------------------------------------------------------------------------------
function compute_norm(d::U, m::U, q::Float64, y::Array{Float64, 2}, E::Array{U, 2}, norm::reg_aniso_iso) where U <: Int64
    if q == Inf
        return maximum(abs(y), dims = 1)
    else
        # lp-norm of each dimension in y
        return (sum(abs(y).^q, dims = 1)) .^ (1/q)
    end
end
#----------------------------------------------------------------------------------------
function compute_norm(d::U, m::U, q::Float64, y::Array{Float64, 2}, E::Array{U, 2}, norm::reg_iso_aniso) where U <: Int64
    # Maximum of each set of neighbors
    normY = zeros(d * m)
    if q == Inf
        for l = 1:d
            for i = 1:m
                normY[E[1, i] + (l - 1) * m] = max(normY[i + (l - 1) * m], abs(y[i + (l - 1) * m]))
            end
        end
    else
        for l = 1:d
            for i = 1:m
                normY[E[1, i] + (l - 1) * m] += abs(y[i + (l - 1) * m]) ^ q
            end
        end
        normy = normY .^ (1/q)
    end
    return normY
end
#----------------------------------------------------------------------------------------
function compute_norm(d::U, m::U, q::Float64, y::Array{Float64, 2}, E::Array{U, 2}, norm::reg_iso_iso) where U <: Int64
    # Maximum over each set of neigbors in all dimensions
    normY = zeros(m)  
    if q == Inf
        for l = 1:d
            for i = 1:m
                normY[E[1, i]] = max(normY[i], abs(y[i + (l - 1) * m]))
            end
        end
    else
        for l = 1:d
            for i = 1:m
                normY[E[1, i]] += abs(y[i + (l - 1) * m]) ^ q
            end
        end
        normy = normY .^ (1/q)
    end
end

###############################################################################
"""
```julia
    proj()
```
Projection onto C
ToDo: Write DocString!
"""
function proj(d::U, m::U, y::Array{T, 2}, E::Array{U, 2}, alpha::T, q::T, p::T, reg::reg_type) where{T<:Real, U<:Int64}
    N = max(E[1, :])
    # Same distinction as in the paper
    if p == 1 && q == 1
        # Each element in absolute value should be smaller or equal alpha
        nz = compute_norm(d, m, y, E, 1, 1, reg)
        proj!(d, m, y, nz, E, alpha, reg_aniso_aniso())
        return y
    elseif p == 1 && q ~= 1
        nz = compute_norm(d, m, y, E, q/(q-1), 1, reg)
        proj!(d, m, y, nz, E, alpha, reg)
        
#             
#         elseif strcmp(type,'aniso-iso') % Norm of each row  <= alpha
#             z = alpha*z./repmat((max(alpha,nz)),1,d);
#         elseif strcmp(type,'iso-aniso') % Norm of a set of neighbors <= alpha
#             nz = reshape(nz,[],d);
#             z =  alpha*z./max(alpha,nz(E(:,1),:));
#         elseif strcmp(type,'iso-iso') % Norm of set of neigbors over all dimensions <= alpha
#             z = alpha*z./max(alpha, repmat(nz(E(:,1)),1,d));
#         end
    elseif p != 1 && q != 1
        # The projection projects the whole array into an alpha p*q*-ball
        #nz = ComputeNorm(z,N,E,type,q/(q-1),p/(p-1));
        z = alpha*z/max(alpha,nz);
    end
end
#----------------------------------------------------------------------------------------
function proj!(d::U, m::U, y::Array{T, 2}, nz::Array{T, 1}, E::Array{U, 2}, alpha::T, norm::reg_aniso_aniso) where{T<:Real, U<:Int64}
    for l = 1:d
        for i = 1:m
            y[l, i] = alpha * y[l, i] / max(alpha, nz[i])
        end
    end
end
#----------------------------------------------------------------------------------------
# function proj!(d::U, m::U, y::Array{T, 2}, nz::Array{T, 1}, E::Array{U, 2}, alpha::T, norm::reg_aniso_aniso) where{T<:Real, U<:Int64}
#     for l = 1:d
#         for i = 1:m
#             y[l, i] = alpha * y[l, i] / max(alpha, nz[i])
#         end
#         z = alpha*z./repmat((max(alpha,nz)),1,d);
#     end
# end

###############################################################################


function huber(nz, epsilon)

end

###############################################################################

function powerit_gradient(E, w, niter)

end

###############################################################################

function compute_primal(f::Array{T, 1}, f0::Array{T, 1}, E::Array{Int64, 1}, 
                        w::Array{T, 1}, N::Int64, alpha::T, epsilon::T, q::T, p::T, reg::reg_type) where {T <:Real}
    D  = 1/2 * sum(f - f0 .^ 2)
    gradf = grad(f, E, w)
    nz = compute_norm(gradf, N, E, q, p, reg)
    # Compute huber if acitve
    if epsilon > 0
        nz = huber(nz, epsilon)
    end
    # Compute the regularizer
    R = sum(nz)
    primal = D + alpha * R
end

###############################################################################

function compute_dual(d::Int64, y::Array{T, 2}, divy::Array{T, 1}, f0, alpha::T, epsilon::T) where {T <:Real}
    dual = - sum(divy .* f0) - 0.5 * sum(divy .^ 2)
    if epsilon > 0
        dual = dual + epsilon/(2 * alpha) * sum(y .^ 2);
    end
end

###############################################################################

function compute_error(gap::T, gap_old::T, m::Int64, err::Int16) where {T <: Real}
    if err == 3
        err = gap
        s = @sprintf(", Primal-Dual-Gap: %12.8f", err)
    elseif err == 4
        err = sqrt(2 * gap/m)
        s = @sprintf(", RMSE: %12.8f\n", err)
    elseif err == 5
        err = abs(gap_old-gap)
        s = @sprintf(", Error between gaps: %12.8f", err)
    elseif err == 1
        err = gap;
        s = @sprintf(", Energy: %12.8f", err)
    elseif err == 2 
        err = abs(gap_old-gap)
        s = @sprintf(", Gap between energies: %012.8f", err)
    end
    return err, s
end

###############################################################################

"""
```julia
    primal_dual(f::Array{T, 1}, par::pd_params) where T <: Real
```
ToDo: Write DocString!
"""
function primal_dual(d::Int64, f::Array{T, 1}, g::VariationalGraph, par::pd_params) where {T <: Real}
    m = div(length(f), d)
    # Original Data
    f0 = f
    iter = 1
    # Set values
    theta = par.theta # Relaxation
    gamma = 0 # Strong convexity for stepsize update
    # If huber is active compute strong convexity constant of regularizer
    if par.epsilon > 0 
        par.delta = par.epsilon/par.alpha
    else 
        par.delta = 0.99
    end
    #
    L = par.normK;
    # norm of K (HAS TO BE APPROXIMATED)
    if par.powerit == true && par.acceleration == 1
        if par.nIter >= 1
            L = powerit_gradient(g, par.niter_power)
            if par.verbose >= 2
                #@sprintf('\t %s \n', ['Estimation of operator norm via power iterations is ' num2str(L)]);
            end
        end
    end 
    #Setup stepsize
    if par.acceleration == 2
        gamma = par.gamma
        tau = sqrt(par.delta)/L
        sigma = 1/(sqrt(par.delta)*L)
    elseif par.acceleration == 1
        #a = par.preconConstant
        #tau   = 1 ./(ColumnSum(w, E, m).^ (2 - a) + eps())
        #sigma = 1 ./(RowSum(w) .^a)
        @printf("not supported yet :(")
    else
        # Just set tau and sigma 
        tau = sqrt(par.delta)/L
        sigma = 1/(sqrt(par.delta) * L)
    end
    # Set some value for gap
    err = 1000
    # Initialize y
    y = grad(f, g)
    divy = diveregnce(d, m, y, g)
    # Set the hubert update parameter h 
    if par.epsilon > 0
        h = 1 - par.delta
    else 
        h = 1
    end
    #
    primal = compute_primal(f, f0, g, m, par.alpha, par.epsilon, par.reg, par.q, par.p)
    dual = compute_dual(y, divy, f0, par.alpha, par.epsilon)
    #
    gap = primal - dual
    primal_old = primal
    #------------------------------------------------------------------------------------
    while iter <= par.niter && err > par.tol
        #----\TR{
        #It seems to me, that the updat of f should be after the y update,
        # I should check this!
        #----}
        # Set f^k 
        f_old = f
        # Primal Update f^k+1
        f = 1 ./(1 + tau .* 1) .*(f + tau .* (f0 + divy))
        # Dual Update y^k + 1
        y = proj(h * y + sigma .* grad(f + theta * (f - f_old), E, w), E, par.alpha, par.q, par.p, par.reg)
        # Update tau, simga and theta
        if par.acceleration == 2
            theta = 1 / (sqrt(1 + 2 * gamma * tau))
            tau = theta * tau
            sigma = sigma / theta
        end
        # Compute divergence of y to save calls
        divy = divergence(d, m, y, E, w);
        #--------------------------------------------------------------------------------
        # Compute primal-dual gap
        if mod(iter, par.update_interval) == 0
            primal = compute_primal(f, f0, g, m, par.alpha, par.epsilon, par.q, par.p, par.reg);
            dual = compute_dual(y,divy,f0,par.alpha,par.epsilon); 
            #
            if par.error == 1 || par.error == 2 
                err, s = compute_error(primal, primal_old, m, par.error)
                primal_old = primal
            else
                gap_old = gap
                gap  = primal - dual
                err, s = compute_error(gap, gap_old, m, par.error)
            end
            #----------------------------------------------------------------------------
            if par.verbose >= 1
                printstyled(@sprintf("Finished iteration %d of at most %d iterations", iter, par.niter); color =:reverse)
            end
        end
        #--------------------------------------------------------------------------------
        iter += 1
    end
    #------------------------------------------------------------------------------------
    return f
end
###############################################################################