###############################################################################

struct pd_params{T <:Real, U<:Int}
    gamma::T
    delta::T
    epsilon::T
    theta::T
    normK::T
    preconConstant::T
    tol::T
    
    p::Float64
    q::Float64
    
    niter::U
    niter_power::U
    
    verbose::Int16
    setSize::Int16
    # Acceleration parameter
    # 2: StepsizeUpdate
    # 1: Preconditioning
    # else: None 
    acceleration::Int16
    
    powerit::Bool
    
    reg::reg_type
end

###############################################################################

function powerit_gradient(E, w, par.niter_power)

end

###############################################################################

function compute_primal(f, f0, E, w, alpha::T, epsilon::T, q::Float64, p::Float64, reg::reg_type) where {T <:Real}

end

###############################################################################

function compute_dual(y, divy, f0, alpha::T, epsilon::T) where {T <:Real}

end

###############################################################################

"""
```julia
    primal_dual(f::Array{T, 1}, par::pd_params) where T <: Real
```
ToDo: Write DocString!
"""
function primal_dual(d::Int64, f::Array{T, 1}, w, E, params::pd_params) where {T <: Real}
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
    if par.powerit == true && ~strcmp(par.acceleration,"Preconditioning")
        if P.nIter >= 1
            L = powerit_gradient(E, w, par.niter_power)
            if par.verbose >= 2
                #fprintf('\t %s \n', ['Estimation of operator norm via power iterations is ' num2str(L)]);
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
    y = grad(f,E,w);
    divy = div(y,E,w,m);
    # Set the hubert update parameter h 
    if par.epsilon > 0
        h = 1 - par.delta
    else 
        h = 1
    end
    #
    primal = compute_primal(f, f0, E, w, par.alpha, par.epsilon, par.reg, par.q, par.p)
    dual = compute_dual(y, divy, f0, par.alpha, par.epsilon)
    #
    gap = primal - dual
    primal_old = primal
    #
    while iter <= par.niter && err > par.tol
        #----\TR{
        #It seems to me, that the updat of f should be after the y update,
        # I should check this!
        #----}
        # Set f^k 
        f_old = f
        # Primal Update f^k+1
        f = 1./(1 + tau .* 1) .*(f + tau.*(f0 + divy))
        # Dual Update y^k+1
        y = proj(h * y + sigma .* grad(f + theta * (f - f_old), E, w), E, par.alpha, par.q, par.p, par.reg)
        # Update tau, simga and theta
        if par.acceleration == 2
            theta = 1/(sqrt(1+2*gamma*tau));
            tau = theta*tau;
            sigma = sigma/theta;
        end
        # Compute divergence of y to save calls
        divy = div(y,E,w,m);
        # Compute primal-dual gap
        if mod(iter, P.update_interval) == 0
            primal = ComputePrimal(f,f0,E,w,P.alpha,P.epsilon,P.regType,P.q,P.p);
            dual = ComputeDual(y,divy,f0,P.alpha,P.epsilon); 
            #
            if strcmp(par.errType, "energy") || strcmp(par.errType, "indifferent_energy")
                [err,s] = computeError(primal, primal_old, m, par.errType)
                primal_old = primal;
            else
                gap_old = gap
                gap  = primal - dual
                [err,s] = computeError(gap,gap_old,m,P.errType)
            end
            
            if P.verbose >= 1
                #output = ['  Iteration %', num2str(numel(num2str(P.nIter))),...
                #'.0f/%.0f', s];
                #fprintf(output, iter, P.nIter);
                #drawnow;
            end
        end
        iter += 1
    end
    
    return f
end
###############################################################################