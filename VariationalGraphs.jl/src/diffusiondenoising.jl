using Printf
function diffusion(steps::T,tau::T,f::Array{T, 2},conf::GraphConfig) where{T<:Float64}
    #printstyled(@sprintf("Hallo\n"); color =:reverse)
    n,m=size(f)
    g=VariationalGraph(f, conf)

    w=zeros(n,m)        #This is a workaround
    w=w+f'              #for a reshape problem

    f=reshape(w, n*m)
    for k=1:steps
        lap=Laplace(f,g,nonlocal())  # nonlocal()
        f=f+tau*lap;
        g=VariationalGraph(reshape(f,n,m), conf)
    end

    w=zeros(n,m)        #This is a workaround
    w=w+reshape(f,n,m)' #for a reshape problem

    return w
end
