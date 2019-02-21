###############################################################################

function generate_flowgraph_aniso(d::Int64, f0::Array{T, 1}, w0::Array{U, 1}, f::Array{T, 1}, g::VariationalGraph, alpha::T) where {T <: Real, U <: Real}
    N = length(f)
    flowgraph = DiGraph(N + 2)
    capacity_matrix = spzeros(N + 2, N + 2)
    for i = 1:d
        off = (i - 1) * g.num_verts
        for u = 1:g.num_verts
            acc = 0
            for v = 1:length(g.edges[u])
                tmp = abs(f[u + off] - f[v + off])
                if tmp >= 1e-10
                    acc += sign(tmp) * w0[u]
                else
                    add_edge!(flowgraph, u + off, v + off)
                    capacity_matrix[u + off, v + off] = 0.5 * alpha * w0[u]
                end
            end
            gradJS = f[u + off] - f0[u + off] + alpha * acc 
            if gradJS >= 0
                add_edge!(flowgraph, N + 1, u + off)
                capacity_matrix[N + 1, u + off] = gradJS
            else
                add_edge!(flowgraph, u + off, N + 2)
                capacity_matrix[u + off, N + 2] = gradJS
            end
        end
    end
    return flowgraph, capacity_matrix
end

###############################################################################

function generate_flowgraph_iso(f0::Array{T, 1}, w0::Array{U, 1}, f::Array{T, 1}, g::VariationalGraph, alpha) where {T <: Real, U <: Real}
    (N,d) = size(f);
#         gradRS = zeros(size(f));
#         
#         % Compute the 2-norm of (f(u)-f(v))
#         R = sqrt(sum((f(E(:,1),:)-f(E(:,2),:)).^2,2));
#         R(abs(R)<1e-10) = 0;
#         
#         % Get the set of diffable parts of the regularizer S
#         S = R~=0;
#         
#         % Compute the gradient of the diffable parts of R
#         for i=1:d
#             gradRS(:,i) = accumarray(E(S,1),w0(S).*(f(E(S,1),i)-f(E(S,2),i))./R(S),[N,1]);
#         end
#         
#         % Gradient of diffable parts of J
#         gradJS = (f-f0) + alpha*gradRS;
#         
#         % Add the nodes to the graph
#         G = digraph;
#         G = G.addnode(N);
#         
#         % Add sink and source to graph
#         G = G.addnode('s'); % Position: N+1
#         G = G.addnode('t'); % Position: N+2
#         
#         % Compute the sum of the gradient of JS over the dimensions
#         sumgradJS = sum(gradJS,2);
#         
#         % Find the positive and negative entries of gradient J_S
#         ind_pos = find(sumgradJS>0);
#         ind_neg = find(sumgradJS<=0);
#         
#         % Add the edges as written in paper 
#         G = addedge(G,N+1,ind_pos,sumgradJS(ind_pos));
#         G = addedge(G,ind_neg,N+2, -sumgradJS(ind_neg));
#         
#         % Add edges between the nodes where R is non-diffable
#         if ~isempty(~S)
#             G = addedge(G,E(~S,1),E(~S,2),1/2*alpha*sqrt(d*w0(~S)));
#         end
end

###############################################################################

"""
    generate_flowgraph
    
This function generates the flow graph for anisotropic and isotropic TV
regularisation as described in the paper

Syntax:  G = GenerateFlowGraph(f0,w0,f,E,alpha,type)

# Inputs:
- f0      Initial input data 
- w0      Initial weights 
- f       Current iterated data f^k
- E       Given edgeset 
- alpha   Regularization parameter
- type    Type of regularization R ['iso', 'aniso'];
"""
function generate_flowgraph(d::Int64, f0::Array{T, 1}, w0::Array{U, 1}, f::Array{T, 1}, 
g::VariationalGraph, alpha::T, type_R::String) where {T <: Real, U <: Real}
    # Generate the flow graph for the anisoropic or isotropic regularizer R
    if type_R == "aniso"
        generate_flowgraph_aniso(d, f0, w0, f, g, alpha)
    elseif type_R == "iso"
        error("Currently unsupported!")
        generate_flowgraph_iso(d, f0, w0, f, g, alpha)
    else
        error("Unsupported type for regularizer")
    end
end