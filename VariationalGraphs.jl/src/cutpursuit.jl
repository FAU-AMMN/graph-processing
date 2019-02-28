###############################################################################

function set_red_graph!(E_red::Array{T, 2}, w_red::Array{T, 1},
                                         ct::Int, E::Array{T, 2}, w::Array{T, 1}, 
                                         bins::Array{Int64, 1}) where T
    p = sortperm(bins[E[1,:]])
    mem = zeros(Int64, maximum(E))
    ct = 0
    while i <= length(p)
        mem .*= 0
        u = bins[E[1, p[i]]]
        while (bins[E[1, p[i]]] == u)
            if (mem[bins[E[2, p[i]]]] == 0)
                E_red[1, ct] = u
                E_red[2, ct] = bins[E[2, p[i]]]
                mem[bins[E[2, p[i]]]] = ct
                ct += 1
            end
            w_red[mem[bins[E[2, p[i]]]]] += w[p[i]]
            i += 1
        end
    end
end

###############################################################################

function set_bins!(bins::Array{Int64, 1}, g::AbstractGraph{T}) where T
    ccmps = strongly_connected_components(g)
    for i = 1:length(ccmps)
        for j = 1:length[ccmps[i]]
           bins[ccmps[i][j]] = i
        end
    end
end 

###############################################################################   

"""
CutPursuit - Procedure for a cut pursuit algorithm with a primal dual approach

 Syntax:  [result] = cutpursuit(f0, E, w, cutType, cutAlpha, nIter, pdType, pdAlpha)

# Inputs:
    - f0          Original input data as a Nxd vector
    - E           Input edge set corresponding to the data f0
    - w           Weights corresponding to the edges in E
    - cutType     Type of the flowgraph. ['iso', 'aniso'] (Also the default
                 for pdType)
    - cutAlpha    Regularisation parameter to generate the flow graph 
                   (Also the default for pdAlpha)
    - nIter       Number of outer iterations for the cut pursuit
    - p           The p of the qp-norm
    - q           The q of the qp-norm
    - pdType      Primal dual type. ['iso-iso', 'aniso-aniso', 'iso-aniso',
                            eit.de           'aniso-iso']
    - pdAlpha     Regularisation parameter of the primal dual algorithm
    - pdIter      Iterations of the primal dual algorithm. Default: 3000
    - pdAcceleration
                Acceleration type for the PD algorithm
 # Outputs:
    - result      Resulting data as a Nxd vector


 Other m-files required: PrimalDual, GenerateFlowGraph
 MAT-files required: none

 See also: PrimalDual,  GenerateFlowGraph
"""
function cutpursuit(d::Int64, f0::Array{T, 1}, g::VariationalGraph, w::Array{T, 1}, cutType::String, cutAlpha, nIter::Int64) where {T <: Real}
    N = g.num_verts
    # Initialize f with the mean value in each channel
    f = similar(f0)
    for i = 1:d
        mean = sum(f0[((i - 1) * N + 1):(i* N)])/N
        for j = ((i - 1) * N + 1):(i* N)
            f[j] = mean
        end
    end
    E = get_alt_edgerep(g)
    E_red = similar(E)
    w_red = similar(w)
    ct = 0
    # Initialize bins and binsizes
    bins = ones(N)
    binsize = N
    ## Main Routine
    if cutType == "aniso"
        # Initialize a directed LightGraph
        G = DiGraph(N)
        G.fadjlist = g.edges
        G.ne = g.num_edges

        # Array for cut_inds
        B = falses(N)
        
        iter = 1
        # Dict for saving history
        history = Dict([("f", zeros(N, d, nIter)), ("bins", zeros(N, nIter)), 
                        ("energy_red", zeros(nIter, 1)), ("energy_full", zeros(nIter, 1))])
        # Treat each dimension independently
        for i = 1:d
            # Object to save cutted graph
            G_cut = G
            # Gathers the edges that were cut
            cutedge = falses(g.num_edges)
            
            iter = 1
            
            while iter <= 1 #&& ~isempty(find(cutedge==0,1))
                # Cut graph
                # Go through decoupled problems
                # Generate the corresponding flow-graph
                g_flow, c_matrix, gradJS, Sc = generate_flowgraph(1, f0[((i - 1) * N + 1):(i * N)], 
                                                                  w, f[((i - 1) * N + 1):(i * N)], 
                                                                  g, cutAlpha, "aniso")
                cs, = mincut(g_flow, N + 1, N + 2, c_matrix, BoykovKolmogorovAlgorithm())
                # Put together for new reduced problem
                # Get the binary partition
                display(gradJS)
                display(collect(edges(g_flow)))
                display(cs)
                # TR: This could possibly be done better
                ctt = 1
                energy_cut = 0
                for j = 1:N
                    if(ctt <= length(cs)) && (j == cs[ctt])
                        B[j] = true
                        ctt += 1
                    else
                        B[j] = false
                        energy_cut += gradJS[j]
                    end
                end
                
                # flag for convergence and value for weight sum
                conv = true
                w_sum = 0
                # Iteration over all edges
                for j = 1:g.num_verts
                    if !cutedge[j] && xor(B[E[1, j]], B[E[2, j]]) 
                        conv = false
                        cutedge[j] = true
                        # remove this edge from the graph
                        rem_edge!(G_cut, E[1, j], E[2, j])
                        if Sc[j]
                            w_sum += w[j]
                        end
                    end
                end
                
                energy_cut += 0.5 * cutAlpha * w_sum
                
                
                display(energy_cut)
                # Convergence?
                if conv
                    display("Algorithm converged, minimal cut found.")
                    break
                end
                
                set_bins!(bins, G_cut)              
                set_red_graph!(E_red, w_red, ct, E[:, cutedge], w[cutedge], bins)
    
                #             % perform debiasing as last step if needed for pdAlpha = 0
                #             if debiasing
                #                 f = PrimalDual(...
                #                     'Data',f0_red,'Edges',E_red,'Weights',w_red,...
                #                     'regType', pdType, 'alpha', 0, 'acceleration',pdAcceleration,...
                #                     'preconConstant',1, 'nIter', pdIter, 'errType', 'gap', 'tol',1e-10, ...
                #                     'p',p,'q',q,...
                #                     'verbose', verbosity,...
                #                     'saveHistory', 'false');
                #                 f = f(bins,:);
                #             end
                #             disp('Algorithm converged. Leaving iteration scheme.');
                #             break;
                #         end
                #         
                #         % Concarnate the edges between the new partitions to the former
                #         E_bet = E(cutedge,:);
                #         
                #         % Remove the edges from the graph
                #         G_cut = G;
                #         G_cut = G_cut.rmedge(E_bet(:,1),E_bet(:,2));
                #         
                #         % Compute the connected components of the cutted graph
                #         bins = conncomp(G_cut);
                #         bins = bins';
                #         
                #         % give some user output if needed
                #         if verbosity >=2
                #             fprintf('\t %s \n', ['Number of bins is ' num2str(max(bins))] );
                #         end
                #         
                #         % Compute the reduced edgeset and weights
                #         [E_red, w_red] = ComputeReducedGraph(E_bet, w0, cutedge, bins);
                #         w(cut_inds) = 0;
                #         % Compute the reduced observation data
                #         f0_red = zeros(max(bins),d);
                #         for i=1:d
                #             f0_red(:,i) = accumarray(bins,f0(:,i))./accumarray(bins,ones(length(f0),1));
                #         end
                #         
                #         pdAlpha_red = length(f0_red)/length(bins)*pdAlpha;
                #         %pdAlpha_red = pdAlpha;
                #         % Compute primal dual on reduced problem
                #         f = PrimalDual(...
                #             'Data',f0_red,'Edges',E_red,'Weights',w_red,...
                #             'regType', pdType, 'alpha', pdAlpha_red, 'acceleration',pdAcceleration,...
                #             'preconConstant',1, 'nIter', pdIter, 'errType', 'gap', 'tol',1e-10, ...
                #             'p',p,'q',q,...
                #             'verbose', verbosity,...
                #             'saveHistory', 'false');
                #         
                #         % Save energy of reduced problem in history
                #         history.energy_red(iter) = ComputeEnergy(f,f0_red,E_red,w_red,pdAlpha_red,pdType,p,q);
                #         
                #         
                #         f = f(bins,:);
                #         history.energy_full(iter) = ComputeEnergy(f,f0,E,w0,pdAlpha,pdType,p,q);
                #         
                #         % show intermediate results if needed
                #         if verbosity >= 2
                #             figure(100);
                #             subplot(1,2,1); plot3(f0(:,1),f0(:,2),f0(:,3), 'r.');
                #             subplot(1,2,2); plot3(f(:,1),f(:,2),f(:,3), 'b.');
                #             pause(0.001);
                #         end
                #         
                #         if verbosity >= 3
                #             figure; scatter3(f0(:,1),f0(:,2),f0(:,3),10,bins,'filled');
                #         end
                #          
                #         % save intermediate results
                #         history.f(:,:,iter) = f;
                #         history.bins(:,iter) = bins;
                #         history.w_red{iter} = w_red;
                #         history.E_red{iter} = E_red;
                #         % increase iteration counter
                #         iter = iter + 1;
                #     end
                #     
                #     % if iteration was stopped before adjust history accordingly
                #     if iter < nIter
                #         history.f(:,:,iter:end) = [];
                #         history.bins(:,iter:end) = [];
                #     end
            end
        end
    elseif strcmp(cutType,"iso")
        error("Currently Unsupported")
    end
end
# % perform debiasing with pdAlpha=0 if needed
# % if debiasing
# %     f = PrimalDual(...
# %             'Data',f0_red,'Edges',E_red,'Weights',w_red,...
# %             'regType', pdType, 'alpha', 0, 'acceleration',pdAcceleration,...
# %             'preconConstant',1, 'nIter', pdIter, 'errType', 'gap', 'tol',1e-10, ...
# %             'p',p,'q',q,...
# %             'verbose', verbosity,...
# %             'saveHistory', 'false');
# % end
# 
# result = f;
# tellem(verbosity,'\nFinished cut pursuit!\n');
# 
# end
# 
# function [E_red, w_red] = ComputeReducedGraph(E_bet, w0, cutedge, bins)
#     [E_red, ~, ic] = unique([bins(E_bet(:,1)),bins(E_bet(:,2))],'rows');
#     w_red = accumarray(ic,w0(cutedge));
# end
# 
# function tellem(verbose, varargin)
#     if ~verbose
#         return;
#     end
#     fprintf(varargin{:}); 
# end
# 
# %% Code to compute energy and stuff
# function energy = ComputeEnergy(f,f0,E,w,alpha,type,p,q)
#     
# % Compute primal
# [D,R] = ComputePrimal(f,f0,E,w,type,p,q);
# energy = D+alpha*R;
# 
# % fprintf('Energy: %d\n Dataterm: %d, Regularization: %d\n',energy, D, R);
# end
# 
# 
# function [D,R] = ComputePrimal(f,f0,E,w,type,p,q)
#     [N,~] = size(f);
#     D  = 1/2*sum((f(:)-f0(:)).^2);
#     gradf = grad(f,E,w);
#     nz = ComputeNorm(gradf,N,E,type,p,q);
#         
#     % Compute the regularizer
#     R = sum(nz);
# end
# 
# function gradf = grad(f,E,w)
#     d = size(f,2);
#     gradf = zeros(length(E),d);
#     
#     for l=1:d
#         % Get rid of sqrt somehow..
#         gradf(:,l) = sqrt(w).*(f(E(:,2),l)-f(E(:,1),l));
#     end
# end
# 
# function normY = ComputeNorm(y,N,E,normType,p,q)
# % This function computes the different norms for different types of
# % minimization problems.
# %
# % Syntax:  normY = ComputeNorm(y,N,E,normType)
# %
# % Inputs:
# %    y          - Input data of the size of a gradient (length(E)xd)
# %    N          - Number of vertices 
# %    E          - Edge set corresponding to N
# %    normType   - Type of the minimization problem. 
# %                   ['aniso-aniso','aniso-iso','iso-aniso','iso-iso'] 
# %   p           - The p of the qp-norm
# %   q           - The q of the qp-norm
# %
# % Outputs:
# %    normY     - Norm of the input y
# %
# 
# if nargin < 5
#     p=2;
# end
# if nargin < 6
#     q=1;
# end
# % Get the dimension of the data
# [~,d] = size(y);
# 
# if p==1 % Every norm is the same as the aniso-aniso norm!
#     normY = abs(y(:)); 
# elseif p==inf
#     switch normType
#         case 'aniso-aniso'
#             % Absolute value of whole y
#             normY = abs(y(:));
#         case 'aniso-iso'
#             % Maximum in dimensions 
#             normY = max(abs(y),[],2);
#         case 'iso-aniso'
#             % Maximum of each set of neighbors
#             normY = accumarray(...
#                 reshape(repmat(E(:,1),1,d)+(0:N:N*(d-1)),[],1),...
#                 abs(y(:)),[d*N,1],@max);
#         case 'iso-iso'
#             % Maximum over each set of neigbors in all dimensions
#             normY = accumarray(repmat(E(:,1),d,1),abs(y(:)),[N,1],@max);  
#     end
# else
#     switch normType
#         case 'aniso-aniso'
#             % Absolute value of each entry
#             normY = abs(y(:));
#         case 'aniso-iso'
#             % lp-norm of each dimension in y
#             normY = (sum(abs(y).^p,2)).^(1/p);
#         case 'iso-aniso'
#             % lp-norm over each set of neighbors
#             normY = (accumarray(...
#                 reshape(repmat(E(:,1),1,d)+(0:N:N*(d-1)),[],1),...
#                 abs(y(:)).^p,[d*N,1])).^(1/p);
#         case 'iso-iso'
#             % lp-norm over array of channels and neighbors
#             normY = zeros(N,1);
#             for l=1:d
#                 normY = normY + accumarray(E(:,1),abs(y(:,l)).^p,[N,1]);
#             end
#             normY = normY.^(1/p);  
#     end
# end
# % When q is infinity the qp-norm is the max over the p norm
# if q == inf
#     normY = max(normY(:));
# elseif q ~=1
#     normY = sum(normY(:).^q).^(1/q);
# end
# % If q is 1 then the qp-norm becomes a sum over inner p-norms
# 
# end
