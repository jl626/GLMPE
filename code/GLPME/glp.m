function [L, S] = glp(X, alpha, beta, G1, G2, param)
% GLPME low rank pursuit part
%
%   Optimise graph regularised low rank pursuit 
%
% Objective function : 
% (1)  \{\Sigma_{i=1}^2 \Vert M_i - U_i \Vert_1 + 
% (2)  \alpha \Vert \nabla_{G_p} U_i^{T} \Vert_1 + 
% (3)  \beta \Vert \nabla_{G_i} U_i \Vert_F^2 \}
%
%   Inputs: 
%       X       : Input Gene matrix size(p,n)
%       alpha   : Regularization parameter for the shared gene graph
%       beta    : Regularization parameter for the sample graph
%       G1      : Graph of p shared genes size(p,p)
%       G2      : Graph of samples of data i size(n,n)
%       param   : Optional parameters
%
%   Outputs:
%       Lr      : Low-rank part 
%       Sp      : Sparse part 
%
%  This is standalone version of Graph regularised low rank pursuit
%  The original code of Shahid's algorithm is available at EPFL and Shahid's personal website. 
%  https://sites.google.com/site/naumanshahid6/research-teaching
%
%
%   References: Shahid et al. (2016) PCA using graph total variation
%   (ICASSP)
%                            Juheon Lee
%                          26 April 2017

%% initialisation

if nargin>5
    param = struct;
end

param.verbose = 1;
param.paramnn = struct;

% maximum eigenvalue
    G1.lmax = compute_lmax(G1.L);
    G2.lmax = compute_lmax(G2.L);
% graph gradient matrix
    G1.Diff = compute_gradient(G1.W);
    G2.Diff = compute_gradient(G2.W);

paraml1.verbose = 2;
paraml1.y = X;

% For data term (1) in the objective function
f1.prox = @(x,T) prox(x,T,paraml1);
f1.eval = @(x) norm(x(:),1);

% Graph regularisation term 1 (Shared gene graph)
D1 = G1.Diff;
param_tvl1.verbose = 1;
f2.prox = @(x,T) prox(x,T*alpha,param_tvl1);
f2.eval = @(x) alpha*sum(sum(abs(D*x)));
f2.L= @(x) D1*x;
f2.Lt = @(x) D1'*x;
f2.norm_L = 2*G1.lmax;

% Graph regularisation term 2 (Sample graph)
f3.grad = @(x) beta*(2*x*G2.L);
f3.eval = @(x) beta*sum(norm_tik(G2,x'));
f3.beta = 2*beta*G2.lmax;

paramsolver = param;
paramsolver.algo = 'fb_based_primal_dual';

% Plug into optimisation scheme (forward-backward based primal dual method)
L = fbpd(X,{f1,f2,f3},paramsolver); % Low rank term

S = X-L; % Sparse term (see details for Candes et al. (2011) Robust Principal Component Analysis?)

end


