function [sol,info] = prox(x, gamma, param)
% Proximal operator
%  Objective Function:
%  argmin_{z} 0.5*||x - z||_2^2 + gamma * ||A z - y ||_1
%
%   Inputs:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Outputs:
%         sol   : Solution.
%         info : Structure summarizing informations at convergence
%
% A modified version from unlocbox package
% http://update-seputar-software.blogspot.com/2014/09/unlocbox-14149-windows.html

% Optional input arguments
if nargin<3, param=struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'At'), param.At = @(x) x; end
if ~isfield(param, 'A'), param.A = @(x) x; end
if ~isfield(param, 'weights'), param.weights = 1; end

% param.y is intitialized below to perform less computations
% test the parameters
if check_gamma(gamma)
    sol = x;
    info.algo=mfilename;
    info.iter=0;
    info.final_eval=0;
    info.crit='--';
    info.time=toc(t1);
    return; 
end

param.weights = 1;

% Projection
    temp = param.A(x);
    if ~isfield(param, 'y'), param.y = zeros(size(temp)); end
    temp2 = param.y + soft_threshold(temp -  param.y , ...
                        gamma*param.nu*param.weights) - temp;
    sol = x + 1/param.nu * param.At(temp2);
    crit = 'REL_OBJ'; iter = 1;
    dummy = temp2 +temp;
    norm_l1 = gamma*sum(param.weights(:).*abs(dummy(:)));

% Log after the prox l1
if param.verbose >= 1
    fprintf(['  prox_L1: ||A x-y||_1 = %e,', ...
        ' %s, iter = %i\n'], norm_l1, crit, iter);
end

end

% soft shrinkage thresholding
function sz = soft_threshold(x,T)
    if T>0
        size_x=size(x);
        x=x(:);
        T=T(:);
        sz = max(abs(x)-T,0)./(max(abs(x)-T,0)+T).*x;
        sz=reshape(sz,size_x);
    else
        sz = x;
    end
end

% check non-negative gamma value
function stop = check_gamma(gamma)

    if gamma<0
        error('gamma can not be negativ!');
    end  
    if gamma==0
        stop = 1;
    else
        stop = 0;
    end
end