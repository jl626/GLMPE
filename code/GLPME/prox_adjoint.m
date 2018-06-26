function [sol]  = prox_adjoint(x, gamma, f)

if gamma == 0
    gamma = eps;
end

if nargin < 3,
    error('not enough input arguments');
end

sol_a = gamma * f.prox(x / gamma, 1 / gamma);
sol = x - sol_a;

end


