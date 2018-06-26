function y = norm_tik(G,x)
%GSP_NORM_TIK Squared L2 norm of the gradient on graph
%   Usage:  y = gsp_norm_tv(G,x);
%
%   Input parameters:
%         G     : Graph structure (or symetric positive matrix)
%         x     : Signal on graph
%   Output parameters:
%         y     : Norm
%
%   Compute the squared L2 norm of the gradient on graph. If x is a matrix
%   a vector of norm is returned.
%
%   This function can also be used for general symetric positive matrices


if isa(x,'single')
   x = double(x); 
end

if ~isnumeric(G)
    L = G.L;
else
    L = G;
end

 y = sum(x .* (L* x) );

end

