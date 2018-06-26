function K = compute_gradient(W)
% Compute graph based gradient matrix
% Originally suggested by Xavier Bresson for graph total variation
% clustering
%
%
%
%
% Reference Bresson et al. (2013) Multi-class total variation clustering, NIPS
%                   Juheon Lee
%                  2017 April 13


n = size(W,2);
triangularW=tril(W); % retain lower triangle parts
[I, J, v]=find(triangularW);
m=length(v); %number of total edges 
GI=[1:m 1:m];
GJ(1:m) = I;
GJ(m+1:2*m) = J;
Gv(1:m)=sqrt(v);
Gv(m+1:2*m)=-sqrt(v);
K=sparse(GI,GJ,Gv,m,n);

