function [G] = compute_weight(data,radii)
%
%      
%       Input: data= Number_of_data x Number_of_features
%       ouput: W = pair-wise similarity matrix
%
%         25/11/2014

N = size(data,1);
[idx,dist] = knnsearch(data,data,'k',radii+1,'IncludeTies',true); % knn graph
%[idx,dist] = rangesearch(data(:,1:2),data(:,1:2),radii); % e-neighbour graph
d = cell2mat(dist);
sigma = mean(d(:))^2;
    ii = [];
    jj = [];
    cc = [];
for i= 1:N,
        j = idx{i};
        % features
        fdist = dist{i};
        c = exp(-(fdist.^2/sigma)); % kNN weights
        % vectorised allocation
        ii = [ii i*ones(1,numel(j))];
        jj = [jj j];
        cc = [cc c'];
end
W = sparse(ii,jj,cc,N,N);
W = (W+W.')/2; % symmetric knn graph
W(1:(N+1):end) = 0;     % Remove diagonal components;    
D=spdiags(sum(W,2),0,length(W),length(W));% diagonal matrix
G.W = W;
G.d = diag(D);
G.sigma = sigma;
G.type = 'knn';
G.L = D-W;
G.A = logical(W);
G.N = N;
