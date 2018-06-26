function [u_m1,u_m2] =  midway_eq(u_1,u_2,option)
%
%       Midway Image equalisation (1-D optimal transport)
%
%
% input
%       u_1 - data 01
%       u_2 - data 02
%       option - histogram bin size (option.bin) and data range (option.range)
% output
%
%   citation. Delon (2004), Midway image equalization 
%
%
%           Juheon Lee      03/11/2016
%
if ~isfield(option,'bound01'); option.bound01 = [min(u_1(:)),max(u_1(:))];end
if ~isfield(option,'bound02'); option.bound02 = [min(u_2(:)),max(u_2(:))];end

option.bin = 50;
%[~,~,initial] = histcounts(u_1(:),'NumBins',option.bin,'Normalization','cdf','BinLimits',option.bound01);
%option.bin = length(unique(initial))+1;



% create normalsed output matrices
u_m1 = zeros ( size ( u_1 ));
u_m2 = zeros ( size ( u_2 ));

% create normalised cumulative histograms for u_1 and u_2
[n1,v1,loc1] = histcounts(u_1(:),'NumBins',option.bin,'Normalization','cdf','BinLimits',option.bound01);
[n2,v2,loc2] = histcounts(u_2(:),'NumBins',option.bin,'Normalization','cdf','BinLimits',option.bound02);


if max(length(unique(n1)),length(unique(n2))) <70
option.bin = 10;
[n1,v1,loc1] = histcounts(u_1(:),'NumBins',option.bin,'Normalization','cdf','BinLimits',option.bound01);
[n2,v2,loc2] = histcounts(u_2(:),'NumBins',option.bin,'Normalization','cdf','BinLimits',option.bound02);
end

% to remove numerical precision errors
n1(n1>0.9999) = 1;
n2(n2>0.9999) = 1;

m1 = zeros(option.bin,1);
m2 = zeros(option.bin,1);

% compute mean values of each bin
for ii = 1:option.bin
    idx1 = loc1==ii;
    if ~isempty(find(idx1, 1))
        m1(ii) = mean(u_1(idx1));
    end
    idx2 = loc2==ii;
    if ~isempty(find(idx2, 1))
        m2(ii) = mean(u_2(idx2));
    end
end
%       matrix f_12
% u1 loc
% u2 loc
% mean

% find transport function f12 for u_1
f12 = zeros(option.bin,1);

for jj = 1:option.bin
    j = 1;
     while n1(jj)>n2(j) 
         j = j+1;
     end
  f12(jj) =  (m1(jj) + m2(j))/2;
end

% find transport function f21 for u_2
f21 = zeros(option.bin,1);

for kk = 1:option.bin
    k = 1;
     while n2(kk)>n1(k) 
         k = k+1;
     end
  f21(kk) =  (m1(k) + m2(kk))/2;
end

u_v1 = m1 - f12;
u_v2 = m2 - f21;

for ii = 1:option.bin
    idx1 = loc1==ii;
    if ~isempty(find(idx1, 1))
        u_m1(idx1) = u_1(idx1) - u_v1(ii);
    end
    idx2 = loc2==ii;
    if ~isempty(find(idx2, 1))
        u_m2(idx2) =u_2(idx2) - u_v2(ii);
    end
end