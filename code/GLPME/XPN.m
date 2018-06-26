function [z1,z2] = XPN(x1,x2,geneClusters,sampleClusters,nRep)

% Example: 
% [x1new, x2new] = XPN(x1, x2);
%
% XPN normalization of data sets x1, x2. In the k-means clustering, it
% uses sampleClusters array clusters and geneClusters row clusters. In the matrices x1 and x2,
% rows are genes and columns are arrays. The number of genes in x1 and x2
% must be equal. The output is average of nRep iterations of the algorithm.




%% size check
[p_ n1] = size(x1);
[p  n2] = size(x2);
if ( p ~= p_);
	error(strcat('Number of genes in x1 and x2 does not match'));
end;
clear p_;

%% Check for NaNs

if(any(isnan(x1(:))) || any(isnan(x2(:))))
	error('Missing values in input data');
end;

%% Check for geneClusters argument

if(nargin<=2)
	geneClusters = min(25, max(floor(p/50),1));
	disp(['geneClusters=' mat2str(geneClusters)]);
end;

if(nargin<=3)
	sampleClusters = min(5, floor(min(n1,n2)/4));
	disp(['sampleClusters=' mat2str(sampleClusters)]);
end;

if(nargin<=4)
	nRep = 30;
	disp(['nRep=' mat2str(nRep)]);
end;

if(sampleClusters <= 0)
	error('geneClusters is too small');
end;

if(sampleClusters <= 0)
	error('sampleClusters is too small');
end;

if(nRep <= 0)
	error('nRep is too small');
end;

if(nRep < 15)
	warning('XPN:nRep','nRep is small, should be at least 15');
end;


%% Check for constant rows

constrows = ( std([x1 x2],[],2)==0 );
if(any(constrows))
	if(all(constrows))
		error('Constant input data');
	end;
	
	x1fix = x1(~constrows,:);
	x2fix = x2(~constrows,:);
	
	[y1,y2] = XPN(x1fix,x2fix,geneClusters,sampleClusters,nRep);
	
	z1 = zeros(size(x1));
	z1(~constrows,:) = y1;
	z1(constrows,:) = x1(constrows,:);	
	
	z2 = zeros(size(x2));
	z2(~constrows,:) = y2;
	z2(constrows,:) = x2(constrows,:);	
	return;
end;
clear constrows;	
%% Check Skewness

z = x1(:);
m = mean(z);
s = std(z);
skewness1 = mean((z - m).^3) / s^3;

z = x2(:);
m = mean(z);
s = std(z);
skewness2 = mean((z - m).^3) / s^3;

if( (abs(skewness1)>8) || (abs(skewness2)>8) )
	warning('XPN:DataSkewness','Input data is excessively skewed.');
	warning('XPN:DataSkewness','The data might require log-transformation.');
end;

clear z m s skewness1 skewness2;

%% transform platforms to have ZERO median
x1med=median(x1,2);
x2med=median(x2,2);

x1primeMed = (x1 - repmat(x1med,1,n1));
x2primeMed = (x2 - repmat(x2med,1,n2));
xxprimeMed = [x1primeMed x2primeMed];

clear x1med x2med x1primeMed x2primeMed;

TheSum1 = zeros(size(x1));
TheSum2 = zeros(size(x2));
%% Begin cycle
for cycle = 1:nRep

%% kmeans clustering
  % correlation vs sqEuclidean
kmeans_display = 'off'; %'iter' 
  IDR = kmeans(xxprimeMed,geneClusters,'distance','correlation','emptyaction','singleton','rep',1,'maxiter',1000,'display',kmeans_display);
  repeat=true;
  while(repeat)
    IDC = kmeans(xxprimeMed',sampleClusters,'distance','correlation','emptyaction','singleton','rep',1,'maxiter',1000,'display',kmeans_display)';
    IDC1 = IDC(1:n1);
    IDC2 = IDC((n1+1):end);
    countL1=accumarray(IDC1',1,[sampleClusters 1])';
    countL2=accumarray(IDC2',1,[sampleClusters 1])';
    if(min([countL1 countL2])>0)
        repeat=false;
    else
        repeat=true;
        disp('A cluster with samples from only one platform occured. Clustering repeated.');
    end;
  end;
 
%% GLP means and GP variance G - gene, L - colCluster, P - study
  % note of the programmer:
  % x(:) - make a column from a matrix (first col; second col; ...)
  % repmat - repeats a column _c_ (or matrix) a number _nn_ of times: repmat(c,nn,1)
  % kron   - repeats each element of a column _c_ _p_ times:
  % kron(_c_,ones(p,1))
  %
  % accumarray([a b],x(:),) - sums elements into elements (a(i),b(i))
  
  % countL =countL1+countL2; $ not used

  sumGL1  =accumarray([repmat((1:p)',n1,1) kron(IDC1',ones(p,1))],x1(:),[p sampleClusters]);
  sumGL2  =accumarray([repmat((1:p)',n2,1) kron(IDC2',ones(p,1))],x2(:),[p sampleClusters]);
  
  sumGL1sq=accumarray([repmat((1:p)',n1,1) kron(IDC1',ones(p,1))],x1(:).^2,[p sampleClusters]);
  sumGL2sq=accumarray([repmat((1:p)',n2,1) kron(IDC2',ones(p,1))],x2(:).^2,[p sampleClusters]);

  xbar1=(sumGL1./repmat(countL1,p,1));
  xbar2=(sumGL2./repmat(countL2,p,1));
  
  sigmaGL1sq= sumGL1sq ./repmat(countL1,p,1) - (xbar1.^2);
  sigmaGL2sq= sumGL2sq ./repmat(countL2,p,1) - (xbar2.^2);

  % apply MLE to each gene cluster
  
  A1=zeros(geneClusters,sampleClusters);
  A2=zeros(geneClusters,sampleClusters);
  b1=zeros(p,1);
  b2=zeros(p,1);
  c1=zeros(p,1);
  c2=zeros(p,1);
  s1=zeros(p,1);
  s2=zeros(p,1);
  for i=1:geneClusters
      [A1_,b1_,c1_,s1_]=XPN_MLE(xbar1(IDR==i,:),sigmaGL1sq(IDR==i,:),countL1);
      A1(i,:)    = A1_;
      b1(IDR==i) = b1_;
      c1(IDR==i) = c1_;
      s1(IDR==i) = s1_;
      
      [A2_,b2_,c2_,s2_]=XPN_MLE(xbar2(IDR==i,:),sigmaGL2sq(IDR==i,:),countL2);
      A2(i,:)    = A2_;
      b2(IDR==i) = b2_;
      c2(IDR==i) = c2_;
      s2(IDR==i) = s2_;
  end; %i=1:geneClusters
  
  % average the estimates across platforms
  
  A = (A1.*repmat(countL1,geneClusters,1) + A2.*repmat(countL2,geneClusters,1))./(repmat(countL1+countL2,geneClusters,1));
  b = (b1*n1 + b2*n2)/(n1+n2);
  c = (c1*n1 + c2*n2)/(n1+n2);
  s = (s1*n1 + s2*n2)/(n1+n2);
  
%% KsiP residual of the model - P - study
  
  Ksi1 = (x1 - A1(IDR,IDC1).*repmat(b1,1,n1) - repmat(c1,1,n1))./repmat(sqrt(s1),1,n1);
  Ksi2 = (x2 - A2(IDR,IDC2).*repmat(b2,1,n2) - repmat(c2,1,n2))./repmat(sqrt(s2),1,n2);

%% xPstar transformed data - P - study
  
  x1star = A(IDR,IDC1).*repmat(b,1,n1) + repmat(c,1,n1) + Ksi1.*repmat(sqrt(s),1,n1);
  x2star = A(IDR,IDC2).*repmat(b,1,n2) + repmat(c,1,n2) + Ksi2.*repmat(sqrt(s),1,n2);

%% Check

  if(  any(any(isnan(x1star))) || any(any(isnan(x2star))) )
	error('NaN values are produced by the XPN code')
  end;
  
%% Aggregate results of several tries in TheSumP.
  TheSum1=TheSum1+x1star;
  TheSum2=TheSum2+x2star;
  disp([num2str(cycle) ' iterations out of  ' num2str(nRep) ' completed.']);
  
%% end cycle  
end; %cycle = 1:nRep

%% take the average over several tries
z1=TheSum1/nRep;
z2=TheSum2/nRep;

return;
end




function [A,b,c,s2] = XPN_MLE(xbar, sigma2, nj)

	% internal MLE procedure for XPN
	% 
	% xbar matrix of averages of X within column clusters (given gene,
	% platform)
	% sigma2 - variance of the elements within column clusters (given gene,
	% platform)
	% nj - number of elements in the column clusters

	% Remark: the platform/gene cluster is fixed

	% n - number of samples in the platform
	n=sum(nj);

	% I - number of genes in the gene cluster
	% J - number of sample clusters.
	[I J]=size(xbar);

	% A, b, c, s2 - peremeters of the model
	A=zeros(1, J);
	b=ones(I,1);
	c=zeros(I,1);
	s2=ones(I,1);

	% previous values values
	old=[A b' c' s2']*0;

	%i=0;

	while ( sum(([A b' c' s2']-old).^2)>1e-16*max(old) )

		old=[A b' c' s2'];
	%    i=i+1;

		% iteratively update the values
		c = sum( (xbar - repmat(b,1,J).*repmat(A,I,1)).*repmat(nj,I,1)   ,2) / n;

		% fix sign of _b_
		if(sum(b)<0)
			b=-b;
		end;

		A = sum( repmat(b,1,J).*(xbar - repmat(c,1,J))./repmat(s2,1,J)   ,1) / ...
			sum(b.^2./s2);

		% enforce constrains on _A_
		A = A - mean(A);
		A = A * sqrt( J / sum(A.^2));


		b = sum( repmat(A,I,1).*(xbar - repmat(c,1,J)).*repmat(nj,I,1)   ,2) / ...
			sum( A.^2.*nj);


		s2= sum( ...
			((xbar - repmat(c,1,J) - repmat(A,I,1).*repmat(b,1,J)).^2 + sigma2).*repmat(nj,I,1) ...
			,2 ...
				)/n;

		s2(s2==0) = realmin('double');

		% disp(i);disp(sum(([A b' c' s2']-old).^2)); % debug output
	end;
end