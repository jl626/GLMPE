%% TCGA

%data01 = readtable('TCGA_BRCA_M.csv','TreatAsEmpty',{'.','NA','N/A'});
%data02 = readtable('TCGA_BRCA_R.csv','TreatAsEmpty',{'.','NA','N/A'});
%data1 = data1';
%data2 = data2';
data01 = readtable('TCGA_BRCA_M.csv','TreatAsEmpty',{'.','NA','N/A'});
data02 = readtable('TCGA_BRCA_R.csv','TreatAsEmpty',{'.','NA','N/A'});
data1 = table2array(data01(:,2:end));
data2 = table2array(data02(:,2:end));
data3 = log2(data2+1);
data1(isnan(data1)) = 0;
data3(isnan(data3)) = 0;
data01 = data1;
data02 = data3;
data1 = data01;
data2= data02;
% filter gene
%idx = find(mean(d2,1)>0);%0.1);

%%
param_data.N = size(data1,1);
param_data.Nx = size(data1,2);
param_data.Ny = 1;

%param_graph.use_flann = 0; % to use the fast C++ library for construction of huge graphs.
%param_graph.use_l1 = 0; 
%param_graph.k = 20;
%c = data1;
%param_data.X = c;
%G1 = gsp_nn_graph(param_data.X,param_graph);
param_graph.use_flann = 0; % to use the fast C++ library for construction of huge graphs.
param_graph.use_l1 = 0; 
param_graph.k = 20;
c = data2;
param_data.X = c;
G3 = gsp_nn_graph(param_data.X,param_graph);
%{

G5 = G3;
G5.W = (G1.W + G3.W)/2;
G5.d = sum(G5.W)';
G5.A = (G3.A + G1.A) > 0;
G5.L = (diag(G5.d)-G5.W);
G5.sigma = (G1.sigma + G3.sigma)/2;
%}
param_graph.k = 5;
% build sample graph
param_data.X = data1';
G2 = gsp_nn_graph(param_data.X,param_graph);
param_data.X = data2';
G4 = gsp_nn_graph(param_data.X,param_graph);
%%
param.maxit = 250;
param.tol = 1e-6;
%for iii = 1:6
% for jjj = 1:6
gamma1 = 0.1;%param1(iii);  % joint graph
gamma2 = 0;%param2(jjj);  % sample distribution
param_data.X = data1;
n1 = gsp_gpcatv_2g(param_data.X, gamma1, gamma2, G3, G2,param)';
%n1 = gsp_frpcaog_2g(param_data.X, gamma1, gamma2, G3, G4,param)';
%gamma1 = 0.3;  % joint graph
%gamma2 = 0.3;  % sample distribution
param_data.X = data2;
n2 = gsp_gpcatv_2g(param_data.X, gamma1, gamma2, G3, G4,param)';
%n2 = gsp_frpcaog_2g(param_data.X, gamma1, gamma2, G3, G4,param)';
%save(sprintf('example1_%03d_%03d.mat',round(param1(iii)*100),round(param2(jjj)*100)),'n1','n2');%,'nn1','nn2');
% end
%end
% applying optimal transport
option.bin = 10;
%option.bound01 = [min(n1(:)) max(n1(:))];
%option.bound02 = [min(n2(:)) max(n1(:))];
new1 = zeros(size(n1));
new2 = zeros(size(n1));
for i = 1:length(n1); [u1,u2] = midway_eq(n1(:,i),n2(:,i),option);new1(:,i) = u1;new2(:,i) = u2;end
%new11 = zeros(size(n1));
%new22 = zeros(size(n1));
%for i = 1:length(n1); [u1,u2] = midway_eq(data1(:,i),data2(:,i),option);new1(:,i) = u1;new2(:,i) = u2;end
%figure
%plot(data1(:),data2(:),'.');
%figure
%plot(n1(:),n2(:),'.');
%figure
%plot(new1(:),new2(:),'.');
%figure
%plot(z1(:),z2(:),'.');
