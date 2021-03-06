%% Single-Cell RNA-Seq clustering
addpath(sprintf('C:/Users/%s/Dropbox/Juheon_work/data_nomalisation/code/code_for_web_upload/Code/GLPME',getenv('username')))


tic;
param_data.N = 4083;
param_data.Nx = 3128;
param_data.Ny = 1;
param_graph.k = 20;

%c = a1;
%param_data.X = c;
%G1 = gsp_nn_graph(param_data.X,param_graph);
c = log2(d22+1);
param_data.X = c;
G3 = compute_weight(param_data.X,param_graph.k);

% build sample graph
param_graph.k = 5;
param_data.X = log2(d11+1);
G2 = compute_weight(param_data.X',param_graph.k);
param_data.X = log2(d22+1);
G4 = compute_weight(param_data.X',param_graph.k);


% Apply low rank pursuit & optimal transport

gamma1 = 0.1; % joint graph
gamma2 = 0.1; % sample distribution

% low rank pursuit
param_data.X = log2(d11+1);
n1 = glp(param_data.X, gamma1, gamma2, G3, G2);
param_data.X = log2(d22+1);
n2 = glp(param_data.X, gamma1, gamma2, G3, G4);

% optimal transport
new1 = zeros(size(n1));new2 = zeros(size(n2));
option.bin = 10;
for i = 1:length(n1); [u1,u2] = midway_eq(n1(i,:),n2(i,:),option);new1(i,:) = u1;new2(i,:) = u2;end

%figure
%plot(new1(:),new2(:),'.');xlim([0 15]); ylim([0 15]); 
toc;
%% evaluation
%compute_acc(a1,b1,label)
%compute_acc(z1,z2,label)
compute_acc(new1,new2,label)
%compute_acc(f1,f2,label)
%compute_acc(e1,e2,label)