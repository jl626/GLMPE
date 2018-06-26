addpath(sprintf('C:/Users/%s/Dropbox/Juheon_work/data_nomalisation/code/code_for_web_upload/Code/GLPME',getenv('username')))

% synthetic data with logistic distribution + gaussian noise
%{
t1 = [1:10;10:-1:1;1:10;10:-1:1];
t2 = [2:2:20;20:-2:2;2:2:20;20:-2:2];
t1 = kron(t1, ones(20,100));
t2 = kron(t2, ones(20,100));
t1_true = t1;
t2_true = power(1.2,t2);
t1_n = t1;
t2_n = t2;
t1 = awgn(t1,1);
t2 = awgn(t2,1);
t2 = power(1.2,t2);
%}
load(sprintf('C:/Users/%s/Dropbox/Juheon_work/data_nomalisation/code/code_for_web_upload/Data/Synthetic/synthetic_data.mat',getenv('username')));

tic;
% build joint gene graph 
param_data.N  = 80;
param_data.Nx = 1000;
k = 20;
c = t1;
param_data.X = c';
G1 = compute_weight(param_data.X,k);

%param_graph.use_flann = 0; % to use the fast C++ library for construction of huge graphs.
%param_graph.use_l1 = 0; 
%param_graph.k = 20;
%c = t2;
%param_data.X = c';
%G3 = gsp_nn_graph(param_data.X,param_graph);

% build sample graph
param_data.N  = 80;
k = 5;
param_data.X = t1;
G2 = compute_weight(param_data.X,k);
param_data.X = t2;
G4 = compute_weight(param_data.X,k);

% optimisation parameters
gamma1 = 0.5; % sample distribution
gamma2 = 0.5; % joint graph

param_data.X = t1';
n1 = glp(param_data.X, gamma1, gamma2, G1, G2)';
%n1 = gsp_frpcaog_2g(param_data.X, gamma1, gamma2, G2, G1)';

param_data.X = t2';
n2 = glp(param_data.X, gamma1, gamma2, G1, G4)';
%n2 = gsp_frpcaog_2g(param_data.X, gamma1, gamma2, G4, G1)';



% conventional method

%[z1,z2] = XPN(t1,log(t2));

new1 = zeros(size(n1));new2 = zeros(size(n2));
option = [];
%new11 = zeros(size(n1));new22 = zeros(size(n2));
for i = 1:size(n1,2); [u1,u2] = midway_eq(n1(:,i),n2(:,i),option);new1(:,i) = u1;new2(:,i) = u2;end
%for i = 1:size(n1,1); [u1,u2] = midway_eq(t1(i,:),t2(i,:),option);new11(i,:) = u1;new22(i,:) = u2;end


figure
plot(t1(:),t2(:),'g.');hold on;xlim([0 15]); ylim([0 70]);
plot(t1_true(:),t2_true(:),'r.','markers',12); 
hold off

figure
plot(new1(:),new2(:),'g.');hold on;xlim([0 40]); ylim([0 70]);
toc;