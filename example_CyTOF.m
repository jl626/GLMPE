data1 = dlmread('Person1_3month_day1.csv');
%data1 = log2(data1+1);
data2 = dlmread('person1_baseline_day1.csv');
%data2 = log2(data2+1);
data3 = dlmread('Person2_3month_day1.csv');
%data3 = log2(data3+1);
data4 = dlmread('person2_baseline_day1.csv');
%data4 = log2(data4+1);
%batch1 = [data1;data2;data3;data4];


data5 = dlmread('Person1_3month_day2.csv');
%data5 = log2(data5+1);
data6 = dlmread('person1_baseline_day2.csv');
%data6 = log2(data6+1);
data7 = dlmread('Person2_3month_day2.csv');
%data7 = log2(data7+1);
data8 = dlmread('person2_baseline_day2.csv');
%data8 = log2(data8+1);
%batch2 = [data5;data6;data7;data8];



addpath('C:/Users/Juheon Lee/Documents/normalisation/RPCAG_FRPCAG/utils/');
addpath('C:/Users/Juheon Lee/Documents/normalisation/RPCAG_FRPCAG/gspbox/');
addpath('C:/Users/Juheon Lee/Documents/normalisation/RPCAG_FRPCAG/utils/fast_kmeans/');
addpath('C:/Users/Juheon Lee/Documents/normalisation/RPCAG_FRPCAG/orl_faces/');
addpath('C:/Users/Juheon Lee/Documents/normalisation/RPCAG_FRPCAG/algorithms/');
gsp_start;
init_unlocbox;
%%
tic;
Data1 = data1;
Data2 = data5;
% batch correction between Data1 and Data 5
param_data1.N = size(Data1,1);
param_data1.Nx = 25;
param_data1.Ny = 1;
param_graph.k = 5;
param_data1.X = Data1;
G2 = gsp_nn_graph(param_data1.X,param_graph);

% batch correction between Data1 and Data 5
param_data2.N = size(Data2,1);
param_data2.Nx = 25;
param_data2.Ny = 1;
param_graph.k = 5;
param_data2.X = Data2;
G3 = gsp_nn_graph(param_data2.X,param_graph);

% common gene
param_data.N = size(Data1,1)+size(Data2,1);
param_data.Nx = 25;
param_graph.k = 5;
c = [Data2]';
param_data.X = c;
G1 = gsp_nn_graph(param_data.X,param_graph);

param_data1.X = Data1';
[n1] = gsp_gpcatv_2g(param_data1.X, 0.1, 0.1, G1, G2)';
param_data2.X = Data2';
[n2] = gsp_gpcatv_2g(param_data2.X,0.1, 0.1, G1, G3)';
%
m1 = zeros(size(n1));
m2 = zeros(size(n2));
for i = 1:size(n1,2); [u1,u2] = midway_eq(n1(:,i),n2(:,i));m1(:,i) = u1;m2(:,i) = u2;end

toc;
%k1 = zeros(size(Data1));
%k2 = zeros(size(Data2));
%for i = 1:size(n1,2); [u1,u2] = midway_eq(Data1(:,i),Data2(:,i));k1(:,i) = u1;k2(:,i) = u2;end

figure;
plot(Data1(:,14),Data1(:,20),'b.');xlim([-3 3]);ylim([-3 3])
hold on
plot(Data2(:,14),Data2(:,20),'r.')

figure;
plot(m1(:,14),m1(:,20),'b.');xlim([-3 3]);ylim([-3 3])
hold on
plot(m2(:,14),m2(:,20),'r.')


%figure;
%plot(k1(:,14),k1(:,20),'b.')
%hold on
%plot(k2(:,14),k2(:,20),'r.')
