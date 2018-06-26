function compute_acc(z1,z2,label)

K = 10;
alpha = 0.3; 
C = 9;

idx1 = find(std(z1,[],2)>1);
Data1 = zscore(z1(idx1,:)');
Dist1 = pdist2(Data1,Data1);

W2 = affinityMatrix(Dist1, K, alpha);
group1 = SpectralClustering(W2,C);

idx2 = find(std(z2,[],2)>1);
Data1 = zscore(z2(idx2,:)');
Dist1 = pdist2(Data1,Data1);
W2 = affinityMatrix(Dist1, K, alpha);
group2 = SpectralClustering(W2,C);
sprintf('\n cluster accuracy 1: %.2f  \\ cluster accuracy 2: %.2f \\ cluster consistency %.2f ',Cal_NMI(group1,label),Cal_NMI(group2,label),Cal_NMI(group1,group2))
