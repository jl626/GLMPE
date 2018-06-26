d1 = tsne([zscore(new1(:,idx));zscore(new2(:,idx))],[label;label]);
d2 = d1(123:end,:);figure;
list = 1:length(d2);
class_plot(d1(1:122,:),label,col,'.',10);hold on;
text(d1(1:122,1),d1(1:122,2),num2str(list'))
class_plot(d2,label,col,'x',6);hold on;
text(d2(1:122,1),d2(1:122,2),num2str(list'),'Color','red')

%%
data1 = [M2(1:30,:);M2(62:91,:)];
data2 = [R2(31:60,:);R2(92:121,:)];
[z1,z2] = XPN(data1',data2');
%%
data3 = [M2(31:60,:);M2(92:121,:)];
data4 = [R2(1:30,:);R2(62:91,:)];
label1 = ones(60,1); label1(31:60) = 2;
d1 = tsne([zscore(data4(:,idx));zscore(new1(:,idx))],[label1;label1]);
d2 = d1(61:120,:);figure;
list = 1:length(d2);
class_plot(d1(1:60,:),label1,col,'.',10);hold on;
text(d1(1:length(list),1),d1(1:length(list),2),num2str(list'))
class_plot(d2,label1,col,'x',6);hold on;
text(d2(1:length(list),1),d2(1:length(list),2),num2str(list'),'Color','red')