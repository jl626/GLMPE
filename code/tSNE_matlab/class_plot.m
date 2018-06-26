function class_plot(data,label,col,sym,sz)


cl = max(label);

for i = 1:cl
    hold on;
    idx = find(label==i);
    plot(data(idx,1),data(idx,2),sym,'MarkerSize',sz,'Color',col(i,:)); hold on;
end