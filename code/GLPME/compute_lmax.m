function lmax = compute_lmax(L)
% maximum eigenvalue compute 
    opts=struct('tol',5e-3,'p',min(size(L,1),10),'disp',0);
    lmax=eigs(L,1,'lm',opts);

    lmax=abs(lmax)*1.01; % increase 1% robust against errors
    

