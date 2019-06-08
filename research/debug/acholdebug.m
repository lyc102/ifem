N = 400;
A = gallery('poisson', N);
[L2,p,Ac] = achol(A); 
% [L2,p,Ac] = acholold(A); 
% save olddata L2 p Ac
disp(size(Ac,1));
    b = ones(size(A,1),1);
    tol = 1e-6; maxit = 100;    
    tic;
%     Ap = A(p,p);
    [x2,fl2,rr2,it2,rv2] = pcg(A,b,tol,maxit,@(r)acholpre(r,A,L2,L2',p,Ac));
    toc;
    fprintf('#dof: %8.0u,  iter: %2.0u\n',size(A,1), it2)
    semilogy(0:it2,rv2./norm(b),'b.');
    hold on