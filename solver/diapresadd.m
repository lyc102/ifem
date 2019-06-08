function [sigmau,relres,itStep] =  diapresadd(M,B,C,f,g,node,elem,bdFlag,option)
%% DIAPRESADD: solve the following saddle point system 
%
% [sigmau,relres,itStep] =  diapresadd(M,B,C,f,g,node,elem,bdFlag,option)
% solves saddle point system discretized from mixed FEM
%
%      |M   B'| |sigma|  = |f|
%      |B  -C | |u|      = |g|
%
% use MINRES method with diagonal preconditioner (M^(-1), A^(-1)) where A =
% B*invM*Bt + C is the Schur complement and its inverse is approximated by
% few V-cycles.
%
% Created by Jie Zhou on July 30, 2014.
%
% Revised by Long Chen on Feb 21, 2014.


t = cputime;
%%
N1 = size(M,1);
N2 = size(B,1);
Bt = B';
if (size(M,2)==size(M,1))
    M = diag(M);   
    bigA =[-M,Bt;B,C];
else
    bigA =[-spdiags(M,0,N1,N1),Bt;B,C];
end
invMv = spdiags(1.0./M,0,N1,N1); 
Abar = B*invMv*Bt + C;
tol = option.tol;
temp = option.solver;
option.solver = 'NO';
[~,~,Ai,Bi,BBi,Res,Pro] = mgHodgeLapE(Abar,g,node,elem,bdFlag,option);
option.solver = temp;
option.setupflag = false;
[sigmau,~,relres,itStep] = minres(bigA,f,tol,N1,@diagpreconditioner);
Ndof  = N1+N2;
time = cputime - t;
fprintf('#dof: %8.0u,  #nnz: %8.0u,  iter: %2.0u,  err = %8.4e,  time = %4.2g s\n',...
Ndof, nnz(bigA), itStep, relres, time);

%% Preconditioner 
 function z = diagpreconditioner(r)
    r1 = r(1:N1);
    r2 = r(N1+1:end);
    z1 = r1./M;
    option.x0 = zeros(N2,1);
    z2 = mgHodgeLapE(Abar,r2,node,elem,bdFlag,option,Ai,Bi,BBi,Res,Pro);
    z  = [z1; z2];
 end
end