function [u,p,info] = stokessolver(A,B,f,g,u,p,elem,option) 
%% STOKESSOLVER solvers for Stokes equations
%
% We compute the solution [u,p] of the saddle point problem:
%     [A B'][u] = [f]
%     [B 0 ][p] = [g]
% by various solvers:
%   * 'direct'  the built in direct solver \ (mldivide)
%   * 'uzawa'   inexact Uzawa's method  
%   * 'mg'      multigrid based on LSC-DGS smoother
%   * 'blkdiag' block diagonal preconditioned MINRES
%   * 'blktri'  block triangular preconditioned GMRES
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


ufreeDof = option.ufreeDof;
pDof = option.pDof;

switch option.solver
    case 'direct'
        t = cputime;
        Nu = length(f);
        Np = length(g);
        bigA = [A, B'; ...
                B, sparse(Np,Np)];
        bigF = [f; g];
        bigu = zeros(Nu+Np);
        bigFreeDof = 1:2*Nu+pDof;
        bigu(bigFreeDof) = bigA(bigFreeDof,bigFreeDof)\bigF(bigFreeDof);
        u = bigu(1:2*Nu);
        p = bigu(2*Nu+1:end);
        residual = norm(bigF - bigA*bigu);
        solverTime = cputime - t;
        info = struct('solverTime',solverTime,'itStep',0,'err',residual,'flag',2,'stopErr',residual);        
    case 'mg'
        option.solver  = 'WCYCLE';
        [u,p,info] = mgstokes(A,B,f,g,u,p,elem,ufreeDof,option);         
    case 'asmg'     
        [u,p,info] = asmgstokes(A,B,f,g,u,p,node,elem,option.bdFlag,ufreeDof,option); 
end