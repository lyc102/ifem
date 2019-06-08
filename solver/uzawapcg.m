function [sigma,u,info] = uzawapcg(M,B,C,f,g,elem,option)
%% UZAWAPCG solves saddle point problem discretized from mixed FEM.
% 
% [u,sigma] = uzawapcg(M,B,C,f,g,elem) solves saddle point problem discretized
% from mixed FEM.
%
%      |M  B'| |sigma|  = |f|
%      |B -C | |u|      = |g|
%
%
% See also: dmg
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

t = cputime;
%% Parameters
if ~exist('option','var'), option = []; end
Ndof = length(f)+length(g); 
option = mgoptions(option,Ndof);    % parameters
x0 = option.x0; tol = option.tol; maxIt = option.solvermaxit; 
printlevel = option.printlevel; 
% additional parameters for interiori pcg 
tolM = 1e-8; pcgmaxIt = 200; 

%% Set up auxiliary matrices
N = size(M,1); Ng = length(g);
Bt = B';
DM = spdiags(diag(M),0,N,N);
if isempty(C)
   C = sparse(Ng,Ng);
end
Abar = B*spdiags(1./diag(M),0,N,N)*Bt + C;

%% Set up matrices in each level
setupOption.solver = 'NO';
[x,info,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(Abar,g,elem,setupOption); %#ok<ASGLU>

%% PCG iteration for solving the Schur complement equation
%
%   (B*M^{-1}*B' + C)u = B*M^{-1}f - g
%
% Preconditioner: B*DM^{-1}*B' + C which is an elliptic matrix on the dual
% grid and can be solved by amg or one Vcycle.

b = B*(1./diag(M).*f) - g; % approximated rhs
nb = norm(b);
% initial residual
u = x0(N+1:end);
[tempr,flag] = pcg(M,f-Bt*u,tolM,pcgmaxIt,DM); % M^{-1}(B'*u)
r = B*tempr - g - C*u;

err = zeros(maxIt,2);
err(1,:) = norm(r)/nb; 
k = 1;
% options for V-cycle of Schur complement
% options for V-cycle of Schur complement
if strcmp(option.solver,'CG') % change default set up in mgoption 
   option.solver = 'Vcycle';
end
if isfield(option,'Vit')
   option.solvermaxIt = option.Vit;
else
   option.solvermaxIt = 1; 
end
option.setupflag = false;
option.printlevel = 0;
option.x0 = zeros(Ng,1);
while (max(err(k,:)) > tol) && (k <= maxIt)
    % given r, compute Pr    
%     Pr = Abar\r;  
%     Pr = amg(Abar,r,option);
    Pr = mg(Abar,r,elem,option,Ai,Bi,BBi,Res,Pro,isFreeDof);
    % update tau, beta, and p
    rho = Pr'*r;  % e'*ABA*e approximates e'*A*e
    if k == 1
        p = Pr;
    else
        beta = rho/rho_old;
        p = Pr + beta*p;
    end
    % update alpha, u, and r
    [tempp,flag] = pcg(M,Bt*p,tolM,pcgmaxIt,DM); %#ok<*NASGU> % M^{-1}(B'*p)
    Ap = B*tempp + C*p;  % A*p = B*M^{-1}*B'*p + C*p; 
    alpha = rho/(Ap'*p);
    r = r - alpha*Ap;
    u = u + alpha*p;
    rho_old = rho;
    % compute err for the stopping criterion
    k = k + 1;
    err(k,1) = sqrt(abs(rho/(u'*b))); % approximate relative error in energy norm
    err(k,2) = norm(r)/nb; % relative error of the residual in L2-norm
end
stopErr = max(err(k,:));
itStep = k-1;
if k> maxIt
    flag = 1;
else
    flag = 0;
end

%% Solve sigma     
[sigma,flag] = pcg(M,f-Bt*u,tol,pcgmaxIt,DM);
time = cputime - t;
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',stopErr);

%% Output
if printlevel >= 1
    fprintf('Uzawa-type MultiGrid Preconditioned PCG \n');
    fprintf('#dof: %8.0u,  #nnz: %8.0u, V-cycle: %2.0u, iter: %2.0u,   err = %8.2e,   time = %4.2g s\n',...
                 Ndof, nnz(M)+2*nnz(B), 1, itStep, stopErr, time)
end
if (flag == 1) && (printlevel>0)
   fprintf('NOTE: the iterative method does not converge! \n');    
end
