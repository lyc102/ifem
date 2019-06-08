function [sigma,u,info] = tripremixPoisson(M,B,C,f,g,elem,option)
%% TRIPREMIXPOISSON: solve the mixed Poisson system by a triangular preconditioner 
%
% [sigma,u] = tripremixPoisson(M,B,C,f,g,elem) solves saddle point problem
% discretized from mixed FEM for Poisson equation.
% 
%      |M   B'| |sigma|  = |f|
%      |B  -C | |u|      = |g|
%
% Use  |D   B'| as the preconditioner in gmres and compute the inverse by
%      |B  -C | 
%
% the factorization
%
% |D B'| |I Dinv*B'| = |D    0       |
% |B -C| |0      -I|   |B B*Dinv*B'+C|
%
% The elliptic matrix B*Dinv*B'+C on the dual grid is inverted by V-cycle
% or W-cycles MG using mg. Type of cycles is given by option.solver and
% number of cycles is given by option.Vit. The default setting is one
% iteration of multigrid V-cycle.
%
% It is faster than uzawapcg since no need to compute Minv by an inner pcg
% iteration.
%
% See also: uzawapcg
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

t = cputime;
%% Parameters
if ~exist('option','var'), option = []; end
Ndof = length(f)+length(g); 
option = mgoptions(option,Ndof);    % parameters
x0 = option.x0; tol = option.tol; maxIt = option.solvermaxit; 
printlevel = option.printlevel; 

%% Set up auxiliary matrices
N = size(M,1); Ng = length(g);
Bt = B';
% DM = spdiags(diag(M),0,N,N);
if isempty(C)
   C = sparse(Ng,Ng);
end
DMinv = spdiags(1./diag(M),0,N,N);
Abar = B*DMinv*Bt + C;

%% Set up matrices in each level
setupOption.solver = 'NO';
[x,info,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(Abar,g,elem,setupOption); %#ok<ASGLU>

%% Form a big matrix equation
bigA = [M Bt; B -C];
bigF = [f; g];

%% Preconditioned GMRES
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
% gmres for the saddle point system
[x,flag,stopErr,itStep,err] = gmres(bigA,bigF,20,tol,maxIt,@DBC,[],x0);
itStep = (itStep(1)-1)*20 + itStep(2); % total iteration
sigma = x(1:N); 
u = x(N+1:end);

%% Output
time = cputime - t;
if printlevel >= 1
    fprintf('Triangular Preconditioner Preconditioned GMRES \n');
    fprintf('#dof: %8.0u,  #nnz: %8.0u, V-cycle: %2.0u, iter: %2.0u,   err = %8.2e,   time = %4.2g s\n',...
                 Ndof, nnz(bigA), option.solvermaxIt, itStep, stopErr, time)
end
if (flag == 1) && (printlevel>0)
   fprintf('NOTE: the iterative method does not converge! \n');    
end
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',stopErr);

%% Preconditioner
    function s = DBC(r)
        rf = r(1:N); 
        rg = r(N+1:end);
        ds = DMinv*rf;
%         du = Abar\(rg - B*ds);
%         du = amg(Abar,rg - B*ds,option);
        du = mg(Abar,rg - B*ds,elem,option,Ai,Bi,BBi,Res,Pro,isFreeDof);
        s = [ds + DMinv*(Bt*du); -du];
    end
end