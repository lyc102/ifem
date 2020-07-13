function [u,p,info] =  diagpreStokes(A,B,C,f,g,elem,option,varargin)
%% diagpreStokes: solve the following saddle point system 
%
% [u,relres,itStep] =  diagpreStokes(A,B,C,f,g,node,elem,bdFlag,option)
% solves saddle point system discretized from mixed FEM
%
%      |A   B'| |u|  = |f|
%      |B  -C | |p|  = |g|
%
% use MINRES method with diagonal preconditioner (M^(-1), A^(-1)) where A =
% B*invM*Bt + C is the Schur complement and its inverse is approximated by
% few V-cycles.
%


t = cputime;
d = size(elem,2) - 1;
% NT = size(elem,1);
Nu = size(A,1);
Nu1 = Nu/d;
Np = size(B,1);
Ndof  = Nu + Np;

%% Form the big matrix
Bt = B';
bigA = [A, Bt; B, C];
bigb = [f; g];

%% Set up the multilevel preconditioner for A
setupOption.solver = 'NO';
% setupOption.freeDof = option.freeDof(1:Nu1);
A1 = A(1:Nu1,1:Nu1);  % get one Poisson
HB = varargin{end};
if isempty(HB)
    [~,~,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(A1,f(1:Nu1),elem,setupOption);
else
    [~,~,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(A1,f(1:Nu1),elem,setupOption,HB);
end

%% Set up diagonal precondition for Schur complement
DAinv = spdiags(1.0./diag(A),0,Nu,Nu); 
S = diag(B*DAinv*Bt + C); % S is a vector

%% Solve by MINRES with diagonal preconditioner
printlevel = option.printlevel;
if isfield(option,'Vit')  % number of Vycles
   option.solvermaxIt = option.Vit;
else
   option.solvermaxIt = 3; 
end
tol = option.tol;
[up,flag,stopErr,itStep] = minres(bigA,bigb,tol,200,@diagpreconditioner);
u = up(1:Nu);
p = up(Nu+1:end);

%% Output
time = cputime - t;
if printlevel >= 1
    fprintf('Diagonal Preconditioner Preconditioned MINRES \n');
    fprintf('#dof: %8.0u,  #nnz: %8.0u, V-cycle: %2.0u, iter: %2.0u,   err = %8.2e,   time = %4.2g s\n',...
                 Ndof, nnz(bigA), option.solvermaxIt, itStep, stopErr, time)
end
if (flag == 1) && (printlevel>0)
   fprintf('NOTE: the iterative method does not converge! \n');    
end
info = struct('solverTime',time,'itStep',itStep,'flag',flag,'stopErr',stopErr);

%% Preconditioner 
 function z = diagpreconditioner(r)
    ru = reshape(r(1:Nu),Nu1,d);
    rp = r(Nu+1:end);
    option.x0 = zeros(Nu1,d);
    option.solver = 'VCYCLE';
    option.printlevel = 0;
    e1 = mg(A1,ru,elem,option,Ai,Bi,BBi,Res,Pro,isFreeDof);
    e2 = rp./S;
    z  = [e1(:); e2];
 end
end