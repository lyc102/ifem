%% COMPARESOLVERS compare the direct solver with mg and amg
% 
% help getfemmatrix fore more options.
%
% For saddle point systems, the output eqn contains blocks of the system
%
%      |A   B'| |u|  = |f|
%      |B  -C | |p|  = |g|
%
% For Stokes and Darcy, C = 0 and the system is singular containing a
% constant kenel of p.

%% Stokes equation on square uniform grids
mesh.shape = 'square';
mesh.type = 'uniform';
mesh.size = 2e4;
pde = 'Stokes';
fem = 'P2P0';
% get the matrix
[eqn,T] = getfemmatrix(mesh,pde,fem);
ufreeDof = eqn.ufreeDof;
pDof = eqn.pDof;
Nu = length(eqn.f);
Np = length(eqn.g);

%% Direct solver
tic; 
disp('Direct solver')
bigA = [eqn.A, (eqn.B)'; ...
        eqn.B, sparse(Np,Np)];
bigF = [eqn.f; eqn.g];    
bigFreeDof = [ufreeDof; Nu+pDof];
bigu = zeros(Nu+Np,1);
bigu(bigFreeDof) = bigA(bigFreeDof,bigFreeDof)\bigF(bigFreeDof);
u1 = bigu(1:Nu);
p1 = bigu(Nu+1:end);
toc;

%% Multigrid solver
tic; 
disp('Multigrid solver')
option.elemType = fem;
option.solver  = 'WCYCLE';
u = zeros(Nu,1);
p = zeros(Np,1);
[u(ufreeDof),p,info] = mgstokes(eqn.A(ufreeDof,ufreeDof),eqn.B(:,ufreeDof),...
                                eqn.f(ufreeDof),eqn.g,u(ufreeDof),p,...
                                T.elem,ufreeDof,option);         
toc;
disp(norm(u-u1));