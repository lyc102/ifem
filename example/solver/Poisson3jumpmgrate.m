%% RATE OF CONVERGENCE OF LINEAR ELEMENT FOR POISSON EQUATION
%
% The purpose of this subroutine is to test the performance of multilevel 
% preconditioners for solving the linear system of algebraic equations
% arising from linear finite element discretization of the elliptic partial
% differential equation with jump coefficients.
%
% $-\nabla \cdot (\omega\nabla u) = f $  in  $\Omega=(-1,1)^3$ 
% $u = 1$ on $x==1$ and $u=0$ on $x==-1$, 
% $\omega\nabla u \cdot n = 0 on other faces.
%
% The diffusion coefficent $\omega$ is piecewise constant with large jump.
%
% jumpmgdata1
%      * $\omega(x) = 1/\epsilon$ if $x\in (0, 1)^3$ and 
%      * $\omega = 1$ otherwise.  
%
% jumpmgdata2
%      * $\omega(x) = 1$ if $x\in (-0.5, 0)^3$ or $x\in (0,0.5)^3$ and 
%      * $\omega = \epsilon$ otherwise.  
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

clear variables

%% Setting
% mesh
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1);
bdFlag = setboundary3(node,elem,'Dirichlet','(x==1) | (x==-1)');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
% option
option.L0 = 3;
option.maxIt = 4;
option.elemType = 'P1';
option.printlevel = 1;
option.plotflag = 0;
option.dispflag = 0;
option.rateflag = 0;
% pde
pde = jumpmgdata2;
global epsilon

%% MGCG (multigrid preconditioned CG) solver
for k = 1:4
    epsilon = 10^(-k);
    [err,time,solver,eqn] = femPoisson3(mesh,pde,option);
    fprintf('\n Table: Solver MGCG for epislon = %0.2e \n',epsilon);
    colname = {'#Dof','Steps','Time','Error'};
    disptable(colname,solver.N,[],solver.itStep,'%2.0u',solver.time,'%4.2g',...
                      solver.stopErr,'%0.4e');
end

%% V-cycle solver
for k = 1:4
    epsilon = 10^(-k);
    option.mgoption.solver = 'Vcycle';
    [err,time,solver,eqn] = femPoisson3(mesh,pde,option);
    fprintf('\n Table: Solver V-cycle for epislon = %0.2e \n',epsilon);
    colname = {'#Dof','Steps','Time','Error'};
    disptable(colname,solver.N,[],solver.itStep,'%2.0u',solver.time,'%4.2g',...
                      solver.stopErr,'%0.4e');    
end

%% Conclusion
%
% For 3-D jump coefficients problem, MGCG (multigrid preconditioned CG)
% solver converges uniformly both to the mesh size and the ratio of the
% jump. 
%
%  V-cycle alone doesn't converge for small epsilon and small h.