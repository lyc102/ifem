function [sigma,u,N,err] = parabolicMixedFEM(dtstr,solver)
%% PARABOLICMIXEDFEM solves the parabolic equation using mixed FEM
%
% parabolicMixedFEM(dtstr,solverstr) discretizes the parabolic equation
% using mixed FEM with step size dt and solves the saddle point equation by
% the solver given by solverstr:
%  'direct': direct solver
%  'uzawapcg': PCG solving the Schur complement equation
%  'dmg': multigrid based on a distributive form
%
% Example
%    parabolicMixedFEM('h^2','uzawapcg');
%    parabolicMixedFEM('1','uzawapcg');
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

[node,elem] = squaremesh([0,1,0,1],0.25);
for k = 1:2
    [node,elem] = uniformrefine(node,elem);
end
bdFlag = setboundary(node,elem,'Dirichlet','all','Neumann','y==1');
pde = mixBCdata;

%% Finite Element Method
maxIt = 4; 
err = zeros(maxIt,2); 
N = zeros(maxIt,1);
option.solver = 'none'; % get the matrix
for k = 1:maxIt
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    h = 1/sqrt(size(elem,1)); %#ok<NASGU>
    dt = eval(dtstr);
    [tempvar1,tempvar2,eqn] = PoissonRT0(node,elem,bdFlag,pde,option);
    v = simplexvolume(node,elem);
    C = spdiags(v/dt,0,size(elem,1),size(elem,1));
    switch lower(solver)
        case 'direct'
            x = eqn.A\[eqn.f;eqn.g];
            sigma = x(1:NE);
            u = x(NE+1:end);
        case 'uzawapcg'
            [sigma,u] = uzawapcg(eqn.M,eqn.B,C,eqn.f,eqn.g,elem);
        case 'dmg'
            [sigma,u] = dmg(eqn.M,eqn.B,C,eqn.f,eqn.g,elem);
    end
    N(k) = size(u,1);
end