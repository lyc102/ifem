function matrix = genMatCurl3D(pde,mesh,fem)
%% Generate global matrices and load vector of FEM for 3D ellipic eq
%     -div(A grad u)  = f,    x\in \Omega
% INPUTS:
% pde --- given data function from equation, e.g. 
%         pde.A --- diffusion coefficient
%         pde.f --- right hand side function
%         pde.gD --- Dirichlet boundary value function 
%         pde.one --- constant function 1.
% mesh --- mesh structure. 
% fem --- global degree of freedom of FEM 
%
% OUTPUTS:
% matrix.S --- stiffness matrix (w/o boundary condition)
% matrix.A --- final FEM matrix (after boundary condition)
% matrix.rhsF --- load vector (w/o boundary condition)
% matrix.f --- final RHS matrix (after boundary condition)

% Last Modified: 08/07/2020 by Xu Zhang

%% 1. Stiffness Matrix
S = globMatrixNed3D(pde.A,1,mesh,fem,fem);
M = globMatrixNed3D(pde.B,0,mesh,fem,fem);

%% 2. Generate the Right Hand Side Vector
rhsF1 = globNedRHS3D(pde.f1, mesh, fem, 0, 1);
rhsF2 = globNedRHS3D(pde.f2, mesh, fem, 0, 2);
rhsF3 = globNedRHS3D(pde.f3, mesh, fem, 0, 3);
rhsF = rhsF1 + rhsF2 + rhsF3;

%% 3. Dirichlet Boundary Conditions
Atotal = S + M;
tu1 = sum(feval(pde.exactu1,fem.gex,fem.gey,fem.gez).*fem.gew,2);
tu2 = sum(feval(pde.exactu2,fem.gex,fem.gey,fem.gez).*fem.gew,2);
tu3 = sum(feval(pde.exactu3,fem.gex,fem.gey,fem.gez).*fem.gew,2);
tgt = mesh.p(mesh.e(:,2),:) - mesh.p(mesh.e(:,1),:);
tgt = tgt./sum(tgt.^2,2).^(1/2);
tu = tu1.*tgt(:,1) + tu2.*tgt(:,2) + tu3.*tgt(:,3);
ub = tu;
ub(fem.mapper) = 0;
rhsB = Atotal*ub;
A = Atotal(fem.mapper,fem.mapper);
f = rhsF(fem.mapper) - rhsB(fem.mapper);

%% Outputs
matrix = struct('A', A, 'f', f, 'S', S, 'rhsF',rhsF,'tu',tu);