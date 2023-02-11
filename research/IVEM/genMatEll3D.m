function matrix = genMatEll3D(pde,mesh,fem)
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
S = globMatrix3DStiff(pde.A,mesh,fem,fem);

%% 2. Generate the Right Hand Side Vector
rhsF = globRHS3D(pde.f, mesh, fem, [0,0,0]);

%% 3. Dirichlet Boundary Conditions
Atotal = S;
tu = feval(pde.gD,fem.p(:,1),fem.p(:,2),fem.p(:,3));
ub = tu;
ub(fem.mapper) = 0;
rhsB = Atotal*ub;
A = Atotal(fem.mapper,fem.mapper);
f = rhsF(fem.mapper) - rhsB(fem.mapper);

%% Outputs
matrix = struct('A', A, 'f', f, 'S', S, 'rhsF',rhsF);